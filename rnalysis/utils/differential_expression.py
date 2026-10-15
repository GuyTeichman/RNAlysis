import abc
import hashlib
import itertools
import re
import time
import unicodedata
from pathlib import Path
from typing import Dict, Iterable, List, Literal, Set, Tuple, Union

import polars as pl

from rnalysis.exceptions import InvalidTypeError, InvalidValueError
from rnalysis.utils import generic, installs, io, parsing, validation


def _fill_template(template: str, values: Dict[str, str]) -> str:
    """
    Replace the given placeholders in an R template in a single pass, \
    so that text inserted for one placeholder (e.g. a design-matrix value that happens to contain \
    '$OUTFILE_NAME') is never itself treated as a placeholder.
    """
    pattern = '|'.join(re.escape(key) for key in sorted(values, key=len, reverse=True))
    return re.sub(pattern, lambda match: values[match.group(0)], template)


class DiffExpRunner(abc.ABC):
    POLY_PATTERN = r'poly\s*\(\s*(\w+)\s*,\s*degree\s*=\s*(\d+)\s*\)'
    TEMPLATES_DIR = Path(__file__).parent.parent.joinpath('data_files/r_templates')
    # prefix of the output file names
    FILENAME_PREFIX: str = NotImplemented
    # how the ':' of an LRT interaction term (e.g. 'factor1:factor2') is spelled in output file names
    LRT_INTERACTION_FILENAME_SEP: str = NotImplemented

    def __init__(
        self,
        data_path: Union[str, Path],
        design_mat_path: Union[str, Path],
        comparisons: Iterable[Tuple[str, str, str]],
        covariates: Iterable[str] = tuple(),
        lrt_factors: Iterable[str] = tuple(),
        model_factors: Iterable[str] = 'auto',
        r_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
    ):
        if not (validation.isiterable(comparisons) and validation.isinstanceiter(comparisons, tuple)):
            raise InvalidTypeError(f'comparisons must be an iterable of tuples, instead got: {comparisons}')
        if not validation.isiterable(covariates):
            raise InvalidTypeError(f'covariates must be an iterable, instead got: {covariates}')
        if not validation.isiterable(lrt_factors):
            raise InvalidTypeError(f'lrt_factors must be an iterable, instead got: {lrt_factors}')
        self.data_path = data_path
        self.design_mat_path = design_mat_path
        self.comparisons = comparisons
        self.covariates = covariates
        self.lrt_factors = lrt_factors
        self.r_installation_folder = r_installation_folder

        self.design_mat: pl.DataFrame = io.load_table(self.design_mat_path)
        if model_factors == 'auto':
            self.model_factors = []
            for factor in itertools.chain(self.design_mat.columns[1:], self.covariates):
                present = factor in self.model_factors
                for lrt in self.lrt_factors:
                    if lrt.startswith(f'poly({factor}'):
                        present = True
                        break
                if not present:
                    self.model_factors.append(factor)
            self.model_factors = self.model_factors + parsing.data_to_list(self.lrt_factors)
        else:
            self.model_factors = model_factors
            for factor in itertools.chain(self.lrt_factors, self.covariates, [comp[0] for comp in self.comparisons]):
                if factor not in self.model_factors:
                    raise InvalidValueError(f'factor {factor} not found in design matrix columns')

    def run(self) -> Path:
        # validate before anything touches R, so a typo fails instantly instead of after a full model fit
        self.validate_design_inputs()
        self.install_required_packages()
        script_path = self.create_script()
        io.run_r_script(script_path, self.r_installation_folder)
        return script_path.parent

    @staticmethod
    def _format_choices(choices: Iterable) -> str:
        return ', '.join(f"'{choice}'" for choice in choices)

    @staticmethod
    def _observed_values(column: pl.Series) -> List:
        return sorted(column.drop_nulls().unique().to_list(), key=str)

    @staticmethod
    def _is_observed_value(value, column: pl.Series) -> bool:
        """
        Check whether a comparison's numerator/denominator is a value of the given design-matrix column. \
        The check is deliberately lenient wherever Python's and R's reading of the design matrix could differ, \
        so that it never rejects an input that R would accept.
        """
        if column.dtype == pl.String:
            # the value is written verbatim into the R script, so it must match exactly (case and whitespace included)
            return str(value) in column.drop_nulls().unique().to_list()
        # other column types (e.g. boolean): leave the judgement to R
        return True

    @staticmethod
    def _find_non_numeric_value(column: pl.Series):
        """
        Return a value of the given design-matrix column that R would not read as a number, \
        or None if R would read the whole column as numeric.
        """
        if column.dtype.is_numeric():
            return None
        if column.dtype != pl.String:
            # e.g. a boolean column, which R reads as logical rather than numeric
            non_null = column.drop_nulls()
            return non_null[0] if len(non_null) > 0 else None
        # R's read.table() reads blank fields of a numeric column as missing values
        values = column.drop_nulls().str.strip_chars()
        values = values.filter(values != '')
        as_numbers = values.cast(pl.Float64, strict=False)
        non_numeric = values.filter(as_numbers.is_null())
        return non_numeric[0] if len(non_numeric) > 0 else None

    def validate_design_inputs(self):
        """
        Validate the comparisons, covariates, and LRT factors against the design matrix, \
        so that a typo fails instantly with a clear message, instead of after a full R model fit.

        :raises InvalidValueError: if a comparison, covariate, or LRT factor does not match the design matrix.
        """
        factors = self.design_mat.columns[1:]
        factors_str = self._format_choices(factors)

        for comparison in self.comparisons:
            if len(comparison) != 3:
                raise InvalidValueError(
                    f'Each comparison must contain exactly three items: '
                    f'(factor, numerator value, denominator value). Instead got: {comparison}'
                )
            factor, numerator, denominator = comparison
            if factor not in factors:
                raise InvalidValueError(
                    f"The factor '{factor}' in the comparison {comparison} is not a column of the design matrix. "
                    f'The available factors are: {factors_str}.'
                )
            column = self.design_mat[factor]
            if column.dtype.is_numeric():
                # R reads an all-numeric column as a continuous variable, which has no levels to compare
                # (DESeq2's results() and limma's contrasts both fail on it)
                raise InvalidValueError(
                    f"The factor '{factor}' in the comparison {comparison} contains only numeric values, "
                    f'so R treats it as a continuous (numeric) variable rather than as groups that can be compared. '
                    f"To compare groups of '{factor}', give its values text labels in the design matrix "
                    f"(for example '{factor}1' and '{factor}2' instead of 1 and 2). "
                    f'To test a continuous variable, add it as a covariate instead.'
                )
            for value in (numerator, denominator):
                if not self._is_observed_value(value, column):
                    raise InvalidValueError(
                        f"'{value}' (in the comparison {comparison}) is not a value of the factor '{factor}' "
                        f"in the design matrix. The values of '{factor}' are: "
                        f'{self._format_choices(self._observed_values(column))}. '
                        f'Note that values must match exactly, including upper/lower case and spaces.'
                    )
            if str(numerator) == str(denominator):
                raise InvalidValueError(
                    f"The numerator and denominator of the comparison {comparison} are both '{numerator}'. "
                    f"Please choose two different values of the factor '{factor}'."
                )

        for covariate in self.covariates:
            if covariate not in factors:
                raise InvalidValueError(
                    f"The covariate '{covariate}' is not a column of the design matrix. "
                    f'The available columns are: {factors_str}.'
                )
            non_numeric = self._find_non_numeric_value(self.design_mat[covariate])
            if non_numeric is not None:
                raise InvalidValueError(
                    f"The covariate '{covariate}' must contain only numeric values "
                    f"(for example age or dose), but it contains the non-numeric value '{non_numeric}'. "
                    f"To compare groups of samples by '{covariate}', use it in a comparison instead of as a covariate."
                )

        for lrt_factor in self.lrt_factors:
            for part in str(lrt_factor).split(':'):
                part = part.strip()
                poly_match = re.fullmatch(self.POLY_PATTERN, part)
                name = poly_match.group(1) if poly_match else part
                if not re.fullmatch(r'[\w.]+', name):
                    # not a plain factor name (some other R expression) - leave the judgement to R
                    continue
                if name not in factors:
                    raise InvalidValueError(
                        f"The factor '{name}' in the likelihood ratio test term '{lrt_factor}' "
                        f'is not a column of the design matrix. The available factors are: {factors_str}.'
                    )

        # make sure that no two analyses would be saved to the same output file
        self._get_export_filenames()

    def _comparison_factor_name(self, factor: str) -> str:
        return factor

    def _get_export_filenames(self) -> Dict[tuple, str]:
        """
        Return the (filename-safe) output file name of every comparison, covariate, and LRT factor, \
        keyed by ('comparison', comparison), ('covariate', covariate), or ('lrt', lrt_factor).

        :raises InvalidValueError: if two different analyses would be saved to the same file \
        (e.g. the values 'KO/1' and 'KO:1', which are both saved as 'KO_1', or values that differ only by case, \
        which share a file on Windows and macOS).
        """
        sanitize = parsing.sanitize_filename_component
        filenames = {}
        descriptions = {}
        for comparison in self.comparisons:
            factor, numerator, denominator = comparison
            key = ('comparison', tuple(comparison))
            filenames[key] = (
                f'{self.FILENAME_PREFIX}_{sanitize(self._comparison_factor_name(factor))}_'
                f'{sanitize(numerator)}_vs_{sanitize(denominator)}.csv'
            )
            descriptions[key] = f'the comparison {comparison}'
        for covariate in self.covariates:
            key = ('covariate', covariate)
            filenames[key] = f'{self.FILENAME_PREFIX}_{sanitize(covariate)}_covariate.csv'
            descriptions[key] = f"the covariate '{covariate}'"
        for lrt_factor in self.lrt_factors:
            key = ('lrt', lrt_factor)
            lrt_name = sanitize(lrt_factor, {':': self.LRT_INTERACTION_FILENAME_SEP})
            filenames[key] = f'{self.FILENAME_PREFIX}_{lrt_name}_LRT.csv'
            descriptions[key] = f"the likelihood ratio test term '{lrt_factor}'"

        seen = {}
        for key, filename in filenames.items():
            # file names are not case-sensitive on Windows and macOS
            normalized = filename.casefold()
            if normalized in seen:
                other = seen[normalized]
                raise InvalidValueError(
                    f'{descriptions[other]} and {descriptions[key]} would both be saved to the same output file '
                    f"('{filenames[other]}' and '{filename}'). "
                    f'This happens when names differ only in characters that cannot be used in file names '
                    f"(< > : \" / \\ | ? * and control characters, which are replaced with '_'), "
                    f'or only in upper/lower case (which Windows and macOS do not distinguish in file names). '
                    f'Please rename these values in your design matrix so that they can be told apart.'
                )
            seen[normalized] = key
        return filenames

    @staticmethod
    def _get_cache_dir():
        cache_dir = io.get_todays_cache_dir().joinpath(hashlib.sha1(str(time.time_ns()).encode('utf-8')).hexdigest())
        if not cache_dir.exists():
            cache_dir.mkdir(parents=True)
        return cache_dir

    @abc.abstractmethod
    def install_required_packages(self):
        pass

    @abc.abstractmethod
    def create_script(self) -> Path:
        raise NotImplementedError

    @abc.abstractmethod
    def _get_templates(self) -> Tuple[str, str, str, str]:
        raise NotImplementedError

    def create_formula(self, without: Union[List[str], Tuple[str, ...]] = tuple()) -> str:
        return '~ ' + ' + '.join([factor for factor in self.model_factors if factor not in without])

    def get_logging_templates(self, log_path: Path) -> Tuple[str, str]:
        with open(self.TEMPLATES_DIR.joinpath('logging.R')) as f:
            # logging.R is shared with featureCounts, so its placeholder keeps its own quotes
            logging_template = f.read().replace('"$LOGFILE"', parsing.r_quote(log_path.as_posix()))
        with open(self.TEMPLATES_DIR.joinpath('sessioninfo_run.R')) as f:
            sessioninfo_template = f.read()
        return logging_template, sessioninfo_template


class DESeqRunner(DiffExpRunner):
    FILENAME_PREFIX = 'DESeq2'
    LRT_INTERACTION_FILENAME_SEP = '_x_'

    def __init__(
        self,
        data_path: Union[str, Path],
        design_mat_path: Union[str, Path],
        comparisons: Iterable[Tuple[str, str, str]],
        covariates: Iterable[str] = tuple(),
        lrt_factors: Iterable[str] = tuple(),
        model_factors: Iterable[str] = 'auto',
        r_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
        scale_factor_path: Union[str, Path, None] = None,
        cooks_cutoff: bool = True,
        scale_factors_ndims: int = 1,
    ):
        super().__init__(
            data_path, design_mat_path, comparisons, covariates, lrt_factors, model_factors, r_installation_folder
        )
        self.scale_factor_path = scale_factor_path
        self.cooks_cutoff = cooks_cutoff
        if scale_factors_ndims not in [1, 2]:
            raise InvalidValueError('scale_factors_ndims must be 1 or 2!')
        self.scaling_factors_ndims = scale_factors_ndims

    def install_required_packages(self):
        super().install_required_packages()
        installs.install_deseq2(self.r_installation_folder)

    def _get_templates(self) -> Tuple[str, str, str, str]:
        with open(self.TEMPLATES_DIR.joinpath('deseq2_run_parametric.R')) as f:
            run_template = f.read()
        with open(self.TEMPLATES_DIR.joinpath('deseq2_comparison_parametric.R')) as f:
            comparison_template = f.read()
        with open(self.TEMPLATES_DIR.joinpath('deseq2_covariate_parametric.R')) as f:
            covariate_template = f.read()
        with open(self.TEMPLATES_DIR.joinpath('deseq2_lrt_parametric.R')) as f:
            lrt_template = f.read()
        return run_template, comparison_template, covariate_template, lrt_template

    def create_script(self):
        filenames = self._get_export_filenames()
        cache_dir = self._get_cache_dir()
        save_path = cache_dir.joinpath('deseq2_run.R')
        log_path = cache_dir.joinpath('logfile.log')
        run_template, comparison_template, covariate_template, lrt_template = self._get_templates()
        logging_template, sessioninfo_template = self.get_logging_templates(log_path)

        with open(save_path, 'w', encoding='utf-8') as outfile:
            outfile.write(logging_template)

            if self.scale_factor_path is None:
                normfactors_key, normfactors_code = '$NORMFACTORS\n', ''
            else:
                normfactors_key = '$NORMFACTORS'
                scale_factor_path = parsing.r_quote(Path(self.scale_factor_path).as_posix())
                if self.scaling_factors_ndims == 1:
                    normfactors_code = (
                        f'sizeFactors <- read.csv({scale_factor_path}, row.names=1, header=TRUE)\n'
                        'sizeFactors(dds) <- as.numeric(unlist(sizeFactors))\n'
                    )
                else:
                    normfactors_code = (
                        f'normFactors <- read.csv({scale_factor_path}, row.names=1, header=TRUE)\n'
                        'normFactors <- normFactors / exp(rowMeans(log(normFactors)))\n'
                        'normalizationFactors(dds) <- data.matrix(normFactors)\n'
                    )

            run_template = _fill_template(
                run_template,
                {
                    '$COUNT_MATRIX': parsing.r_quote(Path(self.data_path).as_posix()),
                    '$DESIGN_MATRIX': parsing.r_quote(Path(self.design_mat_path).as_posix()),
                    '$FORMULA': self.create_formula(),
                    normfactors_key: normfactors_code,
                },
            )
            outfile.write(run_template)

            cooks = 'TRUE' if self.cooks_cutoff else 'FALSE'
            # pairwise comparisons
            for contrast in self.comparisons:
                export_path = cache_dir.joinpath(filenames[('comparison', tuple(contrast))]).as_posix()
                this_comparison = _fill_template(
                    comparison_template,
                    {
                        '$CONTRAST': ', '.join(parsing.r_quote(item) for item in contrast),
                        '$OUTFILE_NAME': parsing.r_quote(export_path),
                        '$COOKS': cooks,
                    },
                )
                outfile.write(this_comparison)
            # covariates
            for covariate in self.covariates:
                export_path = cache_dir.joinpath(filenames[('covariate', covariate)]).as_posix()
                this_covar = _fill_template(
                    covariate_template,
                    {
                        '$COVARIATE': parsing.r_quote(covariate),
                        '$OUTFILE_NAME': parsing.r_quote(export_path),
                        '$COOKS': cooks,
                    },
                )
                outfile.write(this_covar)
            # likelihood ratio tests
            for lrt_factor in self.lrt_factors:
                export_path = cache_dir.joinpath(filenames[('lrt', lrt_factor)]).as_posix()
                this_lrt = _fill_template(
                    lrt_template,
                    {
                        '$REDUCED': self.create_formula([lrt_factor]),
                        '$OUTFILE_NAME': parsing.r_quote(export_path),
                        '$COOKS': cooks,
                    },
                )
                outfile.write(this_lrt)

            outfile.write(sessioninfo_template)
        return save_path


class LimmaVoomRunner(DiffExpRunner):
    FILENAME_PREFIX = 'LimmaVoom'
    LRT_INTERACTION_FILENAME_SEP = 'X'

    def __init__(
        self,
        data_path: Union[str, Path],
        design_mat_path: Union[str, Path],
        comparisons: Iterable[Tuple[str, str, str]],
        covariates: Iterable[str] = tuple(),
        lrt_factors: Iterable[str] = tuple(),
        model_factors: Iterable[str] = 'auto',
        r_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
        random_effect: Union[str, None] = None,
        quality_weights: bool = False,
    ):
        super().__init__(
            data_path, design_mat_path, comparisons, covariates, lrt_factors, model_factors, r_installation_folder
        )
        self.random_effect = random_effect
        self.quality_weights = quality_weights

    def install_required_packages(self):
        super().install_required_packages()
        installs.install_limma(self.r_installation_folder)

    def _comparison_factor_name(self, factor: str) -> str:
        return generic.sanitize_variable_name(factor)

    @staticmethod
    def _r_syntactic_name(name: str) -> str:
        """
        Mirror R's make.names() (as applied to the design-matrix column names in limma_run_parametric.R) \
        for a coefficient name: letters and digits of any script, '.', and '_' are kept; \
        any other character becomes '.'.
        """
        name = ''.join(char if (char in '._' or unicodedata.category(char)[0] in {'L', 'N'}) else '.' for char in name)
        if not re.match(r'([^\W\d_]|\.(?!\d))', name):
            name = 'X' + name
        return name

    def _comparison_coef_name(self, factor: str, value) -> str:
        return self._r_syntactic_name(f'{self._comparison_factor_name(factor)}{value}')

    def validate_design_inputs(self):
        super().validate_design_inputs()
        # R's make.names() can map two different values of a factor to the same limma coefficient name
        # (e.g. 'KO-1' and 'KO.1' both become 'conditionKO.1'), which would make the contrast ambiguous
        for factor, numerator, denominator in self.comparisons:
            values = [str(val) for val in self.design_mat[factor].drop_nulls().unique().to_list()]
            for value in (numerator, denominator):
                coef = self._comparison_coef_name(factor, value)
                clashes = [
                    val for val in values if val != str(value) and self._comparison_coef_name(factor, val) == coef
                ]
                if clashes:
                    raise InvalidValueError(
                        f"The values '{value}' and '{clashes[0]}' of the factor '{factor}' cannot be told apart "
                        f"by limma, since R converts both of them to the same coefficient name ('{coef}'). "
                        f'Please rename one of these values in your design matrix.'
                    )

    def _generate_interaction_terms(self, interaction):
        # Split the interaction string into individual variables
        variables = interaction.split(':')

        # Check if each variable is numeric or factor (non-numeric) and store this information
        is_numeric = {var: self.design_mat[var].dtype.is_numeric() for var in variables}

        # Initialize a list to hold the names of the interaction terms
        interaction_terms = []

        # Generate interaction terms
        for var in variables:
            if is_numeric[var]:
                # Numeric variables are used as they are
                interaction_terms.append([var])
            else:
                # For non-numeric (factor) variables, create a term for each level except the first (reference)
                levels = sorted(self.design_mat[var].unique())  # Sort levels to ensure consistent order
                interaction_terms.append([f'{var}{level}' for level in levels[1:]])

        # Generate all possible combinations of interaction terms
        # If a variable is numeric, it appears as is; if it's a factor, we get one term per level (except the reference level)
        all_combinations = list(itertools.product(*interaction_terms))

        # Convert combinations into coefficient names
        coefficient_names = {':'.join(comb) for comb in all_combinations}

        return coefficient_names

    def _get_coef_names(self) -> Dict[str, Set[str]]:
        coefs = {}
        for factor in self.model_factors:
            coefs[factor] = set()
            if re.match(self.POLY_PATTERN, factor):
                match = re.match(self.POLY_PATTERN, factor)
                base_factor = match.group(1)
                degree = int(match.group(2))
                for i in range(degree):
                    coefs[factor].add(f'poly({base_factor}, degree = {degree}){i + 1}')
            elif factor in self.design_mat and self.design_mat[factor].dtype.is_numeric():
                coefs[factor].add(factor)
            elif ':' in factor:
                coefs[factor] = self._generate_interaction_terms(factor)
            else:
                values = sorted(self.design_mat[factor].unique())
                for value in values[1:]:
                    coefs[factor].add(f'{factor}{value}')
        coef_makenames = {factor: set(parsing.r_make_names(list(names))) for factor, names in coefs.items()}
        return coef_makenames

    @staticmethod
    def _get_random_effect_code(random_effect: Union[str, None]) -> str:
        if random_effect is None:
            random_effect_fit_code = 'fit <- lmFit(voom_object, design)'
        else:
            random_effect_fit_code = (
                'cor <- duplicateCorrelation(voom_object, design, '
                f'block=design_matrix${random_effect})\n'
                'if (cor$consensus.correlation > 0) { '
                '#only include random effect if the correlation is positive\n'
                f'  fit <- lmFit(voom_object, design, block=design_matrix${random_effect}, '
                'correlation=cor$consensus.correlation)\n'
                '}'
                'else {'
                '   fit <- lmFit(voom_object, design)'
                '}'
            )
        return random_effect_fit_code

    def _get_templates(self) -> Tuple[str, str, str, str]:
        with open(self.TEMPLATES_DIR.joinpath('limma_run_parametric.R')) as f:
            run_template = f.read()
        with open(self.TEMPLATES_DIR.joinpath('limma_comparison_parametric.R')) as f:
            comparison_template = f.read()
        with open(self.TEMPLATES_DIR.joinpath('limma_covariate_parametric.R')) as f:
            covariate_template = f.read()
        with open(self.TEMPLATES_DIR.joinpath('limma_lrt_parametric.R')) as f:
            lrt_template = f.read()
        return run_template, comparison_template, covariate_template, lrt_template

    def create_script(self):
        filenames = self._get_export_filenames()
        cache_dir = self._get_cache_dir()
        save_path = cache_dir.joinpath('limma_run.R')
        log_path = cache_dir.joinpath('logfile.log')
        run_template, comparison_template, covariate_template, lrt_template = self._get_templates()
        logging_template, sessioninfo_template = self.get_logging_templates(log_path)
        coef_names = self._get_coef_names()

        with open(save_path, 'w', encoding='utf-8') as outfile:
            outfile.write(logging_template)

            baselines = {}

            random_effect = None
            factor_names = {}
            factors_str = ''
            for factor in self.design_mat.columns[1:]:
                factor_name = generic.sanitize_variable_name(factor)

                if factor == self.random_effect or factor_name == self.random_effect:
                    random_effect = factor_name
                else:
                    factor_names[factor] = factor_name

                if self.design_mat[factor].dtype.is_numeric():
                    factors_str += f'design_matrix${factor_name} <- as.numeric(design_matrix${factor})\n'
                    baselines[factor] = 0
                else:
                    values = sorted(self.design_mat[factor].unique())
                    baselines[factor] = values[0]
                    values_str = ', '.join([parsing.r_quote(val) for val in values])
                    factors_str += (
                        f'design_matrix${factor_name} <- factor(design_matrix${factor}, levels=c({values_str}))\n'
                    )

            formula = self.create_formula(without=tuple() if random_effect is None else [random_effect])
            random_effect_fit_code = self._get_random_effect_code(random_effect)

            run_template = _fill_template(
                run_template,
                {
                    '$COUNT_MATRIX': parsing.r_quote(Path(self.data_path).as_posix()),
                    '$DESIGN_MATRIX': parsing.r_quote(Path(self.design_mat_path).as_posix()),
                    '$DEFINE_FACTORS': factors_str,
                    '$FORMULA': formula,
                    '$VOOMFUNC': 'voomWithQualityWeights' if self.quality_weights else 'voom',
                    '$RANDOM_EFFECT_FIT': random_effect_fit_code,
                },
            )

            outfile.write(run_template)
            # pairwise comparisons
            for comparison in self.comparisons:
                factor, num, denom = comparison
                if factor not in factor_names:
                    raise InvalidValueError(
                        f"Cannot run the comparison {comparison}: '{factor}' is the random effect of the model, "
                        f'so its values cannot be compared to each other.'
                    )
                export_path = cache_dir.joinpath(filenames[('comparison', tuple(comparison))]).as_posix()
                num_not_baseline = num != baselines[factor]
                denom_not_baseline = denom != baselines[factor]
                # the coefficient names, as R's make.names() spells them in the design matrix
                num_coef = self._comparison_coef_name(factor, num)
                denom_coef = self._comparison_coef_name(factor, denom)
                contrast = f'{num_coef * num_not_baseline}{(" - " + denom_coef) * denom_not_baseline}'
                this_comparison = _fill_template(
                    comparison_template,
                    {'$CONTRAST': parsing.r_quote(contrast), '$OUTFILE_NAME': parsing.r_quote(export_path)},
                )
                outfile.write(this_comparison)
            # covariates
            for covariate in self.covariates:
                export_path = cache_dir.joinpath(filenames[('covariate', covariate)]).as_posix()
                coef = covariate if covariate in self.model_factors else f'poly({covariate}, degree = 1)1'
                this_covar = _fill_template(
                    covariate_template,
                    {'$COEF': parsing.r_quote(coef), '$OUTFILE_NAME': parsing.r_quote(export_path)},
                )
                outfile.write(this_covar)
            # likelihood ratio tests
            for lrt_factor in self.lrt_factors:
                export_path = cache_dir.joinpath(filenames[('lrt', lrt_factor)]).as_posix()
                coefs = 'c(' + ', '.join(parsing.r_quote(coef) for coef in sorted(coef_names[lrt_factor])) + ')'
                this_lrt = _fill_template(
                    lrt_template, {'$COEFS': coefs, '$OUTFILE_NAME': parsing.r_quote(export_path)}
                )
                outfile.write(this_lrt)

            outfile.write(sessioninfo_template)
        return save_path
