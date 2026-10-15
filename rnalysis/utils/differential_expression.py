import abc
import hashlib
import itertools
import re
import time
from pathlib import Path
from typing import Dict, Iterable, List, Literal, Set, Tuple, Union

import polars as pl

from rnalysis.exceptions import InvalidTypeError, InvalidValueError
from rnalysis.utils import generic, installs, io, parsing, validation


class DiffExpRunner(abc.ABC):
    POLY_PATTERN = r'poly\s*\(\s*(\w+)\s*,\s*degree\s*=\s*(\d+)\s*\)'
    TEMPLATES_DIR = Path(__file__).parent.parent.joinpath('data_files/r_templates')

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
        observed = column.drop_nulls().unique().to_list()
        if column.dtype == pl.String:
            # the value is written verbatim into the R script, so it must match exactly (case and whitespace included)
            if str(value) in observed:
                return True
            # Polars reads 'nan'/'NaN' as missing values, while R's read.table() keeps them as literal strings
            return column.null_count() > 0 and str(value) in {'nan', 'NaN'}
        if column.dtype.is_numeric():
            if str(value) in {str(item) for item in observed}:
                return True
            try:
                num_value = float(value)
            except (TypeError, ValueError):
                return False
            return any(num_value == float(item) for item in observed)
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
            logging_template = f.read().replace('$LOGFILE', log_path.as_posix())
        with open(self.TEMPLATES_DIR.joinpath('sessioninfo_run.R')) as f:
            sessioninfo_template = f.read()
        return logging_template, sessioninfo_template


class DESeqRunner(DiffExpRunner):
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
        cache_dir = self._get_cache_dir()
        save_path = cache_dir.joinpath('deseq2_run.R')
        log_path = cache_dir.joinpath('logfile.log')
        run_template, comparison_template, covariate_template, lrt_template = self._get_templates()
        logging_template, sessioninfo_template = self.get_logging_templates(log_path)

        with open(save_path, 'w') as outfile:
            outfile.write(logging_template)

            formula = self.create_formula()

            run_template = run_template.replace('$COUNT_MATRIX', Path(self.data_path).as_posix())
            run_template = run_template.replace('$DESIGN_MATRIX', (Path(self.design_mat_path).as_posix()))
            run_template = run_template.replace('$FORMULA', formula)

            if self.scale_factor_path is None:
                run_template = run_template.replace('$NORMFACTORS\n', '')
            else:
                if self.scaling_factors_ndims == 1:
                    run_template = run_template.replace(
                        '$NORMFACTORS',
                        f"sizeFactors <- read.csv('{Path(self.scale_factor_path).as_posix()}', "
                        'row.names=1, header=TRUE)\n'
                        'sizeFactors(dds) <- as.numeric(unlist(sizeFactors))\n',
                    )
                elif self.scaling_factors_ndims == 2:
                    run_template = run_template.replace(
                        '$NORMFACTORS',
                        f"normFactors <- read.csv('{Path(self.scale_factor_path).as_posix()}', "
                        'row.names=1, header=TRUE)\n'
                        'normFactors <- normFactors / exp(rowMeans(log(normFactors)))\n'
                        'normalizationFactors(dds) <- data.matrix(normFactors)\n',
                    )

            outfile.write(run_template)

            cooks = 'TRUE' if self.cooks_cutoff else 'FALSE'
            # pairwise comparisons
            for contrast in self.comparisons:
                export_path = cache_dir.joinpath(f'DESeq2_{contrast[0]}_{contrast[1]}_vs_{contrast[2]}.csv').as_posix()
                this_comparison = comparison_template.replace('$CONTRAST', str(contrast))
                this_comparison = this_comparison.replace('$OUTFILE_NAME', export_path).replace('$COOKS', cooks)
                outfile.write(this_comparison)
            # covariates
            for covariate in self.covariates:
                export_path = cache_dir.joinpath(f'DESeq2_{covariate}_covariate.csv').as_posix()
                this_covar = covariate_template.replace('$COVARIATE', covariate)
                this_covar = this_covar.replace('$OUTFILE_NAME', export_path).replace('$COOKS', cooks)
                outfile.write(this_covar)
            # likelihood ratio tests
            for lrt_factor in self.lrt_factors:
                export_path = cache_dir.joinpath(f'DESeq2_{lrt_factor.replace(":", "_x_")}_LRT.csv').as_posix()
                reduced_model = self.create_formula([lrt_factor])
                this_lrt = lrt_template.replace('$REDUCED', reduced_model)
                this_lrt = this_lrt.replace('$OUTFILE_NAME', export_path).replace('$COOKS', cooks)
                outfile.write(this_lrt)

            outfile.write(sessioninfo_template)
        return save_path


class LimmaVoomRunner(DiffExpRunner):
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
        cache_dir = self._get_cache_dir()
        save_path = cache_dir.joinpath('limma_run.R')
        log_path = cache_dir.joinpath('logfile.log')
        run_template, comparison_template, covariate_template, lrt_template = self._get_templates()
        logging_template, sessioninfo_template = self.get_logging_templates(log_path)
        coef_names = self._get_coef_names()

        with open(save_path, 'w') as outfile:
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
                    values_str = ', '.join([f'"{val}"' for val in values])
                    factors_str += (
                        f'design_matrix${factor_name} <- factor(design_matrix${factor}, levels=c({values_str}))\n'
                    )

            formula = self.create_formula(without=tuple() if random_effect is None else [random_effect])
            random_effect_fit_code = self._get_random_effect_code(random_effect)

            run_template = run_template.replace('$COUNT_MATRIX', Path(self.data_path).as_posix())
            run_template = run_template.replace('$DESIGN_MATRIX', (Path(self.design_mat_path).as_posix()))
            run_template = run_template.replace('$DEFINE_FACTORS', factors_str)
            run_template = run_template.replace('$FORMULA', formula)
            run_template = run_template.replace(
                '$VOOMFUNC', 'voomWithQualityWeights' if self.quality_weights else 'voom'
            )
            run_template = run_template.replace('$RANDOM_EFFECT_FIT', random_effect_fit_code)

            outfile.write(run_template)
            # pairwise comparisons
            for factor, num, denom in self.comparisons:
                factor_name = factor_names[factor]
                export_path = cache_dir.joinpath(f'LimmaVoom_{factor_name}_{num}_vs_{denom}.csv').as_posix()
                num_not_baseline = num != baselines[factor]
                denom_not_baseline = denom != baselines[factor]
                contrast = f'"{(factor_name + num) * num_not_baseline}{(" - " + factor + denom) * denom_not_baseline}"'
                this_comparison = comparison_template.replace('$CONTRAST', contrast)
                this_comparison = this_comparison.replace('$OUTFILE_NAME', export_path)
                outfile.write(this_comparison)
            # covariates
            for covariate in self.covariates:
                export_path = cache_dir.joinpath(f'LimmaVoom_{covariate}_covariate.csv').as_posix()
                coef = covariate if covariate in self.model_factors else f'poly({covariate}, degree = 1)1'
                this_covar = covariate_template.replace('$COEF', coef)
                this_covar = this_covar.replace('$OUTFILE_NAME', export_path)
                outfile.write(this_covar)
            # likelihood ratio tests
            for lrt_factor in self.lrt_factors:
                export_path = cache_dir.joinpath(f'LimmaVoom_{lrt_factor.replace(":", "X")}_LRT.csv').as_posix()
                coefs = 'c(' + ', '.join(sorted([f'"{coef}"' for coef in coef_names[lrt_factor]])) + ')'
                this_lrt = lrt_template.replace('$COEFS', coefs)
                this_lrt = this_lrt.replace('$OUTFILE_NAME', export_path)
                outfile.write(this_lrt)

            outfile.write(sessioninfo_template)
        return save_path
