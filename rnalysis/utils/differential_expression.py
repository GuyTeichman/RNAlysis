import abc
import hashlib
import itertools
import time
from pathlib import Path
from typing import Union, Iterable, Tuple, Literal, List

import pandas as pd

from rnalysis.utils import io, generic, installs, validation, parsing


class DiffExpRunner(abc.ABC):
    TEMPLATES_DIR = Path(__file__).parent.parent.joinpath('data_files/r_templates')

    def __init__(self, data_path: Union[str, Path], design_mat_path: Union[str, Path],
                 comparisons: Iterable[Tuple[str, str, str]],
                 covariates: Iterable[str] = tuple(),
                 lrt_factors: Iterable[str] = tuple(),
                 model_factors: Iterable[str] = 'auto',
                 r_installation_folder: Union[str, Path, Literal['auto']] = 'auto'):
        assert validation.isiterable(comparisons) and validation.isinstanceiter(comparisons, tuple), \
            f"comparisons must be an iterable of tuples, instead got: {comparisons}"
        assert validation.isiterable(covariates), f"covariates must be an iterable, instead got: {covariates}"
        assert validation.isiterable(lrt_factors), f"lrt_factors must be an iterable, instead got: {lrt_factors}"
        self.data_path = data_path
        self.design_mat_path = design_mat_path
        self.comparisons = comparisons
        self.covariates = covariates
        self.lrt_factors = lrt_factors
        self.r_installation_folder = r_installation_folder

        self.design_mat = io.load_table(self.design_mat_path, index_col=0)
        if model_factors == 'auto':
            self.model_factors = parsing.data_to_set(
                list(self.design_mat.columns) + parsing.data_to_list(self.lrt_factors) + parsing.data_to_list(
                    self.covariates))
        else:
            self.model_factors = model_factors
        for factor in itertools.chain(self.lrt_factors, self.covariates, [comp[0] for comp in self.comparisons]):
            assert factor in self.model_factors, f"factor {factor} not found in design matrix columns"

    def run(self) -> Path:
        self.install_required_packages()
        script_path = self.create_script()
        io.run_r_script(script_path, self.r_installation_folder)
        return script_path.parent

    @staticmethod
    def _get_cache_dir():
        cache_dir = io.get_todays_cache_dir().joinpath(hashlib.sha1(str(time.time_ns()).encode('utf-8')).hexdigest())
        if not cache_dir.exists():
            cache_dir.mkdir(parents=True)
        return cache_dir

    @abc.abstractmethod
    def install_required_packages(self):
        raise NotImplementedError

    @abc.abstractmethod
    def create_script(self) -> Path:
        raise NotImplementedError

    @abc.abstractmethod
    def _get_templates(self) -> Tuple[str, str, str, str]:
        raise NotImplementedError

    def create_formula(self, without: Union[List[str], Tuple[str, ...]] = tuple()) -> str:
        return "~ " + " + ".join([factor for factor in self.model_factors if factor not in without])

    def expand_factors(self, factors: str) -> List[str]:
        return factors.split(' + ')  # TODO


class DESeqRunner(DiffExpRunner):
    def install_required_packages(self):
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
        run_template, comparison_template, covariate_template, lrt_template = self._get_templates()

        with open(save_path, 'w') as outfile:
            formula = self.create_formula()
            print(formula)

            run_template = run_template.replace("$COUNT_MATRIX", Path(self.data_path).as_posix())
            run_template = run_template.replace("$DESIGN_MATRIX", (Path(self.design_mat_path).as_posix()))
            run_template = run_template.replace("$FORMULA", formula)

            outfile.write(run_template)
            # pairwise comparisons
            for contrast in self.comparisons:
                export_path = cache_dir.joinpath(f"DESeq2_{contrast[0]}_{contrast[1]}_vs_{contrast[2]}.csv").as_posix()
                this_comparison = comparison_template.replace("$CONTRAST", str(contrast))
                this_comparison = this_comparison.replace("$OUTFILE_NAME", export_path)
                outfile.write(this_comparison)
            # covariates
            for covariate in self.covariates:
                export_path = cache_dir.joinpath(f"DESeq2_{covariate}_covariate.csv").as_posix()
                this_covar = covariate_template.replace("$COVARIATE", covariate)
                this_covar = this_covar.replace("$OUTFILE_NAME", export_path)
                outfile.write(this_covar)
            # likelihood ratio tests
            for lrt_factor in self.lrt_factors:
                export_path = cache_dir.joinpath(f"DESeq2_{lrt_factor.replace(':', '_x_')}_LRT.csv").as_posix()
                reduced_model = self.create_formula([lrt_factor])
                print(reduced_model)
                this_lrt = lrt_template.replace("$REDUCED", reduced_model)
                this_lrt = this_lrt.replace("$OUTFILE_NAME", export_path)
                outfile.write(this_lrt)

        return save_path


class LimmaVoomRunner(DiffExpRunner):
    def __init__(self, data_path: Union[str, Path], design_mat_path: Union[str, Path],
                 comparisons: Iterable[Tuple[str, str, str]],
                 covariates: Iterable[str] = tuple(),
                 lrt_factors: Iterable[str] = tuple(),
                 model_factors: Iterable[str] = 'auto',
                 r_installation_folder: Union[str, Path, Literal['auto']] = 'auto',
                 random_effect: Union[str, None] = None):
        super().__init__(data_path, design_mat_path, comparisons, covariates, lrt_factors, model_factors,
                         r_installation_folder)
        self.random_effect = random_effect

    def install_required_packages(self):
        installs.install_limma(self.r_installation_folder)

    @staticmethod
    def _get_random_effect_code(random_effect: Union[str, None]) -> str:
        if random_effect is None:
            random_effect_fit_code = "fit <- lmFit(voom_object, design)"
        else:
            random_effect_fit_code = (f"cor <- duplicateCorrelation(voom_object, design, block={random_effect})\n"
                                      "if (cor$consensus.correlation > 0) { "
                                      "#only include random effect if the correlation is positive\n"
                                      f"  fit <- lmFit(voom_object, design, block={random_effect}, "
                                      f"correlation=cor$consensus.correlation)\n"
                                      "}"
                                      "else {"
                                      "   fit <- lmFit(voom_object, design)"
                                      "}")
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
        run_template, comparison_template, covariate_template, lrt_template = self._get_templates()

        with open(save_path, 'w') as outfile:
            baselines = {}

            random_effect = None
            factor_names = {}
            factors_str = ''
            for factor in self.design_mat.columns:
                factor_name = generic.sanitize_variable_name(factor)

                if factor == self.random_effect or factor_name == self.random_effect:
                    random_effect = factor_name
                else:
                    factor_names[factor] = factor_name

                if pd.api.types.is_numeric_dtype(self.design_mat[factor]):
                    factors_str += f'design_matrix${factor_name} <- as.numeric(design_matrix${factor})\n'
                    baselines[factor] = 0
                else:
                    values = sorted(self.design_mat[factor].unique())
                    baselines[factor] = values[0]
                    values_str = ', '.join([f'"{val}"' for val in values])
                    factors_str += f'design_matrix${factor_name} <- factor(design_matrix${factor}, levels=c({values_str}))\n'

            formula = self.create_formula()
            random_effect_fit_code = self._get_random_effect_code(random_effect)

            run_template = run_template.replace("$COUNT_MATRIX", Path(self.data_path).as_posix())
            run_template = run_template.replace("$DESIGN_MATRIX", (Path(self.design_mat_path).as_posix()))
            run_template = run_template.replace("$DEFINE_FACTORS", factors_str)
            run_template = run_template.replace("$FORMULA", formula)
            run_template = run_template.replace("$RANDOM_EFFECT_FIT", random_effect_fit_code)

            outfile.write(run_template)
            # pairwise comparisons
            for factor, num, denom in self.comparisons:
                factor_name = factor_names[factor]
                export_path = cache_dir.joinpath(f"LimmaVoom_{factor_name}_{num}_vs_{denom}.csv").as_posix()
                num_not_baseline = num != baselines[factor]
                denom_not_baseline = denom != baselines[factor]
                contrast = f'"{(factor_name + num) * num_not_baseline}{(" - " + factor + denom) * denom_not_baseline}"'
                this_comparison = comparison_template.replace("$CONTRAST", contrast)
                this_comparison = this_comparison.replace("$OUTFILE_NAME", export_path)
                outfile.write(this_comparison)
            # covariates
            for covariate in self.covariates:
                export_path = cache_dir.joinpath(f"LimmaVoom_{covariate[0]}_covariate.csv").as_posix()
                coef = covariate
                this_covar = covariate_template.replace("$COEF", coef)
                this_covar = this_covar.replace("$OUTFILE_NAME", export_path)
                outfile.write(this_covar)
            # likelihood ratio tests
            for lrt_factor in self.lrt_factors:
                export_path = cache_dir.joinpath(f"LimmaVoom_{lrt_factor[0]}_LRT.csv").as_posix()
                # coefficient count starts at 1, and the first coefficient is always "intercept", so start count from 2
                coefs = ""  # TODO
                this_lrt = lrt_template.replace("$COEFS", coefs)
                this_lrt = this_lrt.replace("$OUTFILE_NAME", export_path)
                outfile.write(this_lrt)

        return save_path
