import hashlib
import time
from pathlib import Path
from typing import Union, Iterable, Tuple

from rnalysis.utils import io, generic

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def install_limma(r_installation_folder: Union[str, Path, Literal['auto']] = 'auto'):
    script_path = Path(__file__).parent.parent.joinpath('data_files/r_templates/limma_install.R')
    try:
        io.run_r_script(script_path, r_installation_folder)
    except AssertionError:
        raise AssertionError("Failed to install limma. "
                             "Please make sure you have write premission to R's library folder, "
                             "or try to install limma manually.")


def create_limma_script(data_path: Union[str, Path], design_mat_path: Union[str, Path],
                        comparisons: Iterable[Tuple[str, str, str]]):
    cache_dir = io.get_todays_cache_dir().joinpath(hashlib.sha1(str(time.time_ns()).encode('utf-8')).hexdigest())
    if not cache_dir.exists():
        cache_dir.mkdir(parents=True)
    save_path = cache_dir.joinpath('limma_run.R')

    with open(Path(__file__).parent.parent.joinpath('data_files/r_templates/limma_run_parametric.R')) as f:
        run_template = f.read()
    with open(Path(__file__).parent.parent.joinpath('data_files/r_templates/limma_export_parametric.R')) as f:
        export_template = f.read()

    with open(save_path, 'w') as outfile:
        baselines = {}
        design_mat_df = io.load_csv(design_mat_path, index_col=0)

        factor_names = {}
        factors_str = ''
        for factor in design_mat_df.columns:
            factor_name = generic.sanitize_variable_name(factor)
            factor_names[factor] = factor_name

            values = sorted(design_mat_df[factor].unique())
            values_str = ', '.join([f'"{val}"' for val in values])
            factors_str += f'{factor_name} <- factor(design_matrix${factor}, levels=c({values_str}))\n'
            baselines[factor] = values[0]

        formula = "~ " + " + ".join(factor_names.values())

        run_template = run_template.replace("$COUNT_MATRIX", Path(data_path).as_posix())
        run_template = run_template.replace("$DESIGN_MATRIX", (Path(design_mat_path).as_posix()))
        run_template = run_template.replace("$DEFINE_FACTORS", factors_str)
        run_template = run_template.replace("$FORMULA", formula)

        outfile.write(run_template)

        for factor, num, denom in comparisons:
            factor_name = factor_names[factor]
            export_path = cache_dir.joinpath(f"LimmaVoom_{factor_name}_{num}_vs_{denom}.csv").as_posix()
            num_not_baseline = num != baselines[factor]
            denom_not_baseline = denom != baselines[factor]
            contrast = f'"{(factor_name + num) * num_not_baseline}{(" - " + factor + denom) * denom_not_baseline}"'
            this_export = export_template.replace("$CONTRAST", contrast)
            this_export = this_export.replace("$OUTFILE_NAME", export_path)

            outfile.write(this_export)

    return save_path


def run_limma_analysis(data_path: Union[str, Path], design_mat_path: Union[str, Path],
                       comparisons: Iterable[Tuple[str, str, str]],
                       r_installation_folder: Union[str, Path, Literal['auto']] = 'auto'):
    install_limma(r_installation_folder)
    script_path = create_limma_script(data_path, design_mat_path, comparisons)
    io.run_r_script(script_path, r_installation_folder)
    return script_path.parent


def install_deseq2(r_installation_folder: Union[str, Path, Literal['auto']] = 'auto'):
    script_path = Path(__file__).parent.parent.joinpath('data_files/r_templates/deseq2_install.R')
    try:
        io.run_r_script(script_path, r_installation_folder)
    except AssertionError:
        raise AssertionError("Failed to install DESeq2. "
                             "Please make sure you have write premission to R's library folder, "
                             "or try to install DESeq2 manually.")


def create_deseq2_script(data_path: Union[str, Path], design_mat_path: Union[str, Path],
                         comparisons: Iterable[Tuple[str, str, str]]):
    cache_dir = io.get_todays_cache_dir().joinpath(hashlib.sha1(str(time.time_ns()).encode('utf-8')).hexdigest())
    if not cache_dir.exists():
        cache_dir.mkdir(parents=True)
    save_path = cache_dir.joinpath('deseq2_run.R')

    with open(Path(__file__).parent.parent.joinpath('data_files/r_templates/deseq2_run_parametric.R')) as f:
        run_template = f.read()
    with open(Path(__file__).parent.parent.joinpath('data_files/r_templates/deseq2_export_parametric.R')) as f:
        export_template = f.read()

    with open(save_path, 'w') as outfile:
        design_mat_df = io.load_csv(design_mat_path, index_col=0)
        formula = "~ " + " + ".join(design_mat_df.columns)

        run_template = run_template.replace("$COUNT_MATRIX", Path(data_path).as_posix())
        run_template = run_template.replace("$DESIGN_MATRIX", (Path(design_mat_path).as_posix()))
        run_template = run_template.replace("$FORMULA", formula)

        outfile.write(run_template)

        for contrast in comparisons:
            export_path = cache_dir.joinpath(f"DESeq2_{contrast[0]}_{contrast[1]}_vs_{contrast[2]}.csv").as_posix()
            this_export = export_template.replace("$CONTRAST", str(contrast))
            this_export = this_export.replace("$OUTFILE_NAME", export_path)

            outfile.write(this_export)

    return save_path


def run_deseq2_analysis(data_path: Union[str, Path], design_mat_path: Union[str, Path],
                        comparisons: Iterable[Tuple[str, str, str]],
                        r_installation_folder: Union[str, Path, Literal['auto']] = 'auto'):
    install_deseq2(r_installation_folder)
    script_path = create_deseq2_script(data_path, design_mat_path, comparisons)
    io.run_r_script(script_path, r_installation_folder)
    return script_path.parent
