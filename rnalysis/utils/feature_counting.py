from pathlib import Path
from typing import Union

from typing_extensions import Literal

from rnalysis.utils import io, parsing


def run_featurecounts_analysis(kwargs: dict, output_dir: Path,
                               r_installation_folder: Union[str, Path, Literal['auto']] = 'auto'):
    install_rsubread(r_installation_folder)
    script_path = create_featurecounts_script(kwargs, output_dir)
    io.run_r_script(script_path, r_installation_folder)


def install_rsubread(r_installation_folder: Union[str, Path, Literal['auto']] = 'auto'):
    script_path = Path.joinpath(Path(__file__).parent, '../data_files/r_templates/rsubread_install.R')
    try:
        io.run_r_script(script_path, r_installation_folder)
    except AssertionError:
        raise AssertionError("Failed to install RSubread. "
                             "Please make sure you have write premission to R's library folder, "
                             "or try to install RSubread manually.")


def create_featurecounts_script(kwargs: dict, output_dir: Path):
    save_path = output_dir.joinpath('featurecounts_run.R')

    with open(
        Path.joinpath(Path(__file__).parent, '../data_files/r_templates/rsubread_featurecounts_run_parametric.R')) as f:
        run_template = f.read()

    with open(save_path, 'w') as outfile:
        run_template = run_template.replace("$KWARGS", parsing.python_to_r_kwargs(kwargs))
        run_template = run_template.replace("$OUTPUT_DIR", output_dir.as_posix())
        outfile.write(run_template)

    return save_path
