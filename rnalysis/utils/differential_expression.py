import hashlib
import time
from pathlib import Path
from typing import Union, Iterable, Tuple

from rnalysis.utils import io

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def install_deseq2(r_installation_folder: Union[str, Path, Literal['auto']] = 'auto'):
    script_path = Path.joinpath(Path(__file__).parent, '../data_files/r_templates/deseq2_install.R')
    try:
        io.run_r_script(script_path, r_installation_folder)
    except AssertionError:
        raise AssertionError("Failed to install DESeq2. "
                             "Please make sure you have write premission to R's library folder, "
                             "or try to install DESeq2 manually.")


def create_deseq2_script(data: Union[str, Path], design_matrix: Union[str, Path],
                         comparisons: Iterable[Tuple[str, str, str]]):
    cache_dir = io.get_todays_cache_dir().joinpath(hashlib.sha1(str(time.time_ns()).encode('utf-8')).hexdigest())
    if not cache_dir.exists():
        cache_dir.mkdir(parents=True)
    save_path = cache_dir.joinpath('deseq2_run.R')

    with open(Path.joinpath(Path(__file__).parent, '../data_files/r_templates/deseq2_run_parametric.R')) as f:
        run_template = f.read()
    with open(Path.joinpath(Path(__file__).parent, '../data_files/r_templates/deseq2_export_parametric.R')) as f:
        export_template = f.read()

    with open(save_path, 'w') as outfile:
        design_mat_df = io.load_csv(design_matrix, index_col=0)
        formula = f"~ " + " + ".join(design_mat_df.columns)

        run_template = run_template.replace("$COUNT_MATRIX", Path(data).as_posix())
        run_template = run_template.replace("$DESIGN_MATRIX", (Path(design_matrix).as_posix()))
        run_template = run_template.replace("$FORMULA", formula)

        outfile.write(run_template)

        for contrast in comparisons:
            export_path = cache_dir.joinpath(f"DESeq2_{contrast[0]}_{contrast[1]}_vs_{contrast[2]}.csv").as_posix()
            this_export = export_template.replace("$CONTRAST", str(contrast))
            this_export = this_export.replace("$OUTFILE_NAME", export_path)

            outfile.write(this_export)

    return save_path


def run_deseq2_analysis(data: Union[str, Path], design_matrix: Union[str, Path],
                        comparisons: Iterable[Tuple[str, str, str]],
                        r_installation_folder: Union[str, Path, Literal['auto']] = 'auto'):
    install_deseq2(r_installation_folder)
    script_path = create_deseq2_script(data, design_matrix, comparisons)
    io.run_r_script(script_path, r_installation_folder)
    return script_path.parent
