from pathlib import Path
from typing import Literal, Union

from rnalysis.utils import installs, io, parsing


class FeatureCountsRunner:
    TEMPLATES_DIR = Path(__file__).parent.parent.joinpath('data_files/r_templates')

    def __init__(self, kwargs: dict, output_dir: Path,
                 r_installation_folder: Union[str, Path, Literal['auto']] = 'auto'):
        self.kwargs = kwargs
        self.output_dir = output_dir
        self.r_installation_folder = r_installation_folder

    def run(self):
        self._install_rsubread()
        script_path = self._create_script()
        self._run_script(script_path)

    def _install_rsubread(self):
        installs.install_rsubread(self.r_installation_folder)

    def _create_script(self):
        save_path = self.output_dir.joinpath('featurecounts_run.R')
        log_path = self.output_dir.joinpath('logfile.log')

        with open(self.TEMPLATES_DIR.joinpath('rsubread_featurecounts_run_parametric.R')) as f:
            run_template = f.read()
        with open(self.TEMPLATES_DIR.joinpath('logging.R')) as f:
            logging_template = f.read().replace("$LOGFILE", log_path.as_posix())
        with open(self.TEMPLATES_DIR.joinpath('sessioninfo_run.R')) as f:
            sessioninfo_template = f.read()

        with open(save_path, 'w') as outfile:
            outfile.write(logging_template)
            run_template = run_template.replace("$KWARGS", parsing.python_to_r_kwargs(self.kwargs))
            run_template = run_template.replace("$OUTPUT_DIR", self.output_dir.as_posix())
            outfile.write(run_template)
            outfile.write(sessioninfo_template)

        return save_path

    def _run_script(self, script_path: Path):
        io.run_r_script(script_path, self.r_installation_folder)
