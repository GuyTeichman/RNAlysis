import shutil
import subprocess
import unittest
from pathlib import Path

TEMPLATE_DIR = Path(__file__).parents[1] / 'rnalysis' / 'data_files' / 'r_templates'
R_INSTALL_TEMPLATES = [
    ('limma_install.R', 'limma'),
    ('deseq2_install.R', 'DESeq2'),
    ('rsubread_install.R', 'Rsubread'),
]
RSCRIPT = shutil.which('Rscript')


class RInstallTemplateTests(unittest.TestCase):
    def test_install_templates_fail_if_package_is_still_unavailable(self):
        for script_name, package_name in R_INSTALL_TEMPLATES:
            with self.subTest(script_name=script_name):
                script = TEMPLATE_DIR.joinpath(script_name).read_text(encoding='utf-8')
                install_call = f'BiocManager::install("{package_name}"'
                load_check = f'if (!requireNamespace("{package_name}", quietly = TRUE))'

                self.assertIn(install_call, script)
                self.assertIn(load_check, script)
                self.assertLess(script.index(install_call), script.index(load_check))
                self.assertIn(f'stop("Failed to install {package_name}. The package is still unavailable.")', script)

    def test_install_templates_raise_download_timeout_and_remain_balanced(self):
        for script_name, _ in R_INSTALL_TEMPLATES:
            with self.subTest(script_name=script_name):
                script = TEMPLATE_DIR.joinpath(script_name).read_text(encoding='utf-8')

                self.assertIn('options(timeout = max(300, getOption("timeout")))', script)
                self.assertEqual(script.count('{'), script.count('}'))
                self.assertEqual(script.count('('), script.count(')'))

    @unittest.skipUnless(RSCRIPT, 'Rscript is required to parse R install templates')
    def test_install_templates_parse_as_r(self):
        for script_name, _ in R_INSTALL_TEMPLATES:
            with self.subTest(script_name=script_name):
                subprocess.run(
                    [
                        RSCRIPT,
                        '-e',
                        'parse(file=commandArgs(trailingOnly = TRUE)[1])',
                        TEMPLATE_DIR / script_name,
                    ],
                    check=True,
                    capture_output=True,
                    text=True,
                )


if __name__ == '__main__':
    unittest.main()
