#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['numpy>=1.18.1', 'pandas>=1.0.1', 'matplotlib', 'seaborn>=0.10.0', 'tissue_enrichment_analysis',
                'statsmodels', 'scikit-learn>=0.22', 'ipyparallel', 'grid_strategy', 'pyyaml', 'UpSetPlot',
                'matplotlib-venn', 'scipy', 'goatools', 'wget', 'hdbscan>=0.8', 'scikit-learn-extra', 'xlmhg>=2.5.4',
                'numba', 'pairwisedist>=1.1.0']

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    author="Guy Teichman",
    author_email='guyteichman@gmail.com',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
    description='RNAlysis provides a modular analysis pipeline for RNA sequencing data. '
                'RNAlysis includes various methods for filtering, data visualisation, exploratory analyses, '
                'enrichment anslyses and clustering.',
    install_requires=requirements,
    python_requires='>3.6',
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='RNAlysis',
    name='RNAlysis',
    packages=find_packages(),
    # packages=find_packages(include=['RNAlysis']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/GuyTeichman/RNAlysis',
    version='1.3.4',
    zip_safe=False,
)
