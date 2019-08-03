#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['numpy', 'pandas', 'matplotlib', 'seaborn', 'tissue_enrichment_analysis', 'statsmodels', 'scikit-learn']
# requirements = ['numpy', 'pandas', 'matplotlib', 'seaborn', 'tissue_enrichment_analysis', 'statsmodels',
# 'scikit-learn', 'matplotlib-venn', 'simple-venn']

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    author="Guy Teichman",
    author_email='guyteichman@gmail.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description='RNAlysis is a python package providing modular analysis pipeline for RNA sequencing data from '
                'C. elegans. The package includes various filtering methods, data visualisation, clustering analyses, '
                'enrichment anslyses and exploratory analyses.',
    install_requires=requirements,
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
    version='1.1.0',
    zip_safe=False,
)
