#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import find_packages, setup


def get_extra_requires(path, add_all=True):
    """Parse an optional-dependencies file into a setuptools ``extras_require`` mapping.

    ``path`` is a *standard* pip requirements file: one requirement per line, with blank lines and
    ``#`` comments ignored. A requirement may carry a trailing ``# feature: <tag>[, <tag> ...]``
    comment that groups it under one or more named extras (e.g. ``pip install RNAlysis[fastq]``).
    Every package is additionally exposed as an extra under its own name, and (when ``add_all`` is
    set) an ``all`` extra collects every optional dependency.

    Keeping the file valid pip syntax lets any tool that expects a real requirements file (pip,
    Dependabot's dependency-graph parser, ...) read it without choking on a custom grammar.
    """
    import re
    from collections import defaultdict

    extra_deps = defaultdict(set)
    with open(path) as fp:
        for line in fp:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            requirement, _, comment = line.partition('#')
            requirement = requirement.strip()
            if not requirement:
                continue
            tags = set()
            match = re.search(r'feature\s*:\s*(.+)', comment)
            if match:
                tags.update(tag.strip() for tag in match.group(1).split(',') if tag.strip())
            # expose each package under its own name as an extra too (e.g. RNAlysis[hdbscan])
            tags.add(re.split(r'[<>=!~;\[ ]', requirement)[0].strip())
            for tag in tags:
                extra_deps[tag].add(requirement)

    # add tag `all` at the end
    if add_all:
        extra_deps['all'] = set(dep for deps in extra_deps.values() for dep in deps)

    return extra_deps


if __name__ == '__main__':
    with open('README.rst', encoding='utf8') as readme_file:
        readme = readme_file.read()

    with open('HISTORY.rst') as history_file:
        history = history_file.read()

    with open('requirements.txt') as requirements_file:
        requirements = requirements_file.read().split('\n')

    test_requirements = ['pytest', 'pytest-qt']

    extras_require = get_extra_requires('requirements_extra.txt')

    setup(
        author="Guy Teichman",
        author_email="guyteichman@gmail.com",
        classifiers=[
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Natural Language :: English",
            "Programming Language :: Python :: 3 :: Only",
            "Programming Language :: Python :: 3.11",
            "Programming Language :: Python :: 3.12",
            "Programming Language :: Python :: 3.13",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Scientific/Engineering :: Information Analysis",
            "Topic :: Scientific/Engineering :: Visualization",
        ],
        description="RNAlysis is an analysis software for RNA sequencing data. "
        "RNAlysis can help to filter, visualize, explore, analyze, and share your data. ",
        install_requires=requirements,
        python_requires=">=3.10, <3.15",
        license="MIT license",
        long_description=readme + "\n\n" + history,
        include_package_data=True,
        keywords="RNAlysis",
        name="RNAlysis",
        packages=find_packages(exclude=["tests", "packaging"]),
        test_suite="tests",
        tests_require=test_requirements,
        extras_require=extras_require,
        url="https://github.com/GuyTeichman/RNAlysis",
        version="4.2.0",
        zip_safe=False,
        entry_points={
            "console_scripts": [
                "rnalysis-gui = rnalysis.gui:run_gui",
            ]
        },
    )
