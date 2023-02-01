#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages


def get_extra_requires(path, add_all=True):
    import re
    from collections import defaultdict

    with open(path) as fp:
        extra_deps = defaultdict(set)
        for k in fp:
            if k.strip() and not k.startswith('#'):
                tags = set()
                if ':' in k:
                    k, v = k.split(':')
                    tags.update(vv.strip() for vv in v.split(','))
                tags.add(re.split('[<=>]', k)[0])
                for t in tags:
                    extra_deps[t].add(k)

        # add tag `all` at the end
        if add_all:
            extra_deps['all'] = set(vv for v in extra_deps.values() for vv in v)

    return extra_deps


with open('README.rst', encoding='utf8') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

with open('requirements.txt') as requirements_file:
    requirements = requirements_file.read().split('\n')

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', 'pytest-qt']

extras_require = get_extra_requires('requirements_extra.txt')

setup(
    author="Guy Teichman",
    author_email='guyteichman@gmail.com',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
    description='RNAlysis is an analysis software for RNA sequencing data. '
                'RNAlysis can help to filter, visualize, explore, analyze, and share your data. ',
    install_requires=requirements,
    python_requires='>=3.7, <3.11',
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='RNAlysis',
    name='RNAlysis',
    packages=find_packages(exclude=['tests', 'packaging']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    extras_require=extras_require,
    url='https://github.com/GuyTeichman/RNAlysis',
    version='3.4.2',
    zip_safe=False,
    entry_points={
        'console_scripts': [
            'rnalysis-gui = rnalysis.gui:run_gui',
        ]}
)
