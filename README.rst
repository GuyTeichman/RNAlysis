.. figure:: https://raw.githubusercontent.com/GuyTeichman/RNAlysis/master/docs/source/logo.png
    :target: https://guyteichman.github.io/RNAlysis
    :alt: logo


**Useful links:** `Documentation <https://guyteichman.github.io/RNAlysis>`_ |
`Source code <https://github.com/GuyTeichman/RNAlysis>`_ |
`Bug reports <https://github.com/GuyTeichman/RNAlysis/issues>`_ | |pipimage| | |versionssupported| | |travisci| | |downloads|

----

What is *RNAlysis?*
--------------------

*RNAlysis* is a Python-based modular analysis pipeline for RNA sequencing data.
You can use it to normalize, filter and visualize your data, cluster genes based on their expression patterns,
and perform enrichment analysis for both Gene Ontology terms and user-defined attributes.

*RNAlysis* allows you to perform filtering operations and analyses at any order you wish.
You can save or load your progress at any given point; the operations you performed on your data and their order
will be reflected in saved file's name.

RNAlysis works with gene expression matrices and differential expression tables in general, and integrates in particular with Python's *HTSeq-count* and R's *DESeq2*.

----

What can I do with *RNAlysis*?
---------------------------------

* Filter your gene expression matrices, differential expression tables, fold change data, and tabular data in general.
* Normalize your gene expression matrices
* Visualise, explore and describe your sequencing data
* Find global relationships between sample expression profiles with clustering analyses and dimensionality reduction
* Create and share analysis pipelines
* Perform enrichment analysis with pre-determined Gene Ontology terms, or with used-defined attributes
* Perform enrichment analysis on a single ranked list, instead of a test set and a background set

----

How do I install it?
---------------------
You can install *RNAlysis* via PyPI.
Use the following command in the python prompt::

    pip install RNAlysis


----

Dependencies
------------
All of *RNAlysis*'s dependencies can be installed automatically via PyPI.

* numpy
* pandas
* scipy
* matplotlib
* numba
* scikit-learn
* scikit-learn-extra
* hdbscan
* seaborn
* statsmodels
* ipyparallel
* grid_strategy
* pyyaml
* UpSetPlot
* matplotlib-venn
* tissue_enrichment_analysis
* xlmhg
* goatools
* wget

----

Credits
-------

Development Lead
******************

* Guy Teichman: guyteichman@gmail.com

Contributors
*************

* Or Ganon
* Netta Dunsky
* Shachar Shani

----

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage



.. |pipimage| image:: https://img.shields.io/pypi/v/rnalysis.svg
    :target: https://pypi.python.org/pypi/rnalysis
    :alt: PyPI version
.. |downloads| image:: https://pepy.tech/badge/rnalysis
    :target: https://pepy.tech/project/rnalysis
    :alt: downloads
.. |versionssupported| image:: https://img.shields.io/pypi/pyversions/RNAlysis.svg
    :target: https://pypi.python.org/pypi/rnalysis
    :alt: Python versions supported

..  |travisci| image:: https://travis-ci.org/GuyTeichman/RNAlysis.svg?branch=master
    :target: https://travis-ci.org/GuyTeichman/RNAlysis
    :alt: Travis-CI build status
