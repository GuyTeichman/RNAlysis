.. figure:: https://raw.githubusercontent.com/GuyTeichman/RNAlysis/master/docs/source/logo.png
        :target: https://guyteichman.github.io/RNAlysis

.. |pipimage| image:: https://img.shields.io/pypi/v/rnalysis.svg
        :target: https://pypi.python.org/pypi/rnalysis
.. |downloads| image:: https://pepy.tech/badge/rnalysis
        :target: https://pepy.tech/project/rnalysis

..  |travisci| image:: https://travis-ci.org/GuyTeichman/RNAlysis.svg?branch=master
    :target: https://travis-ci.org/GuyTeichman/RNAlysis

**Useful links:** `Documentation <https://guyteichman.github.io/RNAlysis>`_ |
`Source code <https://github.com/GuyTeichman/RNAlysis>`_ |
`Bug reports <https://github.com/GuyTeichman/RNAlysis/issues>`_ | |pipimage| | |travisci| | |downloads|

----

What is RNAlysis?
------------------

RNAlysis is a python package providing modular analysis pipeline for RNA sequencing data.
The package includes various filtering methods, data visualisation, clustering analyses, enrichment anslyses and
exploratory analyses.

RNAlysis allows you to perform filtering and analyses at any order you wish.
It has the ability to save or load your progress at any given phase,
Wand document the order of operations you performed in the saved file names.

RNAlysis works with sequencing count matrices and differential expression output files in general, and integrates in particular with python's HTSeq-count and R's DESeq2.

----

Main Features
-------------

* Filter your read count matrices, fold change data, differential expression tables, and tabular data in general.
* Normalize your read count data
* Visualise, explore and describe your sequencing data
* Find global relationships between sample expression profiles with clustering and dimensionality reduction
* Create and share analysis pipelines
* Perform enrichment analysis with pre-determined Gene Ontology terms, or with used-defined attributes

----

Dependencies
------------

* numpy
* pandas
* scipy
* matplotlib
* seaborn
* tissue_enrichment_analysis
* statsmodels
* scikit-learn
* ipyparallel
* grid_strategy
* Distance
* pyyaml
* UpSetPlot
* matplotlib-venn

----

Where to get it
---------------
Use the following command in the python prompt::

    pip install RNAlysis


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
