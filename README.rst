====================================================
RNAlysis - *C. elegans* RNA sequencing data analysis
====================================================


.. image:: https://img.shields.io/pypi/v/rnalysis.svg
        :target: https://pypi.python.org/pypi/rnalysis

What is it?
-----------

RNAlysis is a python package providing modular analysis pipeline for RNA sequencing data from *C. elegans*.
The package includes various filtering methods, data visualisation, clustering analyses, enrichment anslyses and
exploratory analyses.

RNAlysis allows you to perform filtering and analyses at any order you wish.
It has the ability to save or load your progress at any given phase,
Wand document the order of operations you performed in the saved file names.

RNAlysis work with HTCount files from python's HTSeq and with DESeq files from R's DESeq2.

* Documentation: https://guyteichman.github.io/RNAlysis


Main Features
-------------

* Filtering of DESeq and HTCount files by counts, biotype, change direction, significance and more.
* Normalization of HTCount files
* Exploratory plots
* GO, Tissue and Phenotype enrichment analysis
* Enrichment analysis for custom feature characteristics

Dependencies
------------

* numpy
* pandas
* matplotlib
* seaborn
* tissue_enrichment_analysis
* statsmodels
* sklearn

Where to get it
---------------
Use the following command in the python prompt::

    pip install RNAlysis


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
