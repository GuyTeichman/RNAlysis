.. figure:: https://raw.githubusercontent.com/GuyTeichman/RNAlysis/master/docs/source/logo.png
    :target: https://guyteichman.github.io/RNAlysis
    :alt: logo


**Useful links:** `Documentation <https://guyteichman.github.io/RNAlysis>`_ |
`Source code <https://github.com/GuyTeichman/RNAlysis>`_ |
`Bug reports <https://github.com/GuyTeichman/RNAlysis/issues>`_ | |pipimage| | |versionssupported| | |githubactions| | |coveralls| | |downloads|

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

* `numpy <https://numpy.org/>`_
* `pandas <https://pandas.pydata.org/>`_
* `scipy <https://www.scipy.org/>`_
* `matplotlib <https://matplotlib.org/>`_
* `numba <http://numba.pydata.org/>`_
* `requests <https://github.com/psf/requests/>`_
* `scikit-learn <https://scikit-learn.org/>`_
* `scikit-learn-extra <https://github.com/scikit-learn-contrib/scikit-learn-extra>`_
* `hdbscan <https://github.com/scikit-learn-contrib/hdbscan>`_
* `seaborn <https://seaborn.pydata.org/>`_
* `statsmodels <https://www.statsmodels.org/>`_
* `joblib <https://joblib.readthedocs.io/en/latest/>`_
* `tqdm <https://github.com/tqdm/tqdm>`_
* `appdirs <https://github.com/ActiveState/appdirs>`_
* `grid_strategy <https://github.com/matplotlib/grid-strategy>`_
* `pyyaml <https://github.com/yaml/pyyaml>`_
* `UpSetPlot <https://github.com/jnothman/UpSetPlot>`_
* `matplotlib-venn <https://github.com/konstantint/matplotlib-venn>`_
* `xlmhg <https://github.com/flo-compbio/xlmhg>`_
* `pairwisedist <https://github.com/GuyTeichman/pairwisedist/>`_

----

Credits
-------

How do I cite *RNAlysis*?
**************************
Teichman, G. (2021) RNAlysis: RNA Sequencing analysis pipeline (Python package version 2.0.0).

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
    :alt: Downloads
.. |versionssupported| image:: https://img.shields.io/pypi/pyversions/RNAlysis.svg
    :target: https://pypi.python.org/pypi/rnalysis
    :alt: Python versions supported

..  |githubactions| image:: https://github.com/guyteichman/RNAlysis/actions/workflows/python-package.yml/badge.svg?branch=V1.4.0_unstable
    :target: https://github.com/GuyTeichman/RNAlysis/actions/workflows/python-package.yml
    :alt: Build status

.. |coveralls| image:: https://coveralls.io/repos/github/GuyTeichman/RNAlysis/badge.svg?branch=V1.4.0_unstable
    :target: https://coveralls.io/github/GuyTeichman/RNAlysis?branch=master
    :alt: Coverage
