.. image:: https://raw.githubusercontent.com/GuyTeichman/RNAlysis/master/docs/source/logo.png
    :target: https://guyteichman.github.io/RNAlysis
    :width: 400
    :alt: logo

**Useful links:** `Documentation <https://guyteichman.github.io/RNAlysis>`_ |
`Source code <https://github.com/GuyTeichman/RNAlysis>`_ |
`Bug reports <https://github.com/GuyTeichman/RNAlysis/issues>`_ | |pipimage| | |versionssupported| | |githubactions| | |coveralls| | |downloads|

----

What is *RNAlysis?*
--------------------

*RNAlysis* is a Python-based software for analyzing RNA sequencing data.
*RNAlysis* allows you to build customized analysis pipelines suiting your specific research questions,
going all the way from exploratory data analysis and data visualization through clustering analysis and gene-set enrichment analysis.

----

What can I do with *RNAlysis*?
---------------------------------

* Filter your gene expression matrices, differential expression tables, fold change data, and tabular data in general.
* Normalize your gene expression matrices
* Visualise, explore and describe your sequencing data
* Find global relationships between sample expression profiles with clustering analyses and dimensionality reduction
* Create and share analysis pipelines
* Perform enrichment analysis with pre-determined Gene Ontology terms/KEGG pathways, or with used-defined attributes
* Perform enrichment analysis on a single ranked list of genes, instead of a test set and a background set

To get an overview of what *RNAlysis* can do, read the `user guide <https://guyteichman.github.io/RNAlysis/build/user_guide.html>`_.

*RNAlysis* supports gene expression matrices and differential expression tables in general, and integrates in particular with Python's *HTSeq-count* and R's *DESeq2*.

----

How do I install it?
---------------------
You can install *RNAlysis* via PyPI.

To install the full version of *RNAlysis* (includes additional features that might not work out-of-the-box on all machines), you should first install `Microsoft Visual C++ 14.0 <https://visualstudio.microsoft.com/visual-cpp-build-tools/>`_ or greater (on Windows computers only), `GraphViz <https://graphviz.org/download/>`_, `R <https://cran.r-project.org/bin/>`_, and `kallisto <https://pachterlab.github.io/kallisto/download>`_.
Then use the following command in your terminal window::

    pip install RNAlysis[all]


To install the basic version of *RNAlysis*, use the following command in your terminal window::

    pip install RNAlysis


You can also install *RNAlysis* with only some of the following additional features:

* `fastq` - adapter trimming and RNA-seq transcript quantification of Fastq files
* `hdbscan` - clustering analysis using the HDBSCAN method
* `single-set` - single-set enrichment analysis using the XL-mHG test
* `randomization` - improved performance for randomization tests

by calling the install command with one or more additional features inside the square brackets, separated by commas. For example::

    pip install RNAlysis[fastq,single-set]


will install the basic version of *RNAlysis*, along with the `fastq` and `single-set` additional features.

----


How do I use it?
---------------------
You can launch the *RNAlysis* software by typing the following command::

    rnalysis-gui

Or through a python console::

    >>> from rnalysis import gui
    >>> gui.run_gui()

Alternatively, you can write Python code that uses *RNAlysis* functions as described in the `user guide <https://guyteichman.github.io/RNAlysis/build/user_guide.html>`_.

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
* `typing_extensions <https://github.com/python/typing_extensions>`_
* `PyQt5 <https://www.riverbankcomputing.com/software/pyqt/>`_
* `qdarkstyle <https://github.com/ColinDuquesnoy/QDarkStyleSheet>`_

----

Credits
-------

How do I cite *RNAlysis*?
**************************
Teichman, G. (2022) RNAlysis: RNA Sequencing analysis software (Python package version 3.1.0).

Development Lead
******************

* Guy Teichman: guyteichman@gmail.com

Contributors
*************

* Dror Cohen
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

..  |githubactions| image:: https://github.com/guyteichman/RNAlysis/actions/workflows/python-package.yml/badge.svg
    :target: https://github.com/GuyTeichman/RNAlysis/actions/workflows/python-package.yml
    :alt: Build status

.. |coveralls| image:: https://coveralls.io/repos/github/GuyTeichman/RNAlysis/badge.svg?branch=master
    :target: https://coveralls.io/github/GuyTeichman/RNAlysis?branch=master
    :alt: Coverage
