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

To get an overview of what *RNAlysis* can do, read the `tutorial <https://guyteichman.github.io/RNAlysis/build/tutorial.html>`_ and the `user guide <https://guyteichman.github.io/RNAlysis/build/user_guide.html>`_.

*RNAlysis* supports gene expression matrices and differential expression tables in general, and integrates in particular with Python's *HTSeq-count* and R's *DESeq2*.

----

How do I install it?
---------------------
You can either install *RNAlysis* as a stand-alone app, or via PyPI.
To learn how to install *RNAlysis*, visit the `Installation page <https://guyteichman.github.io/RNAlysis/build/installation.html>`_.

----


How do I use it?
---------------------
If you installed *RNAlysis* as a stand-alone app, simply open the app ("RNAlysis.exe" on Windows, "RNAlysis.dmg" on MacOS) and wait for it to load (it may take a minute or two, so be patient!).

If you installed *RNAlysis* from PyPi, you can launch *RNAlysis* by typing the following command::

    rnalysis-gui

Or through a python console::

    >>> from rnalysis import gui
    >>> gui.run_gui()

In addition, you can write Python code that uses *RNAlysis* functions as described in the `programmatic interface user guide <https://guyteichman.github.io/RNAlysis/build/user_guide.html>`_.

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
* `defusedxml <https://https://github.com/tiran/defusedxml>`_
* `cutadapt <https://github.com/marcelm/cutadapt>`_

----

Credits
-------

How do I cite *RNAlysis*?
**************************
If you use *RNAlysis* in your research, please cite::

    Teichman, G., Cohen, D., Ganon, O., Dunsky, N., Shani, S., Gingold, H., and Rechavi, O. (2022).
    RNAlysis: analyze your RNA sequencing data without writing a single line of code. BioRxiv 2022.11.25.517851.
    https://doi.org/10.1101/2022.11.25.517851

If you use the *CutAdapt* adapter trimming tool in your research, please cite::

    Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads.
    EMBnet.journal, 17(1), pp. 10-12.
    https://doi.org/10.14806/ej.17.1.200

If you use the *kallisto* RNA sequencing quantification tool in your research, please cite::

    Bray, N., Pimentel, H., Melsted, P. et al.
    Near-optimal probabilistic RNA-seq quantification.
    Nat Biotechnol 34, 525–527 (2016).
    https://doi.org/10.1038/nbt.3519

If you use the *DESeq2* differential expression tool in your research, please cite::

    Love MI, Huber W, Anders S (2014).
    “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.”
    Genome Biology, 15, 550.
    https://doi.org/10.1186/s13059-014-0550-8

If you use the *HDBSCAN* clustering feature in your research, please cite::

     L. McInnes, J. Healy, S. Astels, hdbscan:
    Hierarchical density based clustering In:
    Journal of Open Source Software, The Open Journal, volume 2, number 11. 2017
    https://doi.org/10.1371/journal.pcbi.0030039

If you use the *XL-mHG* single-set enrichment test in your research, please cite::

    Eden, E., Lipson, D., Yogev, S., and Yakhini, Z. (2007).
     Discovering Motifs in Ranked Lists of DNA Sequences. PLOS Comput. Biol. 3, e39.
    https://doi.org/10.1371/journal.pcbi.0030039>doi.org/10.1371/journal.pcbi.0030039</a>

    Wagner, F. (2017). The XL-mHG test for gene set enrichment. ArXiv.
    https://doi.org/10.48550/arXiv.1507.07905


Development Lead
******************

* Guy Teichman: guyteichman@gmail.com

Contributors
*************

* Dror Cohen
* Or Ganon
* Netta Dunsky
* Shachar Shani
* `Mintxoklet <https://github.com/Mintxoklet>`_
* `Bipin Kumar <https://github.com/kbipinkumar>`_
* Matthias Wilm

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

..  |githubactions| image:: https://github.com/guyteichman/RNAlysis/actions/workflows/build_ci.yml/badge.svg
    :target: https://github.com/GuyTeichman/RNAlysis/actions/workflows/build_ci.yml
    :alt: Build status

.. |coveralls| image:: https://coveralls.io/repos/github/GuyTeichman/RNAlysis/badge.svg?branch=master
    :target: https://coveralls.io/github/GuyTeichman/RNAlysis?branch=master
    :alt: Coverage
