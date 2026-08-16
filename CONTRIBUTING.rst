.. highlight:: shell

====================================
How can I contribute to *RNAlysis*?
====================================

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/GuyTeichman/rnalysis/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

*RNAlysis* could always use more documentation, whether as part of the
official *RNAlysis* docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/GuyTeichman/rnalysis/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)


Implement Features
~~~~~~~~~~~~~~~~~~

You are welcome to implement new features that you think would be useful!
If you are looking for ideas, check out the Ideas/To-Do list below:

Ideas/To-Do
~~~~~~~~~~~~~~~~~~

This is a rather unsorted list of features that would be nice to have,
of things that could be improved in the source code, and of possible algorithmic improvements:

* Interactive scatter plots (for example - click on a point to show/hide it's gene/sample name)
* Implement the GSEA single-list enrichment algorithm
* Accept per-gene scaling factors in addition to per-sample scaling factors
* Distplot for count matrices
* Support for single-cell transcriptomics analysis
* Support for ATAC-seq analysis
* Generate sashimi plots
* Generate genome browser-like plots
* Merge tables row/column-wise
* Support for additional formats of customized annotations
* Correlation plots
* Rank-sum test comparing data from two different columns/tables

*RNAlysis* Design Philosophy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As a guiding principle for future development of *RNAlysis*, see this short list of the general design principles *RNAlysis* should follow:

* **Universality** - *RNAlysis* should be usable on most computers. This means new features should be supported by commonly-used operating systems (Windows, Linux, MacOS) and currently supported Python versions. In the best case scenario, users who open *RNAlysis* should never ask themselves "which of these functions/modules are supported by my operating system?"
* **Graphical and Programmatic support** - The majority of *RNAlysis* features should be usable both in the Graphical User Interface, and in Python scripts.
* **Customization and clarity** - Functions and algorithms implemented in *RNAlysis* should give user a high degree of control over the function/algorithm parameters. This serves two purposes: firstly, it allows more advanced users to modify analysis parameters to suit their needs. Second, and perhaps most importantly, specifying and documenting all possible analysis parameters provides a high degree of transparency. Knowing exactly what parameters can be changed, and what are the default values, makes accurate and transparent reporting of analysis much easier. For example, when running GO enrichment, users can know and report exactly how the GO annotations are fetched, filtered, and annotated, which statistical test is used for analysis, etc.

Get Started!
------------

Ready to contribute? Here's how to set up `rnalysis` for local development.

1. Fork the `rnalysis` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/rnalysis.git

3. Install your local copy into a virtual environment::

    $ cd rnalysis/
    $ python -m venv venv
    $ source venv/bin/activate    # on Windows: venv\Scripts\activate
    $ pip install -e .[all]
    $ pip install -r requirements_dev.txt
    $ pre-commit install    # enable the auto-formatting git hooks (Ruff); recommended

4. Create a branch for local development, based on the ``development`` branch
   (``development`` is the integration branch; ``master`` only receives ``development``
   at a version release)::

    $ git checkout development
    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that the tests still pass::

    $ pytest tests/

   While iterating you can run a single module, e.g. ``pytest tests/test_filtering.py``.
   Note that the full suite is long and some tests require R, kallisto, bowtie2, or network
   access. New code should come with tests.

6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website, targeting the ``development`` branch.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. Open the pull request against the ``development`` branch.
4. The pull request should work across the supported Python versions (currently 3.10 - 3.14)
   on Windows, macOS, and Linux. GitHub Actions runs this matrix automatically on every pull
   request; make sure all jobs pass. Coverage is reported at
   https://coveralls.io/github/GuyTeichman/RNAlysis.
5. The files your pull request touches must pass the pre-commit hooks (Ruff formatting/import
   sorting, end-of-file/whitespace fixups, YAML/TOML validity). GitHub Actions checks this
   automatically, scoped to just the files the PR changed. Installing the git hook with
   ``pre-commit install`` (see "Get Started!" above) catches this before you push.

Tips
----

To run a subset of tests::

$ pytest tests/test_filtering.py


Deploying
---------

A reminder for the maintainers on how to deploy.
Make sure all your changes are committed (including an entry in HISTORY.rst), and that
``development`` has been merged into ``master``.
Then, from ``master``, run::

$ bumpversion patch # possible: major / minor / patch
$ git push
$ git push --tags

GitHub Actions will then build and publish the release (PyPI package and standalone installers)
when the tests pass.
