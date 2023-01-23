.. highlight:: shell

============
Installation
============

.. |pipimage| image:: https://img.shields.io/pypi/v/rnalysis.svg
        :target: https://pypi.python.org/pypi/rnalysis

Latest version: |pipimage|

*RNAlysis* stand-alone app (most beginner-friendly)
-----------------------------------------------------
A stand-alone version of *RNAlysis* is available for both Windows and MacOS.
This version requires the least amount of setup #TODO

How to install it
^^^^^^^^^^^^^^^^^
You can download the latest stand-alone of *RNAlysis* from the
`GitHub Releases page<https://github.com/GuyTeichman/RNAlysis/releases/latest>`_ ('RNAlysis-X.Y.Z_windows.zip' for Windows, and 'RNAlysis-X.Y.Z_macos.zip' for MacOS).

If you use the *RNAlysis* stand-alone app,


How to run it
^^^^^^^^^^^^^
**On Windows:**


.. image:: /docs/source/installation_screenshots/01b01_open_windows.png
  :width: 600
  :alt: Open *RNAlysis* stand-alone app on Windows - click on "More info"


.. image:: /docs/source/installation_screenshots/01b02_open_windows.png
  :width: 600
  :alt: Open *RNAlysis* stand-alone app on Windows - click on "Run anyway"

The *RNAlysis* app should launch now - this may take a minute or two, so be patient!

**On MacOS:**

The *RNAlysis* app should launch now - this may take a minute or two, so be patient!


Install as a Python package with *pip* (best performance)
----------------------------------------------------------


You can install *RNAlysis* via `pip`_.

How to install it
^^^^^^^^^^^^^^^^^

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

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

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/

How to run it
^^^^^^^^^^^^^

If you installed *RNAlysis* with *pip*, you can open the *RNAlysis* app by executing the command `rnalysis-gui` from your terminal.

Alternatively, you can open the *RNAlysis* app by typing the following code into a Python console::

    >>> from rnalysis import gui
    >>> gui.run_gui()



From sources
------------

The source code for RNAlysis can be downloaded from the `Github repository`_.

How to install it
^^^^^^^^^^^^^^^^^

First, clone the public repository:

.. code-block:: console

    $ git clone git://github.com/GuyTeichman/rnalysis


Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repository: https://github.com/GuyTeichman/RNAlysis


How to run it
^^^^^^^^^^^^^

If you installed *RNAlysis* from source, you can open the *RNAlysis* app by executing the command `rnalysis-gui` from your terminal.

Alternatively, you can open the *RNAlysis* app by typing the following code into a Python console::

    >>> from rnalysis import gui
    >>> gui.run_gui()



