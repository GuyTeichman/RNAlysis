.. highlight:: shell

============
Installation
============


Stable release
--------------
.. |pipimage| image:: https://img.shields.io/pypi/v/rnalysis.svg
        :target: https://pypi.python.org/pypi/rnalysis

current version: |pipimage|

You can install *RNAlysis* via PyPI.
To install the basic version of *RNAlysis*, use the following command in your terminal window::

    pip install RNAlysis


To install the full version of *RNAlysis* (includes additional features that might not work out-of-the-box on all machines), you should first install `GraphViz <https://graphviz.org/download/>`_, and `Microsoft Visual C++ 14.0 <https://visualstudio.microsoft.com/visual-cpp-build-tools/>`_ or greater (on Windows computers only).
Then use the following command in your terminal window::

    pip install RNAlysis[all]



If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


From sources
------------

The source code for RNAlysis can be downloaded from the `Github repository`_.

First, clone the public repository:

.. code-block:: console

    $ git clone git://github.com/GuyTeichman/rnalysis


Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repository: https://github.com/GuyTeichman/RNAlysis
