####################################
Frequently asked questions
####################################

*************
Installation
*************


I get an error message when I try to install *RNAlysis*. what can I do?
====================================================================================
Try one of the following methods to solve your problem:

1. Make sure that your computer has an internet connection while you are installing *RNAlysis*.
2. Instead of installing the full version of *RNAlysis*, try installing the basic version first by using the command:
   ::

        pip install RNAlysis
   If that worked, your system might have trouble installing one of the extended features of *RNAlysis*.
   Now you can find out which feature is causing issues by installing them one by one, by running the following commands one after the other:
    ::

        pip install RNAlysis[fastq]
        pip install RNAlysis[hdbscan]
        pip install RNAlysis[single-set]
        pip install RNAlysis[randomization]
   See the next FAQ items to solve issues with a specific extended feature.
3. Make sure your version of Python is supported by *RNAlysis*. Currently, any Python version between 3.7-3.10 should work, but 3.8 would probably work best.
4. Try `creating a new Python environment <https://towardsdatascience.com/virtual-environments-104c62d48c54?gi=40d0a7444922>`_ for *RNAlysis*, or `re-installing Python <http://docs.python-guide.org/en/latest/starting/installation/>`_.

If you are still having trouble installing *RNAlysis*, please submit a `bug report <https://github.com/GuyTeichman/RNAlysis/issues>`_ or `contact me <mailto:guyteichman@gmail.com>`_.


I get an error message when I try to install the *RNAlysis* `hdbscan` feature. what can I do?
=========================================================================================================
The HDBSCAN module makes heavy use of pre-compiled code for performance, but that unfortunately means that installation can cause issues, particularly on computers running Windows.
To successfully install HDBSCAN, try one of the following solutions:

1. If your computer has a Windows operating system, try installing `Microsoft Visual C++ 14.0 <https://visualstudio.microsoft.com/visual-cpp-build-tools/>`_.
After the installation is done, restart your computer, and try to install the HDBSCAN module again.
2. Instead of installing Python normally, try installing the `Anaconda or Miniconda <https://www.edureka.co/blog/python-anaconda-tutorial/>`_ Python distribution.
You can then install HDBSCAN via anaconda and *RNAlysis* via pip, like so:
   ::

        conda install hdbscan
        pip install RNAlysis[all]

I get an error message when I try to install the *RNAlysis* `randomization` feature. what can I do?
=========================================================================================================
The randomization module makes heavy use of the Python package `numba` for performance, but that unfortunately means that installation can cause issues, particularly on computers running Windows.
To successfully install `numba`, try one of the following solutions:

1. Make sure your Python version is supported by *RNAlysis*. Currently *RNAlysis* supports Python versions 3.7-3.10, but version 3.8 usually works best.
2. Try `creating a new Python environment <https://towardsdatascience.com/virtual-environments-104c62d48c54?gi=40d0a7444922>`_ for *RNAlysis*, or `re-installing Python <http://docs.python-guide.org/en/latest/starting/installation/>`_.
3. Instead of installing Python normally, try installing the `Anaconda or Miniconda <https://www.edureka.co/blog/python-anaconda-tutorial/>`_ Python distribution.
You can then install `numba` via anaconda and *RNAlysis* via pip, like so:
   ::

        conda install numba
        pip install RNAlysis[all]



I am trying to open *RNAlysis* for the first time and it takes a long time to load/doesn't load at all. What should I do?
===========================================================================================================================
When you run *RNAlysis* for the first time, it will attempt to download the tutorial videos from the official *RNAlysis* repository.
This might take a long time if you have a slow internet connection, and might not work at all if your computer is not connected to the internet.
To solve this problem, try one of the following solutions:

1. Connect your computer to the internet the first time you open *RNAlysis*, to allow it to download the tutorial videos.
2. If the automatic download of the videos doesn't work, you can download them manually from the `following link<https://github.com/GuyTeichman/RNAlysis/tree/master/rnalysis/gui/videos>`_, and copy them to the *RNAlysis* appdata directory.
On a Windows computer, it would look something like this: C:\Users\<user>\AppData\Roaming\RNAlysis\RNAlysis\videos
On a Linux computer, it would look something like this:  /home/<user>/.local/share/RNAlysis/videos
On a MacOS computer, it would look something like this: /Users/<user>>/Library/Application Support/RNAlyisis/videos


*************
Usage
*************

