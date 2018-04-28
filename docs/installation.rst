============
Installation
============

IRISpy requires Python 3.5+, SunPy 0.9+, ndcube 1.0.1+, astropy and
matplotlib.

.. warning::
    
    IRISPy is still under heavy development and has not yet seen its first
    release.  The API can change at any time and so should not be
    relied upon.  However, we are striving towards releasing a stable
    version.  If you would like to help by providing user feedback,
    reporting bugs, or contributing code, install the development
    version as outlined below.

.. _dev_install:

Installing the Development Version
----------------------------------

This section outlines how to install the development version of
IRISpy. The two primary packages on which IRISpy relies are ndcube and
SunPy. Both of these have stable released versions that work with
IRISpy. However, some developers may want to use the latest updates of
these packages in their work on IRISpy. Below we will first outline
how to install IRISpy with its stable dependencies, and then with the
development versions of ndcube and SunPy.

To install these packages we will use a combination of conda, conda
environments, pip and git. We will assume these are all installed on
your current system. If you do not have anaconda installed, see the
`anaconda website`_ for instructions.


Stable Dependencies Install
---------------------------

The first step is to create a conda environment, called irispy-dev in
which to install the development version of IRISpy.  This will allow
you to keep your root environment clean of development packages.  From
the command line, type:

.. code-block:: console

		conda config --append channels conda-forge
		conda create -n irispy-dev sunpy hypothesis pytest-mock pip sphinx coverage ipython jupyter ndcube

The first line opens a conda channel so that IRISpy and its
dependencies can be installed. The second line creates the irispy-dev
conda environment with a list of dependencies. Next, you must activate
that environment, i.e. switch into it.  Windows users should type:

.. code-block:: console

		activate irispy-dev

while Linux and MacOS users should type:

.. code-block:: console

		source activate irispy-dev

The second step is to clone the IRISpy repository from `GitHub`_ into
a directory called ``irispy-git``. From the directory you want
``irispy-git`` to reside, do:

.. code-block:: console

		$ git clone https://github.com/sunpy/irispy.git irispy-git

Finally, we can install the IRISpy development version.

.. code-block:: console

		$ cd irispy-git
		$ pip install -e .

You should now be ready to use IRISpy. To check it's installed, open
an Python/IPython/Jupyter Notebook session from any directory and try:
.. code-block:: python

		>>> import irispy

To make sure you have the latest updates, regularly do

.. code-block:: console

		git pull origin master

.. _ndcube: https://github.com/sunpy/ndcube
.. _SunPy: https://github.com/sunpy/sunpy
.. _anaconda website: https://docs.anaconda.com/anaconda/install.html
.. _GitHub: https://github.com/
.. _ndcube GitHub repository: https://github.com/sunpy/ndcube
