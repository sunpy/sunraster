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
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create Conda Environment
""""""""""""""""""""""""
The first step is to create a conda environment (let's call it
``irispy-dev``) in which to install the development version of IRISpy.
This will allow you to keep your root environment clean of development
packages.  From the command line, type:

.. code-block:: console

		conda config --append channels conda-forge
		conda create -n irispy-dev sunpy hypothesis pytest-mock pip sphinx coverage ipython jupyter ndcube

The first line opens a conda channel so that IRISpy and its
dependencies can be installed. The second line creates the
``irispy-dev`` conda environment with a list of dependencies. Next,
you must activate that environment, i.e. switch into it.  Windows
users should type:

.. code-block:: console

		activate irispy-dev

while Linux and MacOS users should type:

.. code-block:: console

		source activate irispy-dev

Clone IRISpy Repository
"""""""""""""""""""""""

The second step is to clone the IRISpy repository from `GitHub`_ into
a directory called ``irispy-git``. From the directory in which you
want ``irispy-git`` to reside, type:

.. code-block:: console

		git clone https://github.com/sunpy/irispy.git irispy-git

Install IRISpy
""""""""""""""
Finally, we can install the IRISpy development version.

.. code-block:: console

		cd irispy-git
		pip install -e .

You should now be ready to use IRISpy. To check it's installed, open
an Python/IPython/Jupyter Notebook session from any directory and try:

.. code-block:: python

		>>> import irispy

To make sure you have the latest updates, regularly do

.. code-block:: console

		git pull origin master

Development Dependencies Install
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create Conda Environment
""""""""""""""""""""""""
The first step is to create a conda environment (let's call it
``irispy-dev``) in which to install the development version of IRISpy.
This will allow you to keep your root environment clean of development
packages.  From the command line, type:

.. code-block:: console

		conda config --append channels conda-forge
		conda create -n irispy-dev sunpy hypothesis pytest-mock pip sphinx coverage ipython jupyter ndcube

The first line opens a conda channel so that IRISpy and its
dependencies can be installed. The second line creates the
``irispy-dev`` conda environment with a list of dependencies. Next,
you must activate that environment, i.e. switch into it.  Windows
users should type:

.. code-block:: console

		activate irispy-dev

while Linux and MacOS users should type:

.. code-block:: console

		source activate irispy-dev

Remove Stable Versions of SunPy and ndcube
""""""""""""""""""""""""""""""""""""""""""

We installed the stable versions of SunPy and ndcube above in
order to get get all its dependencies. Now that is done, the second
step is to remove the stable versions of SunPy and ndcube, leaving the
dependencies intact.
CAUTION: Make sure you are in (have activated) the ``irispy-dev``
conda environment otherwise the next step will remove SunPy and ndcube
from the wrong conda environment. From the command line in any
directory, type:

.. code-block:: console

		conda remove ndcube
		conda remove sunpy

.. _clone_repos:

Clone Development Versions of SunPy, ndcube and IRISpy
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Let's make a directory and then clone (download) the
development versions of the SunPy, ndcube and IRISpy from GitHub into
subdirectories.  Let's call them ``sunpy-git``, ``ndcube-git``,
``irispy-git``.  On the command line from the directory in which you
want your repos to live, type:

.. code-block:: console

		mkdir github_repos
		cd github_repos
		git clone https://github.com/sunpy/sunpy.git sunpy-git
		git clone https://github.com/sunpy/ndcube.git ndcube-git
		git clone https://github.com/sunpy/irispy.git irispy-git

If you already have these repos cloned, make sure they are up-to-date
but by pulling the latest version of the master branches. For example,
for sunpy, do:

.. code-block:: console

		cd ~/github_repos/sunpy-git
		git pull origin master

assuming that ``origin`` is the remote pointing to the main sunpy
repo, i.e. https://github.com/sunpy/sunpy.git. The same should be done
for ndcube and irispy. To determine the correct remote name , ``cd``
into the repo's directory and do

.. code-block:: console

		git remote -v

Install the Development Versions of SunPy, ndcube and IRISpy
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. code-block:: console

		cd ~/github_repos/sunpy-git
		pip install -e .
		cd ~/github_repos/ndcube-git
		pip install -e .
		cd ~/github_repos/irispy-got
		pip install -e .

You should now be ready to use IRISpy. To check it's installed, open
an Python/IPython/Jupyter Notebook session from any directory and try:

.. code-block:: python

		import irispy

N.B. To ensure you continue to have the latest version of IRISpy, be
sure to regularly update the sunpy, ndcube and irispy git repos as
discussed at the end of :ref:`clone_repos`.

.. _ndcube: https://github.com/sunpy/ndcube
.. _SunPy: https://github.com/sunpy/sunpy
.. _anaconda website: https://docs.anaconda.com/anaconda/install.html
.. _GitHub: https://github.com/
.. _ndcube GitHub repository: https://github.com/sunpy/ndcube
