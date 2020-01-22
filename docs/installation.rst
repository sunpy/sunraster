============
Installation
============

SlitSpectrographPy requires Python 3.5+, SunPy 0.9+, ndcube 1.0.1+, astropy and
matplotlib.

.. warning::
    
    SlitSpectrographPy is still under heavy development and has not yet seen its first
    release.  The API can change at any time and so should not be
    relied upon.  However, we are striving towards releasing a stable
    version.  If you would like to help by providing user feedback,
    reporting bugs, or contributing code, install the development
    version as outlined below.

.. _dev_install:

Installing the Development Version
----------------------------------

This section outlines how to install the development version of
SlitSpectrographPy. The two primary packages on which SlitSpectrographPy relies are `ndcube`_
and `SunPy`_. Both of these have stable released versions that work
with SlitSpectrographPy. However, some developers may want to use the latest
updates of these packages in their work on SlitSpectrographPy. Below we will first
outline how to install SlitSpectrographPy with its stable dependencies, and then
with the development versions of ndcube and SunPy.

To install these packages we will use a combination of conda, conda
environments, pip and git. We will assume these are all installed on
your current system. If you do not have anaconda installed, see the
`anaconda website`_ for instructions.


Stable Dependencies Install
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create Conda Environment
""""""""""""""""""""""""
The first step is to create a conda environment (let's call it
``slitspectrographpy-dev``) in which to install the development version of SlitSpectrographPy.
This will allow you to keep your root environment clean of development
packages.  From the command line, type:

.. code-block:: console

		conda config --append channels conda-forge
		conda create -n slitspectrographpy-dev sunpy hypothesis pytest-mock pip sphinx coverage ipython jupyter ndcube

The first line opens a conda channel so that SlitSpectrographPy and its
dependencies can be installed. The second line creates the
``slitspectrographpy-dev`` conda environment with a list of dependencies. Next,
you must activate that environment, i.e. switch into it.  Windows
users should type:

.. code-block:: console

		activate slitspectrographpy-dev

while Linux and MacOS users should type:

.. code-block:: console

		source activate slitspectrographpy-dev

Clone SlitSpectrographPy Repository
"""""""""""""""""""""""

The second step is to clone the `SlitSpectrographPy repository`_ from `GitHub`_ into
a directory called ``slitspectrographpy-git``. From the directory in which you
want ``slitspectrographpy-git`` to reside, type:

.. code-block:: console

		git clone https://github.com/sunpy/slitspectrographpy.git slitspectrographpy-git

Install SlitSpectrographPy
""""""""""""""
Finally, we can install the SlitSpectrographPy development version.

.. code-block:: console

		cd slitspectrographpy-git
		pip install -e .

You should now be ready to use SlitSpectrographPy. To check it's installed, open
an Python/IPython/Jupyter Notebook session from any directory and try:

.. code-block:: python

		import slitspectrographpy

To make sure you have the latest updates, regularly do

.. code-block:: console

		git pull origin master

Development Dependencies Install
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create Conda Environment
""""""""""""""""""""""""
The first step is to create a conda environment (let's call it
``slitspectrographpy-dev``) in which to install the development version of SlitSpectrographPy.
This will allow you to keep your root environment clean of development
packages.  From the command line, type:

.. code-block:: console

		conda config --append channels conda-forge
		conda create -n slitspectrographpy-dev sunpy hypothesis pytest-mock pip sphinx coverage ipython jupyter ndcube

The first line opens a conda channel so that SlitSpectrographPy and its
dependencies can be installed. The second line creates the
``slitspectrographpy-dev`` conda environment with a list of dependencies. Next,
you must activate that environment, i.e. switch into it.  Windows
users should type:

.. code-block:: console

		activate slitspectrographpy-dev

while Linux and MacOS users should type:

.. code-block:: console

		source activate slitspectrographpy-dev

Remove Stable Versions of SunPy and ndcube
""""""""""""""""""""""""""""""""""""""""""

We installed the stable versions of SunPy and ndcube above in
order to get get all their dependencies. Now that is done, the second
step is to remove the stable versions of SunPy and ndcube, leaving the
dependencies intact.
CAUTION: Make sure you are in (have activated) the ``slitspectrographpy-dev``
conda environment otherwise the next step will remove SunPy and ndcube
from the wrong conda environment. From the command line in any
directory, type:

.. code-block:: console

		conda remove ndcube
		conda remove sunpy

.. _clone_repos:

Clone Development Versions of SunPy, ndcube and SlitSpectrographPy
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Let's make a directory and then clone (download) the
development versions of `SunPy,`_ `ndcube,`_ and `SlitSpectrographPy,`_ from
`GitHub`_ into subdirectories.  Let's call them ``sunpy-git``,
``ndcube-git``, ``slitspectrographpy-git``.  On the command line from the
directory in which you want your repos to live, type:

.. code-block:: console

		mkdir github_repos
		cd github_repos
		git clone https://github.com/sunpy/sunpy.git sunpy-git
		git clone https://github.com/sunpy/ndcube.git ndcube-git
		git clone https://github.com/sunpy/slitspectrographpy.git slitspectrographpy-git

If you already have these repos cloned, make sure they are up-to-date
but by pulling the latest version of the master branches. For example,
for sunpy, do:

.. code-block:: console

		cd ~/github_repos/sunpy-git
		git pull origin master

assuming that ``origin`` is the remote pointing to the main sunpy
repo, i.e. https://github.com/sunpy/sunpy.git. The same should be done
for ndcube and slitspectrographpy. To determine the correct remote name , ``cd``
into the repo's directory and do

.. code-block:: console

		git remote -v

Install the Development Versions of SunPy, ndcube and SlitSpectrographPy
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. code-block:: console

		cd ~/github_repos/sunpy-git
		pip install -e .
		cd ~/github_repos/ndcube-git
		pip install -e .
		cd ~/github_repos/slitspectrographpy-got
		pip install -e .

You should now be ready to use SlitSpectrographPy. To check it's installed, open
an Python/IPython/Jupyter Notebook session from any directory and try:

.. code-block:: python

		import slitspectrographpy

N.B. To ensure you continue to have the latest version of SlitSpectrographPy, be
sure to regularly update the sunpy, ndcube and slitspectrographpy git repos as
discussed at the end of :ref:`clone_repos`.

.. _ndcube: http://docs.sunpy.org/projects/ndcube/en/stable/
.. _SunPy: http://sunpy.org
.. _anaconda website: https://docs.anaconda.com/anaconda/install.html
.. _SlitSpectrographPy repository: https://github.com/sunpy/slitspectrographpy
.. _GitHub: https://github.com/
.. _SunPy,: https://github.com/sunpy/sunpy
.. _ndcube,: https://github.com/sunpy/ndcube
.. _SlitSpectrographPy,: https://github.com/sunpy/slitspectrographpy
