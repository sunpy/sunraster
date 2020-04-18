============
Installation
============

``sunraster`` requires Python 3.7+, SunPy 1.0.0+, ndcube 1.3.0+, astropy and
matplotlib.

Below we will outline how to install ``sunraster``.
The stable version of ``sunraster`` should be reliable and mostly bug free.
(If you do find a bug in the stable version please let us know.
See the :ref:`reporting_bugs` section of this guide to learn how.)
However, if users would like to get new features at soon as possible or
help to develop ``sunraster``, they will have to install the development version.
The cost of using the development version is a less stable code-base. We will outline
how to perform both types of install below.

.. _stable_install:

Installing the Stable Version
-----------------------------

There are two options for installing the stable version of ``sunraster``. The first is
via the anaconda distribution using the conda-forge channel.  For more
information on installing the anaconda distribution, see the `anaconda
website`_.

.. code-block:: console

  conda install --channel conda-forge sunraster

To update ``sunraster`` do:

.. code-block:: console

  conda update sunraster

The second option for installing the stable version of ``sunraster`` is via
pip.

.. code-block:: console

  pip install sunraster

Then to update ``sunraster`` do:

.. code-block:: console

  pip install sunraster --upgrade

.. _dev_install:

Installing the Development Version
----------------------------------

This section outlines how to install the development version of
``sunraster``. The two primary packages on which ``sunraster`` relies are `ndcube`_
and `SunPy`_. Both of these have stable released versions that work
with ``sunraster``. However, some developers may want to use the latest
updates of these packages in their work with ``sunraster``. Below we will first
outline how to install ``sunraster`` with its stable dependencies, and then
with the development versions of ``ndcube`` and ``sunpy``.

To install these packages we will use a combination of conda, conda
environments, pip and git. We will assume these are all installed on
your current system. If you do not have anaconda installed, see the
`anaconda website`_ for instructions.


Stable Dependencies Install
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create Conda Environment
""""""""""""""""""""""""
The first step is to create a conda environment (let's call it
``sunraster-dev``) in which to install the development version of ``sunraster``.
This will allow you to keep your root environment clean of development
packages.  From the command line, type:

.. code-block:: console

  conda config --append channels conda-forge
  conda create -n sunraster-dev sunpy hypothesis pytest-mock pip sphinx coverage ipython jupyter ndcube

The first line opens a conda channel so that ``sunraster`` and its
dependencies can be installed. The second line creates the
``sunraster-dev`` conda environment with a list of dependencies. Next,
you must activate that environment, i.e. switch into it.  Windows
users should type:

.. code-block:: console

  activate sunraster-dev

whereas Linux and MacOS users should type:

.. code-block:: console

  conda activate sunraster-dev

Clone ``sunraster`` Repository
""""""""""""""""""""""""""""""

The second step is to clone the `sunraster repository`_ from `GitHub`_ into
a directory.  Let's call it ``sunraster-git``. From the directory in which you
want ``sunraster-git`` to reside, type:

.. code-block:: console

  git clone https://github.com/sunpy/sunraster.git sunraster-git

Install ``sunraster``
"""""""""""""""""""""
Finally, we can install the ``sunraster`` development version:

.. code-block:: console

  cd sunraster-git
  pip install -e .

You should now be ready to use ``sunraster``. To check it's installed, open
an Python/IPython/Jupyter Notebook session from any directory and try:

.. code-block:: python

  >>> import sunraster

To make sure you have the latest updates, regularly do

.. code-block:: console

  git pull origin master

Development Dependencies Install
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create Conda Environment
""""""""""""""""""""""""
The first step is to create a conda environment (let's call it
``sunraster-dev``) in which to install the development version of ``sunraster``.
This will allow you to keep your root environment clean of development
packages.  From the command line, type:

.. code-block:: console

  conda config --append channels conda-forge
  conda create -n sunraster-dev sunpy hypothesis pytest-mock pip sphinx coverage ipython jupyter ndcube

The first line opens a conda channel so that ``sunraster`` and its
dependencies can be installed. The second line creates the
``sunraster-dev`` conda environment with a list of dependencies. Next,
you must activate that environment, i.e. switch into it. Windows
users should type:

.. code-block:: console

  activate sunraster-dev

while Linux and MacOS users should type:

.. code-block:: console

  conda activate sunraster-dev

Remove Stable Versions of ``sunpy`` and ``ndcube``
""""""""""""""""""""""""""""""""""""""""""""""""""

We installed the stable versions of SunPy and ndcube above in
order to get get all their dependencies. Now that is done, the second
step is to remove the stable versions of SunPy and ndcube, leaving the
dependencies intact.

.. note::

    **CAUTION**: Make sure you are in (have activated) the ``sunraster-dev``
    conda environment otherwise the next step will remove SunPy and ndcube
    from the wrong conda environment. From the command line in any
    directory, type:

.. code-block:: console

  conda remove ndcube
  conda remove sunpy

.. _clone_repos:

Clone Development Versions of ``sunpy``, ``ndcube`` and ``sunraster``
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Let's make a directory and then clone (download) the
development versions of `SunPy,`_ `ndcube,`_ and `sunraster,`_ from
`GitHub`_ into subdirectories.  Let's call them ``sunpy-git``,
``ndcube-git``, ``sunraster-git``.  On the command line from the
directory in which you want your repos to live, type:

.. code-block:: console

  mkdir github_repos
  cd github_repos
  git clone https://github.com/sunpy/sunpy.git sunpy-git
  git clone https://github.com/sunpy/ndcube.git ndcube-git
  git clone https://github.com/sunpy/sunraster.git sunraster-git

If you already have these repos cloned, make sure they are up-to-date
but by pulling the latest version of the master branches. For example,
for ``sunpy``, do:

.. code-block:: console

  cd ~/github_repos/sunpy-git
  git pull origin master

assuming that ``origin`` is the remote pointing to the main ``sunpy``
repo, i.e. https://github.com/sunpy/sunpy.git. The same should be done
for ``ndcube`` and ``sunraster``. To determine the correct remote name , ``cd``
into the repo's directory and do

.. code-block:: console

  git remote -v

Install the Development Versions of ``sunpy``, ``ndcube`` and ``sunraster``
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. code-block:: console

  cd ~/github_repos/sunpy-git
  pip install -e .
  cd ~/github_repos/ndcube-git
  pip install -e .
  cd ~/github_repos/sunraster-got
  pip install -e .

You should now be ready to use ``sunraster``. To check it's installed, open
an Python/IPython/Jupyter Notebook session from any directory and try:

.. code-block:: python

  import sunraster

N.B. To ensure you continue to have the latest version of ``sunraster``, be
sure to regularly update the ``sunpy``, ``ndcube`` and ``sunraster`` git repos as
discussed at the end of :ref:`clone_repos`.

.. _ndcube: http://docs.sunpy.org/projects/ndcube/en/stable/
.. _SunPy: http://sunpy.org
.. _anaconda website: https://docs.anaconda.com/anaconda/install.html
.. _sunraster repository: https://github.com/sunpy/sunraster
.. _GitHub: https://github.com/
.. _SunPy,: https://github.com/sunpy/sunpy
.. _ndcube,: https://github.com/sunpy/ndcube
.. _sunraster,: https://github.com/sunpy/sunraster
