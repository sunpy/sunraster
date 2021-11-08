============
Installation
============

Below we will outline how to install ``sunraster``.
The stable version of ``sunraster`` is what most people will want to install.
If you do find a bug or a behavior you think is incorrect please let us know.

However, if users would like to get new features at soon as possible or help to develop ``sunraster``, they will have to install the development version.

.. _stable_install:

Installing the stable version
-----------------------------

There are two options for installing the stable version of ``sunraster``.
The first is via the anaconda distribution using the conda-forge channel.
For more information on installing the anaconda distribution, see the `anaconda website`_.

.. code-block:: console

  conda install --channel conda-forge sunraster

To update ``sunraster`` do:

.. code-block:: console

  conda update sunraster

The second option for installing the stable version of ``sunraster`` is via pip.

.. code-block:: console

  pip install sunraster

Then to update ``sunraster`` do:

.. code-block:: console

  pip install sunraster --upgrade

.. _dev_install:

Installing the development version
----------------------------------

This section outlines how to install the development version of ``sunraster``.
The two primary packages on which ``sunraster`` relies are `ndcube`_ and `sunpy`_.
Both of these have stable released versions that work with ``sunraster``.
However, some developers may want to use the latest updates of these packages in their work with ``sunraster``.

To install these packages we will use a combination of conda, conda environments, pip and git.
We will assume these are all installed on your current system.

Stable dependencies install
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create conda environment
""""""""""""""""""""""""
The first step is to create a conda environment (let's call it ``sunraster-dev``) in which to install the development version of ``sunraster``.
This will allow you to keep your root environment clean of development packages.
From the command line, type:

.. code-block:: console

  conda config --append channels conda-forge
  conda create -n sunraster-dev pip

The first line opens a conda channel so that ``sunraster`` and its dependencies can be installed.
The second line creates the ``sunraster-dev`` conda environment with a list of dependencies.
Next, you must activate that environment, i.e. switch into it.
Windows users should type:

.. code-block:: console

  activate sunraster-dev

whereas Linux and MacOS users should type:

.. code-block:: console

  conda activate sunraster-dev

Clone ``sunraster`` repository
""""""""""""""""""""""""""""""

The second step is to clone the `sunraster repository`_ from `GitHub`_ into a directory.
Let's call it ``sunraster-git``. From the directory in which you want ``sunraster-git`` to reside, type:

.. code-block:: console

  git clone https://github.com/sunpy/sunraster.git sunraster-git

If you want to develop ``sunraster``, you will need to fork the repository and clone your fork instead.

Install ``sunraster``
"""""""""""""""""""""
Finally, we can install the ``sunraster`` development version:

.. code-block:: console

  cd sunraster-git
  pip install -e .\[dev\]

You should now be ready to use ``sunraster``.
To check it's installed, open an Python/IPython/Jupyter Notebook session from any directory and try:

.. code-block:: python

  >>> import sunraster

To make sure you have the latest updates, regularly do

.. code-block:: console

  git pull origin main

.. _ndcube: https://docs.sunpy.org/projects/ndcube/en/stable/
.. _SunPy: https://sunpy.org
.. _anaconda website: https://docs.anaconda.com/anaconda/install.html
.. _sunraster repository: https://github.com/sunpy/sunraster
.. _GitHub: https://github.com/
.. _SunPy,: https://github.com/sunpy/sunpy
.. _ndcube,: https://github.com/sunpy/ndcube
