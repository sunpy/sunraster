*********
sunraster
*********

sunraster is an open-source Python library for something.

Installation
============

There is no easy way to install this (yet).

Developing
==========

If you want to develop sunraster you will need to install from GitHub.
The best way to do this is to create a new conda environment and install the git version of sunraster in it:

.. code:: bash

    $ conda config --append channels conda-forge
    $ conda create -n sunraster-dev sunpy
    $ conda activate sunraster-dev
    $ git clone https://github.com/sunpy/sunraster.git sunraster-git
    $ cd sunraster-git
    $ pip install -e .[all,dev]

For detailed installation instructions, see the `Newcomers' guide`_ in the SunPy docs.

Getting Help
============

For more information or to ask questions about sunraster, check out:

-  `sunraster Documentation`_
-  `SunPy Matrix Channel`_
-  `SunPy Mailing List`_

.. _SunPy Documentation: https://docs.sunpy.org/en/stable/
.. _SunPy Matrix Channel: https://chat.openastronomy.org/#/room/#sunpy:openastronomy.org
.. _SunPy Mailing List: https://groups.google.com/forum/#!forum/sunpy

Contributing
============

If you would like to get involved, start by joining the `SunPy mailing list`_ and check out the `Developers Guide`_ section of the SunPy docs.
Stop by our chat room `#sunpy:openastronomy.org`_ if you have any questions.
Help is always welcome so let us know what you like to work on, or check out the `issues page`_ for the list of known outstanding items.

For more information on contributing to sunraster, please read SunPy's `Newcomers' guide`_.

.. _SunPy mailing list: https://groups.google.com/forum/#!forum/sunpy
.. _Developers Guide: https://docs.sunpy.org/en/latest/dev_guide/index.html
.. _`#sunpy:openastronomy.org`: https://chat.openastronomy.org/#/room/#sunpy:openastronomy.org
.. _issues page: https://github.com/sunpy/sunraster/issues
.. _Newcomers' guide: https://docs.sunpy.org/en/latest/dev_guide/newcomers.html

Code of Conduct
===============

When you are interacting with the SunPy community you are asked to follow our `Code of Conduct`_.

.. _Code of Conduct: https://docs.sunpy.org/en/latest/code_of_conduct.html
