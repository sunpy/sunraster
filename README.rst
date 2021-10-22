*********
sunraster
*********

|Latest Version| |codecov| |matrix| |DOI| |Powered by NumFOCUS| |Powered by SunPy|

.. |Latest Version| image:: https://img.shields.io/pypi/v/sunraster.svg
   :target: https://pypi.python.org/pypi/sunraster/
.. |matrix| image:: https://img.shields.io/matrix/sunpy:openastronomy.org.svg?colorB=%23FE7900&label=Chat&logo=matrix&server_fqdn=openastronomy.modular.im
   :target: https://openastronomy.element.io/#/room/#sunpy:openastronomy.org
.. |codecov| image:: https://codecov.io/gh/sunpy/sunraster/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/sunpy/sunraster
.. |DOI| image:: https://zenodo.org/badge/2165383.svg
   :target: https://zenodo.org/badge/latestdoi/2165383
.. |Powered by NumFOCUS| image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
   :target: https://numfocus.org
.. |Powered by SunPy| image:: https://img.shields.io/badge/powered%20by-SunPy-orange.svg?style=flat
   :target: https://www.sunpy.org

``sunraster`` is an open-source Python library that provides the tools to read in and analyze spectrogram data.

Installation
============

An easy way to install ``sunraster`` is to do so with the anaconda distribution using the conda-forge channel, with the following command at the terminal:

.. code-block:: console

    conda install --channel conda-forge sunraster

Another equally easy way to install ``sunraster`` is with pip:

.. code-block:: console

    pip install sunraster

Developing
==========

If you want to develop ``sunraster`` you will need to install from GitHub.
We suggest you fork ``sunraster`` so you can work on it.
The best way to do this is to create a new python virtual environment (conda/pipenv or others) and then install the git version of ``sunraster``:

.. code:: bash

    $ git clone https://github.com/<your username>/sunraster.git
    $ cd sunraster
    $ pip install -e .\[dev\]


For detailed installation instructions (aimed at installing ``sunpy``), see the `Newcomers' guide`_ in the sunpy docs.

Getting help
============

For more information or to ask questions about ``sunraster``, check out:

-  `sunraster Documentation`_
-  `sunpy Matrix Channel`_
-  `sunpy Mailing List`_

.. _sunraster Documentation: https://docs.sunpy.org/projects/sunraster/en/latest/
.. _sunpy Matrix Channel: https://chat.openastronomy.org/#/room/#sunpy:openastronomy.org
.. _sunpy Mailing List: https://groups.google.com/forum/#!forum/sunpy

Contributing
============

If you would like to get involved, start by joining the `SunPy mailing list`_ and check out the `Developers Guide`_ section of the SunPy docs.
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
