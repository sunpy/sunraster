``sunraster``
=============

sunraster is an open-source Python library that provides the tools to read in and analyze spectrogram data.

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
------------

An easy way to install ``sunraster`` is to do so with the anaconda distribution using the conda-forge channel, with the following command at the terminal:

.. code-block:: console

    conda install --channel conda-forge sunraster

Another equally easy way to install ``sunraster`` is with pip:

.. code-block:: console

    pip install sunraster

Developing
----------

If you want to develop ``sunraster`` you will need to install from GitHub.
We suggest you fork ``sunraster`` so you can work on it.
The best way to do this is to create a new python virtual environment (conda/pipenv or others) and then install the git version of ``sunraster``:

.. code:: bash

    $ git clone https://github.com/<your username>/sunraster.git
    $ cd sunraster
    $ pip install -e .\[dev\]


For detailed installation instructions (aimed at installing ``sunpy``), see the `Newcomers' guide`_ in the sunpy docs.

Getting help
------------

For more information or to ask questions about ``sunraster``, check out:

-  `sunraster Documentation`_
-  `sunpy Matrix Channel`_
-  `sunpy Mailing List`_

.. _sunraster Documentation: https://docs.sunpy.org/projects/sunraster/en/latest/
.. _sunpy Matrix Channel: https://chat.openastronomy.org/#/room/#sunpy:openastronomy.org
.. _sunpy Mailing List: https://groups.google.com/forum/#!forum/sunpy


Usage of Generative AI
----------------------

We expect authentic engagement in our community.
**Do not post the output from Large Language Models or similar generative AI as code, issues or comments on GitHub or any other platform.**
If you use generative AI tools as an aid in developing code or documentation changes, ensure that you fully understand the proposed changes and can explain why they are the correct approach and an improvement to the current state.
For more information see our documentation on fair and appropriate `AI usage <https://docs.sunpy.org/en/latest/dev_guide/contents/ai_usage.html>`__.

Contributing
------------

We love contributions! sunraster is open source,
built on open source, and we'd love to have you hang out in our community.

If you would like to get involved, check out the `Developers Guide`_ section of the SunPy docs.
Stop by our chat room `#sunpy:openastronomy.org`_ if you have any questions.
Help is always welcome so let us know what you like to work on, or check out the `issues page`_ for the list of known outstanding items.

For more information on contributing to SunPy, please read our `Newcomers' guide`_.

.. _Developers Guide: https://docs.sunpy.org/en/latest/dev_guide/index.html
.. _`#sunpy:openastronomy.org`: https://app.element.io/#/room/#sunpy:openastronomy.org
.. _issues page: https://github.com/sunpy/sunraster/issues
.. _Newcomers' guide: https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html

When you are interacting with the SunPy community you are asked at to follow our `code of conduct <https://sunpy.org/coc>`__.
