=========
sunraster
=========

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


License
=======

This project is Copyright (c) The SunPy Community and licensed under
the terms of the BSD 2-Clause license. This package is based upon
the `Openastronomy packaging guide <https://github.com/OpenAstronomy/packaging-guide>`_
which is licensed under the BSD 3-clause licence. See the licenses folder for
more information.

Usage of Generative AI
======================

We expect authentic engagement in our community.
Be wary of posting output from Large Language Models or similar generative AI as comments on GitHub or any other platform, as such comments tend to be formulaic and low quality content.
If you use generative AI tools as an aid in developing code or documentation changes, ensure that you fully understand the proposed changes and can explain why they are the correct approach and an improvement to the current state.

Contributing
============

We love contributions! sunraster is open source,
built on open source, and we'd love to have you hang out in our community.

**Imposter syndrome disclaimer**: We want your help. No, really.

There may be a little voice inside your head that is telling you that you're not
ready to be an open source contributor; that your skills aren't nearly good
enough to contribute. What could you possibly offer a project like this one?

We assure you - the little voice in your head is wrong. If you can write code at
all, you can contribute code to open source. Contributing to open source
projects is a fantastic way to advance one's coding skills. Writing perfect code
isn't the measure of a good developer (that would disqualify all of us!); it's
trying to create something, making mistakes, and learning from those
mistakes. That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can
help out by writing documentation, tests, or even giving feedback about the
project (and yes - that includes giving feedback about the contribution
process). Some of these contributions may be the most valuable to the project as
a whole, because you're coming to the project with fresh eyes, so you can see
the errors and assumptions that seasoned contributors have glossed over.

For more information on contributing to sunraster, please read SunPy's `Newcomers' guide`_.

.. _SunPy mailing list: https://groups.google.com/forum/#!forum/sunpy
.. _Developers Guide: https://docs.sunpy.org/en/latest/dev_guide/index.html
.. _`#sunpy:openastronomy.org`: https://chat.openastronomy.org/#/room/#sunpy:openastronomy.org
.. _issues page: https://github.com/sunpy/sunraster/issues
.. _Newcomers' guide: https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html


Note: This disclaimer was originally written by
`Adrienne Lowe <https://github.com/adriennefriend>`_ for a
`PyCon talk <https://www.youtube.com/watch?v=6Uj746j9Heo>`_, and was adapted by
sunraster based on its use in the README file for the
`MetPy project <https://github.com/Unidata/MetPy>`_.
