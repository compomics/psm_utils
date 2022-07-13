#########
psm_utils
#########

Common utilities for parsing and handling peptide-spectrum matches and search
engine results in Python.


.. image:: https://flat.badgen.net/github/release/compomics/psm_utils
    :alt: GitHub release
    :target: https://github.com/compomics/psm_utils/releases

.. image:: https://flat.badgen.net/pypi/v/psm-utils
    :alt: PyPI
    :target: https://pypi.org/project/psm-utils/

.. image:: https://flat.badgen.net/github/checks/compomics/psm_utils/master
    :alt: GitHub Workflow Status
    :target: https://github.com/compomics/psm_utils/actions/

.. image:: https://img.shields.io/github/license/compomics/psm_utils.svg?style=flat-square
    :alt: GitHub
    :target: https://www.apache.org/licenses/LICENSE-2.0

.. image:: https://flat.badgen.net/twitter/follow/compomics?icon=twitter
    :alt: Twitter
    :target: https://twitter.com/compomics



About
#####
*TODO*


Installation
############

Install with pip:

.. code-block:: sh

    pip install psm-utils


Quickstart
##########
*TODO*

Development
###########

Local install
*************

For development, install with pip in editable mode:

.. code-block:: sh

    pip install --editable .[dev]


Or with Flit:

.. code-block:: sh

    flit install --symlink


Unit tests
**********

Run tests with pytest:

.. code-block:: sh

    pytest ./tests


Documentation
*************

To work on the documentation and get a live preview, install the requirements
and run ``sphinx-autobuild``:

.. code-block:: sh

    pip install .[doc]
    sphinx-autobuild  --watch ./psm_utils ./docs/source/ ./docs/_build/html/