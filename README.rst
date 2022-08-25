#########
psm_utils
#########

Common utilities for parsing and handling peptide-spectrum matches and search
engine results in Python.


.. image:: https://img.shields.io/github/v/release/compomics/psm_utils?sort=semver&style=flat-square
   :alt: GitHub release
   :target: https://github.com/compomics/psm_utils/releases

.. image:: https://img.shields.io/pypi/v/psm-utils?style=flat-square
   :alt: PyPI
   :target: https://pypi.org/project/psm-utils/

.. image:: https://img.shields.io/github/workflow/status/compomics/psm_utils/Test?label=tests&style=flat-square
   :alt: GitHub Actions tests status
   :target: https://github.com/compomics/psm_utils/actions/workflows/test.yml

.. image:: https://img.shields.io/github/workflow/status/compomics/psm_utils/Publish?label=build&style=flat-square
   :alt: GitHub Actions build status
   :target: https://github.com/compomics/psm_utils/actions/workflows/publish.yml

.. image:: https://img.shields.io/codecov/c/github/compomics/psm_utils?style=flat-square
   :alt: Codecov
   :target: https://app.codecov.io/gh/compomics/psm_utils

.. image:: https://img.shields.io/github/license/compomics/psm_utils.svg?style=flat-square
   :alt: GitHub
   :target: https://www.apache.org/licenses/LICENSE-2.0

.. image:: https://img.shields.io/twitter/follow/CompOmics?style=flat-square
   :alt: Twitter
   :target: https://twitter.com/compomics



About
#####

Introduction
************

psm_utils is a Python package with utilities for parsing and
handling peptide-spectrum matches (PSMs) and proteomics search engine results.
It is mainly developed to be used in Python packages developed at
`CompOmics <https://www.compomics.com>`_, such as
`MS²PIP <https://github.com/compomics/ms2pip_c>`_,
`DeepLC <https://github.com/compomics/deeplc>`_, and
`MS²Rescore <https://github.com/compomics/ms2rescore>`_,
but can be useful to anyone dealing with PSMs and PSM files. Moreover, it
provides an easy-to-use CLI and web server to convert search engine results from
one PSM file format into another.


Goals and non-goals
*******************
- To provide an easy-to-use Python API for **handling PSMs**.
- To provide a unified Python API to the plethora of **proteomics search engine
  output formats** that are in existence.
- To follow **community standards**: psm_utils pragmatically adheres to the
  standards developed by the
  `HUPO Proteomics Standards Initiative <http://psidev.info>`_, such as
  `ProForma 2.0 <https://psidev.info/proforma>`_ , the
  `Universal Spectrum Identifier <https://psidev.info/usi>`_, and
  `mzIdentML <https://psidev.info/mzidentml>`_
- To be **open and dynamic**: psm_utils is fully open source, under the
  permissive Apache 2.0 license. New reader and writer modules can easily be
  added, and we welcome everyone to contribute to the project. See
  :doc:`contributing` for more information.
- **NOT to reinvent the wheel**: Instead, psm_utils heavily makes
  use of packages such as `pyteomics <http://pyteomics.readthedocs.io/>`_ and
  `psims <https://github.com/mobiusklein/psims>`_ that have existing
  functionality for reading and/or writing PSM files. ``psm_utils.io``
  provides a unified, higher level Python API build on top of these packages.


Supported file formats
**********************

 ===================================================================================================================== =============== ===============
  File format                                                                                                           Read support    Write support
 ===================================================================================================================== =============== ===============
  `OpenMS idXML <https://www.openms.de/>`_                                                                              ✅              ❌
  `MaxQuant msms.txt <https://www.maxquant.org/>`_                                                                      ✅              ❌
  `Peptide Record <https://psm-utils.readthedocs.io/en/latest/api/psm_utils.io/#module-psm_utils.io.peptide_record>`_   ✅              ✅
  `Percolator tab <https://github.com/percolator/percolator/wiki/Interface>`_                                           ✅              ✅
  `TSV <https://psm-utils.readthedocs.io/en/latest/api/psm_utils.io/#module-psm_utils.io.tsv>`_                         ✅              ✅
  `X!Tandem XML <https://www.thegpm.org/tandem/>`_                                                                      ✅              ❌
 ===================================================================================================================== =============== ===============



Installation
############

Install with pip:

.. code-block:: sh

    pip install psm-utils


Note: In the PyPI package name, a hyphen is used instead of an underscore, as
per `PEP8 convention <https://peps.python.org/pep-0008/#package-and-module-names>`_.



Full documentation
##################

The full documentation, including a quickstart guide and Python API reference
is available on `psm_utils.readthedocs.io <https://psm-utils.readthedocs.io>`_.
