##########
Quickstart
##########

Peptidoform
###########

:py:class:`~psm_utils.peptidoform.Peptidoform` accepts peptidoforms (combination
of peptide, modifications, and — optionally — charge state) in `ProForma 2.0
notation <https://github.com/HUPO-PSI/ProForma/>`_ and supports several peptide-related
operations, e.g.:

.. code-block:: python

   >>> from psm_utils import Peptidoform, PSM, PSMList
   >>> peptidoform = Peptidoform("ACDEK/2")
   >>> peptidoform.theoretical_mass
   564.2213546837


.. code-block:: python

   >>> peptidoform.composition
   Composition({'H': 36, 'C': 21, 'O': 10, 'N': 6, 'S': 1})


.. code-block:: python

   >>> peptidoform.sequential_composition
   [Composition({'H': 1}),
   Composition({'H': 5, 'C': 3, 'O': 1, 'N': 1}),
   Composition({'H': 5, 'C': 3, 'S': 1, 'O': 1, 'N': 1}),
   Composition({'H': 5, 'C': 4, 'O': 3, 'N': 1}),
   Composition({'H': 7, 'C': 5, 'O': 3, 'N': 1}),
   Composition({'H': 12, 'C': 6, 'N': 2, 'O': 1}),
   Composition({'H': 1, 'O': 1})]


PSM
###

:py:class:`~psm_utils.psm.PSM` links a
:py:class:`~psm_utils.peptidoform.Peptidoform` to a specific spectrum where it
was (putatively) identified. A :py:class:`~psm_utils.psm.PSM`
therefore contains the peptidoform, spectrum (meta)data, and peptide-spectrum
match information:

.. code-block:: python

   >>> psm = PSM(
   ...     peptidoform=Peptidoform("VLHPLEGAVVIIFK/2"),
   ...     spectrum_id=17555,
   ...     run="Adult_Frontalcortex_bRP_Elite_85_f09",
   ...     collection="PXD000561",
   ...     is_decoy=False,
   ...     precursor_mz=767.9714,
   ... )
   >>> psm.get_usi()
   mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555:VLHPLEGAVVIIFK/2

The spectrum can be retrieved by the USI through the ProteomeXchange USI aggregator:
http://proteomecentral.proteomexchange.org/usi/?usi=mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555:VLHPLEGAVVIIFK/2


PSMList and psm_utils.io
########################

The :py:mod:`psm_utils.io` subpackage contains readers and writers for various
PSM file formats (see :ref:`Supported file formats`). Each reader parses the
specific PSM file format into a unified :py:class:`~psm_utils.psm_list.PSMList`
object, with peptidoforms parsed into the ProForma notation. Use the high-level
:py:func:`psm_utils.io.read_file`, :py:func:`psm_utils.io.write_file`, and
:py:func:`psm_utils.io.convert` functions to easily read, write, and convert
PSM files:

.. code-block:: python

   >>> from psm_utils.io import read_file
   >>> psm_list = read_file("data/QExHF04054_tandem.idXML", filetype="idxml")
   >>> psm_list[0]
   PSM(
      peptidoform=Peptidoform('QSGD[Ammonium]E[Ammonium]SYC[Carbamidomethyl]E[Ammonium]R/2'),
      spectrum_id='controllerType=0 controllerNumber=1 scan=4941',
      run=None,
      collection=None,
      spectrum=None,
      is_decoy=True,
      score=17.1,
      precursor_mz=624.252254215645,
      retention_time=1197.74208,
      protein_list=['sP06800'],
      source='idXML',
      provenance_data=None,
      metadata={
         'idxml:score_type': 'XTandem',
         'idxml:higher_score_better': 'True',
         'idxml:significance_threshold': '0.0'
      },
      rescoring_features=None
   )


Alternatively, the more low-level file format-specific reader and writer classes can be
used. Each reader has a :py:meth:`read_file` function:

>>> from psm_utils.io.mzid import MzidReader
>>> psm_list = MzidReader("psms.mzid").read_file()
>>> psm_list[0].peptidoform
Peptidoform('GLTEGLHGFHVHEFGDNTAGC[Carbamidomethyl]TSAGPHFNPLSR/4')


And all readers support iteration over PSMs:

>>> for psm in MzidReader("psms.mzid"):
...     print(psm.peptidoform.proforma)
ACDEK
AC[Carbamidomethyl]DEFGR
[Acetyl]-AC[Carbamidomethyl]DEFGHIK
[...]


Similarly, writers can write single PSMs to a file:

>>> from psm_utils.io.tsv import TSVWriter
>>> with TSVWriter("psm_list.tsv", example_psm=psm_list[0]) as writer:
...     writer.write_psm(psm_list[0])


And writers can write entire PSM lists at once:

>>> with TSVWriter("psm_list.tsv", example_psm=psm_list[0]) as writer:
...     writer.write_file(psm_list)


Take a look at the :doc:`Python API Reference <api/psm_utils>` for details, more
examples, and additional information on the supported file formats.



Handling peptide modifications
##############################


Supported notations
*******************

:py:class:`~psm_utils.peptidoform.Peptidoform` accepts all supported
`ProForma 2.0 <https://github.com/HUPO-PSI/ProForma/>`_ modification types and
notations, through the :py:mod:`pyteomics.proforma` module. However, for some
functionality, such as the :py:attr:`~psm_utils.peptidoform.Peptidoform.composition` and
:py:attr:`~psm_utils.peptidoform.Peptidoform.mass` properties, the modification
composition and mass, respectively, should be resolvable. This can be achieved in
multiple ways:

Using a controlled vocabulary identifier or name, such as PSI-MOD or Unimod:

>>> Peptidoform("AC[UNIMOD:4]DEK").theoretical_mass
621.24282637892

>>> Peptidoform("AC[U:4]DEK").theoretical_mass
621.24282637892

>>> Peptidoform("AC[U:Carbamidomethyl]DEK").theoretical_mass
621.24282637892


Using a molecular formula or mass shift:

>>> Peptidoform("AC[Formula:H3C2NO]DEK/2").theoretical_mass
621.24282637892

>>> Peptidoform("AC[+57.021464]DEK/2").theoretical_mass
621.24282637892


A drawback of using the mass shift is that the composition is not be resolvable:

>>> Peptidoform("AC[+57.021464]DEK/2").composition
[...]
ModificationException: Cannot resolve composition for modification 57.021464.


Renaming modifications
**********************

Often search engines use specific, arbitrary names for modifications. In that case,
properties such as their mass or composition will not be resolvable.

>>> from psm_utils.io import read_file
>>> psm_list = read_file("msms.txt")
>>> psm_list["peptidoform"]
array([Peptidoform('AAAAAAALQAK/2'),
       Peptidoform('[ac]-AAAAAEQQQFYLLLGNLLSPDNVVR/3'),
       Peptidoform('[ac]-AAAAAEQQQFYLLLGNLLSPDNVVRK/3'), ...,
       Peptidoform('YYYLPLVSN[de]PK/2'),
       Peptidoform('YYYLTNVERLEELESDLK/3'), Peptidoform('YYYNGFYLLWI/3')],
      dtype=object)

To address this issue, modifications can be renamed:

>>> psm_list.rename_modifications({
    "ac": "U:Acetylation",
    "ox": "U:Oxidation",
    "de": "U:Deamidation",
    "gl": "U:Gln->pyro-Glu",
})
>>> psm_list["peptidoform"]
array([Peptidoform('AAAAAAALQAK/2'),
       Peptidoform('[UNIMOD:Acetylation]-AAAAAEQQQFYLLLGNLLSPDNVVR/3'),
       Peptidoform('[UNIMOD:Acetylation]-AAAAAEQQQFYLLLGNLLSPDNVVRK/3'),
       ..., Peptidoform('YYYLPLVSN[UNIMOD:Deamidation]PK/2'),
       Peptidoform('YYYLTNVERLEELESDLK/3'), Peptidoform('YYYNGFYLLWI/3')],
      dtype=object)


Handling fixed modifications
****************************

Additionally, fixed modifications that are not already part of the search engine output
can be added and applied across the sequence:

>>> psm_list[19].peptidoform
Peptidoform('AAAPAPEEEMDECEQALAAEPK/2')

>>> psm_list.add_fixed_modifications([("Carbamidomethyl", ["C"])])
>>> psm_list[19].peptidoform
Peptidoform('<[Carbamidomethyl]@C>AAAPAPEEEMDECEQALAAEPK/2')

>>> psm_list.apply_fixed_modifications()
>>> psm_list[19].peptidoform
Peptidoform('AAAPAPEEEMDEC[Carbamidomethyl]EQALAAEPK/2')
