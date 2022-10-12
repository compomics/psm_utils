##########
Quickstart
##########

Peptidoform
###########

:py:class:`~psm_utils.peptidoform.Peptidoform` accepts peptidoforms (combination
of peptide, modifications, and — optionally — charge state) in ProForma 2.0
notation and supports several peptide-related operations, e.g.:

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
####################

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
   >>> psm.universal_spectrum_identifier()
   mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555:VLHPLEGAVVIIFK/2


PSMList and psm_utils.io
########################

The :py:mod:`psm_utils.io` subpackage contains readers and writers for various
PSM file formats (see :ref:`Supported file formats`). Each reader parses the
specific PSM file format into a unified :py:class:`~psm_utils.psm_list.PSMList`
object, with peptidoforms parsed into the ProForma notation:

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

Take a look at the :doc:`Python API Reference <api/psm_utils>` for details, more examples, and additional
information on the supported file formats.


