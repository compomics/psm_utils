from __future__ import annotations

from typing import Any, Dict, List, Optional, Union

from pydantic import ConfigDict, BaseModel

from psm_utils.peptidoform import Peptidoform


class PSM(BaseModel):
    """Data class representing a peptide-spectrum match (PSM)."""

    peptidoform: Union[Peptidoform, str]
    spectrum_id: Union[str]
    run: Optional[str] = None
    collection: Optional[str] = None
    spectrum: Optional[Any] = None
    is_decoy: Optional[bool] = None
    score: Optional[float] = None
    qvalue: Optional[float] = None
    pep: Optional[float] = None
    precursor_mz: Optional[float] = None
    retention_time: Optional[float] = None
    ion_mobility: Optional[float] = None
    protein_list: Optional[List[str]] = None
    rank: Optional[int] = None
    source: Optional[str] = None
    provenance_data: Optional[Dict[str, str]] = dict()
    metadata: Optional[Dict[str, str]] = dict()
    rescoring_features: Optional[Dict[str, float]] = dict()
    model_config = ConfigDict(arbitrary_types_allowed=True, coerce_numbers_to_str=True)

    def __init__(self, **data):
        """
        Data class representing a peptide-spectrum match (PSM).

        Links a :class:`~psm_utils.peptidoform.Peptidoform` to an observed spectrum
        and holds the related information. Attribute types are coerced and enforced upon
        initialization.

        Parameters
        ----------
        peptidoform : Peptidoform, str
            Peptidoform object or string in ProForma v2 notation.
        spectrum_id : str, int
            Spectrum identifier as used in spectrum file (e.g., mzML or MGF),
            usually in HUPO-PSI nativeID format (MS:1000767), e.g.,
            ``controllerType=0 controllerNumber=0 scan=423``.
        run : str, optional
            Name of the MS run. Usually the spectrum file filename without
            extension.
        collection : str, optional
            Identifier of the collection of spectrum files. Usually, the
            ProteomeXchange identifier, e.g. ``PXD028735``.
        spectrum : any, optional
            Observed spectrum. Can be freely used, for instance as a
            :py:class:`spectrum_utils.spectrum.MsmsSpectrum` object.
        is_decoy : bool, optional
            Boolean specifying if the PSM is a decoy (``True``) or target hit
            (``False``).
        score : float, optional
            Search engine score.
        qvalue : float, optional
            PSM-level q-value.
        pep : float, optional
            PSM-level posterior error probability.
        precursor_mz : float, optional
            Precursor m/z.
        retention_time : float, optional
            Precursor retention time.
        ion_mobility : float, optional
            Precursor ion mobility.
        protein_list : list[str]
            List of proteins or protein groups associated with peptide.
        rank : int
            rank of a psm
        source : str, optional
            PSM file type where PSM was stored or search engine that generated it. E.g.,
            ``mzid``, or ``X!Tandem``.
        provenance_data : dict[str, str], optional
            Freeform dict to hold data describing the PSM origin, e.g. a search
            engine-specific identifier.
        metadata : dict[str, str], optional
            More data about PSM.
        rescoring_features : dict[str, str], optional
            Dict with features that can be used for PSM rescoring.

        """
        super().__init__(**data)
        # Parse peptidoform
        if isinstance(self.peptidoform, str):
            self.peptidoform = Peptidoform(self.peptidoform)
        elif not isinstance(self.peptidoform, Peptidoform):
            raise TypeError(
                f"Peptidoform or str expected for `peptidoform`, not `{type(self.peptidoform)}`."
            )

    def __getitem__(self, item) -> any:
        return getattr(self, item)

    def __setitem__(self, item, value: any) -> None:
        setattr(self, item, value)

    @property
    def precursor_mz_error(self) -> float:
        """Difference between observed and theoretical m/z in Da."""
        theoretical_mz = self.peptidoform.theoretical_mz
        return self.precursor_mz - theoretical_mz

    def get_precursor_charge(self) -> int:
        """Precursor charge, as embedded in :py:attr:`PSM.peptidoform`."""
        return self.peptidoform.precursor_charge

    def get_usi(self, as_url=False) -> str:
        """
        Compile Universal Spectrum Identifier for :py:class:`~psm_utils.psm.PSM`.

        Parameters
        ----------
        as_url : bool, optional
            Return URL to proteomeXchange.org USI aggregator.

        Notes
        -----
        - The resulting USI will only be valid if the ``collection``, ``run``,
          and ``spectrum_id`` are defined.
        - There is no guarantee that the resulting USI is resolvable at ProteomeXchange.
          This requires that the spectrum has been fully indexed in a ProteomeXchange
          partner repository. For instance, the following USI should be resolvable:
          `mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555:VLHPLEGAVVIIFK/2 <http://proteomecentral.proteomexchange.org/usi/?usi=mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555:VLHPLEGAVVIIFK/2>`_

        """
        # TODO: Support nativeID spectrum flag
        usi = f"mzspec:{self.collection}:{self.run}:scan:{self.spectrum_id}:{self.peptidoform}"
        if as_url:
            usi = "http://proteomecentral.proteomexchange.org/usi/?usi=" + usi
        return usi
