from __future__ import annotations

from typing import Any, Dict, List, Optional, Union

from pydantic import BaseModel

from psm_utils.peptidoform import Peptidoform


class PeptideSpectrumMatch(BaseModel):  # TODO: Rename `PeptideSpectrumMatch` to `PSM`?
    """
    Data class representing a peptide-spectrum match (PSM).

    Links a :class:`~psm_utils.peptidoform.Peptidoform` to an observed spectrum
    and holds the related information. Attribute types are coerced and enforced upon
    initialization.

    Parameters
    ----------
    peptide : Peptidoform, str
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
    precursor_mz : float, optional
        Precursor m/z.
    retention_time : float, optional
        Retention time.
    protein_list : list[str]
        List of proteins or protein groups associated with peptide.
    rank : int
        rank of a psm
    source : str, optional
        PSM file type where PSM was stored. E.g., ``MaxQuant``.
    provenance_data : dict[str, str], optional
        Freeform dict to hold data describing the PSM origin, e.g. a search
        engine-specific identifier.
    metadata : dict[str, str], optional
        More data about PSM.
    rescoring_features : dict[str, str], optional
        Dict with features that can be used for PSM rescoring.

    """

    peptide: Union[Peptidoform, str]
    spectrum_id: Union[int, str]
    run: Optional[str] = None
    collection: Optional[str] = None
    # TODO: Not yet sure if `spectrum` should be included...
    spectrum: Optional[Any] = None
    is_decoy: Optional[bool] = None
    score: Optional[float] = None
    qvalue: Optional[float] = None
    pep: Optional[float] = None
    precursor_mz: Optional[float] = None
    retention_time: Optional[float] = None
    protein_list: Optional[List[str]] = None
    rank: Optional[int] = None
    source: Optional[str] = None
    provenance_data: Optional[Dict[str, str]] = None
    metadata: Optional[Dict[str, str]] = None
    rescoring_features: Optional[Dict[str, str]] = None

    class Config:
        arbitrary_types_allowed = True  # Allow non-pydantic class Peptidoform

    def __init__(self, **data):
        super().__init__(**data)
        # Parse peptidoform
        if isinstance(self.peptide, str):
            self.peptide = Peptidoform(self.peptide)
        elif not isinstance(self.peptide, Peptidoform):
            raise TypeError(
                f"Peptidoform or str expected for `peptide`, not `{type(self.peptide)}`."
            )

    def __getitem__(self, item) -> any:
        # return self.__dict__[item]
        return getattr(self, item)

    def __setitem__(self, item, value: any) -> None:
        setattr(self, item, value)

    def get_precursor_charge(self) -> int:
        """Precursor charge, as embedded in :py:attr:`peptide`."""
        return self.peptide.precursor_charge

    def universal_spectrum_identifier(self, as_url=False) -> str:
        """
        Compile USI for PeptideSpectrumMatch.

        Parameters
        ----------
        as_url : bool, optional
            Return URL to proteomeXchange.org USI aggregator.
        """
        # TODO: Support nativeID spectrum flag
        usi = f"mzspec:{self.collection}:{self.run}:scan:{self.spectrum_id}:{self.peptide}"
        if as_url:
            usi = "http://proteomecentral.proteomexchange.org/usi/?usi=" + usi
        return usi
