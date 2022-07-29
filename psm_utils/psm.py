from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional, Union

from psm_utils.peptidoform import Peptidoform


# TODO: Rename `PeptideSpectrumMatch` to `PSM`?
@dataclass
class PeptideSpectrumMatch:
    """
    Data class representing a peptide-spectrum match (PSM).

    Links a :class:`psm_utils.peptidoform.Peptidoform` to an observed spectrum
    and holds relevant the metadata.

    Parameters
    ----------
    peptide : Peptidoform, str
        Peptidoform object or string in ProForma v2 notation.
    spectrum_id : str
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
        ``spectrum_utils.spectrum.MsmsSpectrum`` object.
    is_decoy : bool, optional
        Boolean specifying if the PSM is a decoy (``True``) or target hit
        (``False``).
    score : float, optional
        Search engine score.
    precursor_charge : int, optional
        Precursor charge.
    precursor_mz : float, optional
        Precursor m/z.
    retention_time : float, optional
        Retention time.
    source: str, optional
        PSM file type where PSM was stored. E.g., ``MaxQuant``.
    provenance_data: dict, optional
        Freefrom dict to hold data describing the PSM origin, e.g. a search
        engine-specific identifier.
    metadata : dict, optional
        More data about PSM.

    """

    peptide: Union[Peptidoform, str]
    spectrum_id: str
    run: Optional[str] = None
    collection: Optional[str] = None
    # TODO: Not yet sure if `spectrum` should be included...
    spectrum: Optional[Any] = None
    is_decoy: Optional[bool] = None
    score: Optional[float] = None
    precursor_charge: Optional[int] = field(default=None, repr=False)
    precursor_mz: Optional[float] = field(default=None, repr=False)
    retention_time: Optional[float] = field(default=None, repr=False)
    protein_list: Optional[str] = field(default=None, repr=False)
    source: Optional[str] = field(default=None, repr=False)
    provenance_data: Optional[dict] = field(default=None, repr=False)
    metadata: Optional[dict] = field(default=None, repr=False)

    def __post_init__(self):
        # Parse peptidoform
        if isinstance(self.peptide, str):
            self.peptide = Peptidoform(self.peptide)
        elif not isinstance(self.peptide, Peptidoform):
            raise TypeError(
                f"Peptidoform or str expected for `peptide`, not `{type(self.peptide)}`."
            )
        # Parse charge state
        if isinstance(self.precursor_charge, str):
            self.precursor_charge = int(self.precursor_charge.strip("+"))
        elif self.peptide.properties["charge_state"]:
            self.precursor_charge = self.peptide.properties["charge_state"].charge

    @property
    def universal_spectrum_identifier(self) -> str:
        """"""
        raise NotImplementedError

    def to_peprec_entry(self) -> dict:
        entry = {
            "spec_id": self.spectrum_id,
            "peptide": self.stripped_sequence,
            "modifications": self.modifications,
            "charge": self.precursor_charge,
            "is_decoy": int(self.is_decoy),
        }
        if self.retention_time:
            entry["observed_retention_time"] = self.retention_time
        if self.metadata:
            entry.update(self.metadata)
        return entry
