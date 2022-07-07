from dataclasses import dataclass, field
from typing import Optional, Union

from psm_utils._exceptions import PeptidoformError
from psm_utils.peptidoform import Peptidoform

@dataclass
class PeptideSpectrumMatch:
    """
    Data Class representing a PeptideSpectrumMatch (PSM).

    Parameters
    ----------
    peptide : Peptidoform, str
        Peptidoform object or string in ProForma v2 notation.
    spectrum_id : str
        Spectrum identifier as used in spectrum file (e.g., mzML or MGF),
        usually in HUPO-PSI nativeID format (MS:1000767), e.g.,
        ``controllerType=0 controllerNumber=0 scan=423``
    run : str
        Name of the MS run. Usually the spectrum file filename without
        extension.
    is_decoy : bool, optional
        Boolean specifying if the PSM is a decoy (``True``) or target hit
        (``False``).
    precursor_charge : int, optional
        PSM precursor charge.
    precursor_mz : float, optional
        PSM precursor m/z.
    retention_time : float, optional
        PSM retention time.
    metadata : dict, optional
        More data about PSM, e.g., search engine score.

    """
    peptide : Union[Peptidoform, str]
    spectrum_id : str
    run : str
    is_decoy : Optional[bool] = None
    precursor_charge : Optional[int] = field(default=None, repr=False)
    precursor_mz : Optional[float] = field(default=None, repr=False)
    retention_time : Optional[float] = field(default=None, repr=False)
    metadata : Optional[dict] = field(default=None, repr=False)

    def __post_init__(self):
        # Parse peptidoform
        if isinstance(self.peptide, str):
            self.peptide = Peptidoform(self.peptide)
        # Parse charge state
        if isinstance(self.precursor_charge, str):
            self.precursor_charge = int(self.precursor_charge.strip("+"))
        elif self.peptide.modifiers["charge_state"]:
            self.precursor_charge = self.peptide.modifiers["charge_state"].charge

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
