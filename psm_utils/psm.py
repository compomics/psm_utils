from __future__ import annotations

from dataclasses import field
from typing import Any, Dict, List, Optional, Union

from pydantic.dataclasses import dataclass
from pyteomics import proforma

from psm_utils.peptidoform import Peptidoform


class _PydanticConfig:
    arbitrary_types_allowed = True


@dataclass(config=_PydanticConfig)
class PeptideSpectrumMatch:  # TODO: Rename `PeptideSpectrumMatch` to `PSM`?
    """
    Data class representing a peptide-spectrum match (PSM).

    Links a :class:`psm_utils.peptidoform.Peptidoform` to an observed spectrum
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
    precursor_mz: Optional[float] = field(default=None)
    retention_time: Optional[float] = field(default=None)
    protein_list: Optional[List[str]] = field(default=None)
    source: Optional[str] = field(default=None)
    provenance_data: Optional[Dict[str, str]] = field(default=None)
    metadata: Optional[Dict[str, str]] = field(default=None)
    rescoring_features: Optional[Dict[str, str]] = field(default=None)

    def __post_init__(self):
        # Parse peptidoform
        if isinstance(self.peptide, str):
            self.peptide = Peptidoform(self.peptide)
        elif not isinstance(self.peptide, Peptidoform):
            raise TypeError(
                f"Peptidoform or str expected for `peptide`, not `{type(self.peptide)}`."
            )

    @property
    def precursor_charge(self) -> int:
        """Precursor charge, as embedded in :py:attr:`peptide`."""
        try:
            return self.peptide.properties["charge_state"].charge
        except (AttributeError, KeyError):
            return None

    @precursor_charge.setter
    def precursor_charge(self, charge: Union[str, int]) -> int:
        if isinstance(charge, str):
            charge = charge.strip("+")
        self.peptide.properties["charge_state"] = proforma.ChargeState(charge)

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
