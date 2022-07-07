from dataclasses import dataclass
from typing import List

import pandas as pd

from psm_utils.psm import PeptideSpectrumMatch
from psm_utils.io.peptide_record import PeptideRecord


@dataclass
class PSMList:
    """
    Data class representing a list of PSMs.

    Parameters
    ----------
    psm_list : list[PeptideSpectrumMatch]
        List of PeptideSpectrumMatch instances.

    """
    psm_list : List[PeptideSpectrumMatch]

    @classmethod
    def from_csv(cls) -> 'PSMList':
        """Read PSMList from comma-separated values file."""
        raise NotImplementedError

    def to_csv(self) -> None:  # Or "to_dataframe"?
        """Write PSMList to comma-separated values file."""
        raise NotImplementedError

    def to_peprec(self) -> PeptideRecord:
        """Convert PSMList to MS²PIP PEPREC format."""
        return PeptideRecord.from_dataframe(
            pd.DataFrame(
                [psm.to_peprec_entry() for psm in self.psm_list]
            )
        )

    def to_pin(self):
        """Write PSMList to Percolator In."""
        raise NotImplementedError


