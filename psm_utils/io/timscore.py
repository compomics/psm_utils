"""Reader for Proteome Discoverer MSF PSM files."""

import logging
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd

from psm_utils import PSM
from psm_utils.io._base_classes import ReaderBase

logger = logging.getLogger(__name__)


class TIMScoreReader(ReaderBase):
    """Reader for TIMScore Parquet files."""

    def __init__(
        self,
        filename: Union[str, Path],
        *args,
        **kwargs,
    ) -> None:
        """
        Reader for TIMScore Parquet files.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to MSF file.

        """
        super().__init__(filename, *args, **kwargs)

        self.data = pd.read_parquet(self.filename)

    def __len__(self):
        """Return number of PSMs in file."""
        return len(self.data)

    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""
        for entry in self.data.itertuples():
            yield PSM(
                peptidoform=_parse_peptidoform(
                    entry.stripped_peptide, entry.ptms, entry.ptm_locations, entry.precursor_charge
                ),
                spectrum_id=entry.ms2_id,
                is_decoy=None,  # TODO: Parse from protein?
                score=entry.tims_score,  # TODO: Correct score?
                precursor_mz=entry.precursor_mz,
                retention_time=entry.rt,
                ion_mobility=entry.corrected_ook0,
                protein_list=list(entry.locus_name),
                rank=entry.rank,
                source="TIMScore",
                provenance_data={
                    "candidate_id": str(entry.candidate_id),
                    "ms2_id": str(entry.ms2_id),
                    "parent_id": str(entry.parent_id),
                },
                metadata={
                    "leading_aa": str(entry.leading_aa),
                    "trailing_aa": str(entry.trailing_aa),
                    "ook0": str(entry.ook0),
                },
                rescoring_features={
                    "x_corr_score": float(entry.x_corr_score),
                    "delta_cn_score": float(entry.delta_cn_score),
                    "ppm_error": float(entry.ppm_error),
                    "number_matched_ions": float(entry.number_matched_ions),
                    "number_expected_ions": float(entry.number_expected_ions),
                    "ion_proportion": float(entry.ion_proportion),
                    "spectrum_total_ion_intensity": float(entry.spectrum_total_ion_intensity),
                },
            )


def _parse_peptidoform(
    stripped_peptide: str, ptms: np.ndarray, ptm_locations: np.ndarray, precursor_charge: int
) -> str:
    """Parse peptide sequence and modifications to ProForma."""
    # TODO: How are terminal modifications handled?
    peptidoform = list(stripped_peptide)
    for ptm, ptm_location in zip(ptms, ptm_locations):
        peptidoform[ptm_location] = f"{peptidoform[ptm_location]}[{ptm}]"
    peptidoform.append(f"/{precursor_charge}")
    return "".join(peptidoform)
