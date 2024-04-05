"""Reader for ProteoScape Parquet files."""

import logging
import re
from pathlib import Path
from typing import Union
from collections import namedtuple

import numpy as np
import pandas as pd

from psm_utils import PSM
from psm_utils.io._base_classes import ReaderBase
from psm_utils.peptidoform import format_number_as_string

logger = logging.getLogger(__name__)

DECOY_PATTERN = re.compile(r"^Reverse_")


class ProteoScapeReader(ReaderBase):
    """Reader for ProteoScape Parquet files."""

    def __init__(
        self,
        filename: Union[str, Path],
        *args,
        **kwargs,
    ) -> None:
        """
        Reader for ProteoScape Parquet files.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to MSF file.

        """
        if isinstance(filename, pd.DataFrame):
            self.data = filename
        else:
            super().__init__(filename, *args, **kwargs)
            self.data = pd.read_parquet(self.filename)

        self._Row = namedtuple("Row", self.data.columns)

    def __len__(self):
        """Return number of PSMs in file."""
        return len(self.data)

    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""
        for entry in self.data.itertuples():
            yield _parse_entry(entry)

    def __getitem__(self, index):
        """Return PSM at index."""
        return _parse_entry(self._Row(*self.data.iloc[index]))

    @classmethod
    def from_dataframe(cls, dataframe: pd.DataFrame, *args, **kwargs):
        """Create a ProteoScapeReader from a DataFrame."""
        return cls(dataframe, *args, **kwargs)


def _parse_peptidoform(
    stripped_peptide: str, ptms: np.ndarray, ptm_locations: np.ndarray, precursor_charge: int
) -> str:
    """Parse peptide sequence and modifications to ProForma."""
    peptidoform = list(stripped_peptide)
    n_term = ""
    c_term = ""
    for ptm, ptm_location in zip(ptms, ptm_locations):
        ptm = format_number_as_string(ptm)
        if ptm_location == -1:
            n_term = f"[{ptm}]-"
        elif ptm_location == len(peptidoform):
            c_term = f"-[{ptm}]"
        else:
            peptidoform[ptm_location] = f"{peptidoform[ptm_location]}[{ptm}]"
    return f"{n_term}{''.join(peptidoform)}{c_term}/{precursor_charge}"


def _parse_entry(entry) -> PSM:
    """Parse a single entry from ProteoScape Parquet file to PSM object."""
    return PSM(
        peptidoform=_parse_peptidoform(
            entry.stripped_peptide, entry.ptms, entry.ptm_locations, entry.precursor_charge
        ),
        spectrum_id=entry.ms2_id,
        run=getattr(entry, "run", None),
        is_decoy=all(DECOY_PATTERN.match(p) for p in entry.locus_name),
        score=entry.x_corr_score,
        precursor_mz=entry.precursor_mz,
        retention_time=entry.rt,
        ion_mobility=entry.ook0,
        protein_list=list(entry.locus_name),
        rank=entry.rank,
        source="ProteoScape",
        provenance_data={
            "candidate_id": str(entry.candidate_id),
            "ms2_id": str(entry.ms2_id),
            "parent_id": str(entry.parent_id),
        },
        metadata={
            "leading_aa": str(entry.leading_aa),
            "trailing_aa": str(entry.trailing_aa),
            "corrected_ook0": str(entry.corrected_ook0),
        },
        rescoring_features={
            "tims_score": float(entry.tims_score),
            "x_corr_score": float(entry.x_corr_score),
            "delta_cn_score": float(entry.delta_cn_score),
            "ppm_error": float(entry.ppm_error),
            "number_matched_ions": float(entry.number_matched_ions),
            "number_expected_ions": float(entry.number_expected_ions),
            "ion_proportion": float(entry.ion_proportion),
            "spectrum_total_ion_intensity": float(entry.spectrum_total_ion_intensity),
        },
    )
