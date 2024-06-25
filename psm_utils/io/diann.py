"""
Reader for PSM files from DIA-NN

Reads the '.tsv' file as defined on the `DIA-NN documentation page <https://github.com/vdemichev/DiaNN>`_.
"""

from __future__ import annotations

import csv
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Iterable, Optional
import re

import pyarrow.parquet as pq
from pyteomics import mass

from psm_utils.io._base_classes import ReaderBase
from psm_utils.io._utils import set_csv_field_size_limit
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

set_csv_field_size_limit()

class DIANNReader(ReaderBase, ABC):
    def __init__(
    self, filename, score_column: str = "CScore", *args, **kwargs
    ) -> None:
        """
        Reader for DIA-NN '.tsv' file.

        Parameters
        ----------
        filename : str or Path
            Path to PSM file.
        score_column: str, optional
            Name of the column that holds the primary PSM score. Default is
            ``CScore``.

        """
        super().__init__(filename, *args, **kwargs)
        self.filename = filename
        self.score_column = score_column

    def __iter__(self) -> Iterable[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        with open(self.filename) as msms_in:
            reader = csv.DictReader(msms_in, delimiter="\t")
            for row in reader:
                yield self._get_peptide_spectrum_match(row)

    def _get_peptide_spectrum_match(self, psm_dict) -> PSM:
        """Parse a single PSM from a DIA-NN PSM file."""
        rescoring_features = {}
        for ft in RESCORING_FEATURES:
            try:
                rescoring_features[ft] = psm_dict[ft]
            except KeyError:
                continue

        return PSM(
            peptidoform=self._parse_peptidoform(
                psm_dict["Modified.Sequence"],
                psm_dict["Precursor.Charge"]),
            spectrum_id='NA', # DIA-NN does not output spectrum ID
            run=psm_dict["Run"],
            is_decoy=False,
            qvalue=psm_dict["Q.Value"],
            pep=float(psm_dict["PEP"]),
            score=float(psm_dict[self.score_column]),
            retention_time=float(psm_dict["RT"]),
            ion_mobility=float(psm_dict["IM"]),
            protein_list=psm_dict["Protein.Names"].split(";"),
            source="diann",
            rank=None, # Leave out?
            provenance_data=({"diann_filename": str(self.filename)}),
            rescoring_features=rescoring_features,
            metadata={},
        )

    @staticmethod
    def _parse_peptidoform(peptide: str, charge: Optional[str]) -> str:
        if charge:
            peptide += f"/{int(float(charge))}"
        pattern = r"\(UniMod:(\d+)\)"
        replacement = r"[UNIMOD:\1]"
        peptide = re.sub(pattern, replacement, peptide)
        # If [UNIMOD:n] occurs before the first amino acid, a hyphen is added before the first amino acid
        if peptide[0] == "[":
            # Hyphen after the closing bracket
            peptide = peptide.replace("]", "]-", 1)
        return peptide

    def _parse_precursor_mz():
        return NotImplementedError("Method not implemented yet. DIA-NN does not yet output precursor m/z.")

    def from_dataframe(cls, dataframe) -> PSMList:
        """Create a PSMList from a DIA-NN Pandas DataFrame."""
        return PSMList(
            ptm_list=[
                cls._get_peptide_spectrum_match(cls(""), entry)
                for entry in dataframe.to_dict(orient="records")
            ]
        )


# TODO: Check
RESCORING_FEATURES = [
    "CScore",
    "RT",
    "Predicted.RT",
    "iRT",
    "Predicted.iRT",
    "Ms1.Profile.Corr",
    "Ms1.Area",
    "IM",
    "iIM"
]
