"""
Reader for PSM files from the Sage search engine.

Reads the ``results.sage.tsv`` file as defined on the
`Sage documentation page <https://github.com/lazear/sage/blob/v0.14.7/DOCS.md#interpreting-sage-output>`_.

"""

from __future__ import annotations

import csv
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Iterable, Optional

import pyarrow.parquet as pq
from pyteomics import mass

from psm_utils.io._base_classes import ReaderBase
from psm_utils.io._utils import set_csv_field_size_limit
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

set_csv_field_size_limit()


class _SageReaderBase(ReaderBase, ABC):
    def __init__(
        self, filename, score_column: str = "sage_discriminant_score", *args, **kwargs
    ) -> None:
        """
        Reader for Sage ``results.sage.tsv`` file.

        Parameters
        ----------
        filename : str or Path
            Path to PSM file.
        score_column: str, optional
            Name of the column that holds the primary PSM score. Default is
            ``sage_discriminant_score``, ``hyperscore`` could also be used.

        """
        super().__init__(filename, *args, **kwargs)
        self.filename = filename
        self.score_column = score_column

    @abstractmethod
    def __iter__(self) -> Iterable[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        raise NotImplementedError("Use `SageTSVReader` or `SageParquetReader` instead.")

    def _get_peptide_spectrum_match(self, psm_dict) -> PSM:
        """Parse a single PSM from a sage PSM file."""
        rescoring_features = {}
        for ft in RESCORING_FEATURES:
            try:
                rescoring_features[ft] = psm_dict[ft]
            except KeyError:
                continue

        return PSM(
            peptidoform=self._parse_peptidoform(
                psm_dict["peptide"],
                psm_dict["charge"],
            ),
            spectrum_id=psm_dict["scannr"],
            run=Path(psm_dict["filename"]).stem,
            is_decoy=psm_dict["is_decoy"],
            qvalue=psm_dict["spectrum_q"],
            score=float(psm_dict[self.score_column]),
            precursor_mz=self._parse_precursor_mz(psm_dict["expmass"], psm_dict["charge"]),
            retention_time=float(psm_dict["rt"]),
            protein_list=psm_dict["proteins"].split(";"),
            source="sage",
            rank=int(float(psm_dict["rank"])),
            provenance_data=({"sage_filename": str(self.filename)}),
            rescoring_features=rescoring_features,
            metadata={},
        )

    @staticmethod
    def _parse_peptidoform(peptide: str, charge: Optional[str]) -> str:
        if charge:
            peptide += f"/{int(float(charge))}"
        return peptide

    @staticmethod
    def _parse_precursor_mz(expmass: str, charge: Optional[str]) -> Optional[float]:
        if charge:
            charge = float(charge)
            expmass = float(expmass)
            return (expmass + (mass.nist_mass["H"][1][0] * charge)) / charge
        else:
            return None

    @classmethod
    def from_dataframe(cls, dataframe) -> PSMList:
        """Create a PSMList from a Sage Pandas DataFrame."""
        return PSMList(
            psm_list=[
                cls._get_peptide_spectrum_match(cls(""), entry)
                for entry in dataframe.to_dict(orient="records")
            ]
        )


class SageTSVReader(_SageReaderBase):
    def __iter__(self) -> Iterable[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        with open(self.filename, "r") as open_file:
            reader = csv.DictReader(open_file, delimiter="\t")
            for row in reader:
                row["is_decoy"] = (
                    True if row["label"] == "-1" else False if row["label"] == "1" else None
                )

                yield self._get_peptide_spectrum_match(row)

SageReader = SageTSVReader  # Alias for backwards compatibility


class SageParquetReader(_SageReaderBase):
    def __iter__(self) -> Iterable[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        with pq.ParquetFile(self.filename) as pq_file:
            for batch in pq_file.iter_batches():
                for row in batch.to_pylist():
                    yield self._get_peptide_spectrum_match(row)


RESCORING_FEATURES = [
    "expmass",
    "calcmass",
    "delta_mass",
    "peptide_len",
    "missed_cleavages",
    "isotope_error",
    "precursor_ppm",
    "fragment_ppm",
    "hyperscore",
    "delta_next",
    "delta_best",
    "delta_rt_model",
    "aligned_rt",
    "predicted_rt",
    "matched_peaks",
    "longest_b",
    "longest_y",
    "longest_y_pct",
    "matched_intensity_pct",
    "scored_candidates",
    "poisson",
    # "ms1_intensity",  # Removed in Sage v0.14
    "ms2_intensity",
]
