from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Iterable, Optional

from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList


class SageReader(ReaderBase):
    def __init__(self, filename, *args, **kwargs) -> None:
        super().__init__(filename, *args, **kwargs)
        self.filename = filename

    def __iter__(self) -> PSMList:
        """Read full Peptide Record PSM file into a PSMList object."""
        psm_list = []
        with open(self.filename) as open_file:
            reader = csv.DictReader(open_file, delimiter="\t")
            for row in reader:
                psm = self._get_peptide_spectrum_match(row)
                yield psm

    def read_file(self) -> PSMList:
        """Read full PSM file into a PSMList object."""
        psm_list = []
        for psm in self.__iter__():
            psm_list.append(psm)
        return PSMList(psm_list=psm_list)

    def _get_peptide_spectrum_match(self, psm_dict) -> PSM:
        """Parse a single PSM from a sage PSM file."""

        psm_dict["delta_mass"] = float(psm_dict["expmass"]) - float(
            psm_dict["calcmass"]
        )
        return PSM(
            peptidoform=self._parse_peptidoform(
                psm_dict["peptide"],
                psm_dict["charge"],
            ),
            spectrum_id=psm_dict["scannr"],
            run=psm_dict["filename"].split(".")[0],
            is_decoy=psm_dict["label"] == "-1",
            qvalue=psm_dict["spectrum_fdr"],
            score=float(psm_dict["hyperscore"]),
            precursor_mz=None,
            retention_time=float(psm_dict["rt"]),
            protein_list=psm_dict["proteins"].split(";"),
            source="sage",
            rank=int(float(psm_dict["rank"])),
            provenance_data=({"sage_filename": str(self.filename)}),
            rescoring_features={
                ft: psm_dict[ft]
                for ft in [
                    "expmass",
                    "calcmass",
                    "delta_mass",
                    "peptide_len",
                    "missed_cleavages",
                    "isotope_error",
                    "precursor_ppm",
                    "fragment_ppm",
                    "delta_hyperscore",
                    "aligned_rt",
                    "matched_peaks",
                    "longest_b",
                    "longest_y",
                    "longest_y_pct",
                    "matched_intensity_pct",
                    "scored_candidates",
                    "poisson",
                    "sage_discriminant_score",
                    "ms1_intensity",
                ]
            },
            metadata={},
        )

    @staticmethod
    def _parse_peptidoform(peptide, charge):
        # Add charge state
        if charge:
            peptide += f"/{int(float(charge))}"

        return peptide
