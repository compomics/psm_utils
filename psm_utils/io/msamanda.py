"""Interface to MS Amanda CSV result files."""

from __future__ import annotations

import logging
import csv
import re
from pathlib import Path
from itertools import compress
from typing import Union

import numpy as np

from psm_utils.exceptions import PSMUtilsException
from psm_utils.io._base_classes import ReaderBase
from psm_utils.psm import PSM, Peptidoform
from psm_utils.psm_list import PSMList

logger = logging.getLogger(__name__)


REQUIRED_COLUMNS = [
    "Title",
    "Sequence",
    "Modifications",
    "Protein Accessions",
    "Amanda Score",
    "m/z",
    "Charge",
    "RT",
    "Filename",
]


class MSAmandaReader(ReaderBase):
    """Reader for psm_utils TSV format."""

    def __init__(self, filename: Union[str, Path], *args, **kwargs) -> None:
        super().__init__(filename, *args, **kwargs)

    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""

        with open(self.filename, "rt") as open_file:
            if not next(open_file).startswith("#"):
                open_file.seek(0)
            reader = csv.DictReader(open_file, delimiter="\t")
            self._evaluate_columns(reader.fieldnames)
            for psm_dict in reader:
                yield self._get_peptide_spectrum_match(psm_dict)

    def read_file(self) -> PSMList:
        """Read full PSM file into a PSMList object."""
        return PSMList(psm_list=[psm for psm in self.__iter__()])

    @staticmethod
    def _evaluate_columns(columns) -> bool:
        """Case insensitive column evaluation MsAmanda file."""

        if "Rank" in columns:
            REQUIRED_COLUMNS.append("Rank")

        column_check = [True if col in columns else False for col in REQUIRED_COLUMNS]
        if not all(column_check):
            raise MSAmandaParsingError(
                f"Missing columns: {list(compress(REQUIRED_COLUMNS, list(~np.array(column_check))))}"
            )

    def _get_peptide_spectrum_match(
        self, psm_dict: dict[str, Union[str, float]]
    ) -> PSM:
        """Return a PSM object from MaxQuant msms.txt PSM file."""

        psm = PSM(
            peptidoform=self._parse_peptidoform(
                psm_dict["Sequence"], psm_dict["Modifications"], psm_dict["Charge"]
            ),
            spectrum_id=psm_dict["Title"],
            run=psm_dict["Filename"]
            .replace(".mgf", "")
            .replace(".raw", "")
            .replace(".mzml", ""),
            is_decoy=psm_dict["Protein Accessions"].startswith("REV_"),
            score=float(psm_dict["Amanda Score"]),
            precursor_mz=float(psm_dict["m/z"]),
            retention_time=float(psm_dict["RT"]),
            protein_list=psm_dict["Protein Accessions"].split(";"),
            source="MSAmanda",
            provenance_data=({"MSAmanda_filename": str(self.filename)}),
            metadata={
                col: str(psm_dict[col])
                for col in psm_dict.keys()
                if col not in REQUIRED_COLUMNS
            },
        )

        try:
            psm["rank"] = psm_dict["Rank"]
        except KeyError:
            pass

        return psm

    @staticmethod
    def _parse_peptidoform(seq, modifications, charge):
        "Parse MSAmanda sequence, modifications and charge to proforma sequence"
        peptide = [""] + [aa.upper() for aa in seq] + [""]

        pattern = re.compile(
            r"(?P<site>[A-Z])(?P<loc>-term|\d+)\((?P<mod_name>[A-Za-z]+)\|([-0-9.]+)\|(variable|fixed)\);?"
        )

        for match in pattern.finditer(modifications):
            if match.group("loc") == "-term":
                if match.group("site") == "N":
                    peptide[0] = peptide[0] + f'[{match.group("mod_name")}]'
                elif match.group("site") == "C":
                    peptide[-1] = peptide[-1] + f'[{match.group("mod_name")}]'
            else:
                peptide[int(match.group("loc"))] = (
                    peptide[int(match.group("loc"))] + f'[{match.group("mod_name")}]'
                )

        peptide[0] = peptide[0] + "-" if peptide[0] else ""
        peptide[-1] = "-" + peptide[-1] if peptide[-1] else ""
        proforma_seq = "".join(peptide)

        return Peptidoform(proforma_seq + f"/{charge}")


class MSAmandaParsingError(PSMUtilsException):
    """Error while parsing MaxQuant msms.txt PSM file."""
