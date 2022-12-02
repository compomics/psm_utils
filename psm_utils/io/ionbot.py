"""
ionbot first csv file parser
"""

from __future__ import annotations

import csv
import re
from collections import namedtuple
from pathlib import Path
from typing import Iterable, NamedTuple, Optional

import pandas as pd

from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

REQUIRED_COLUMNS = [
    "database_peptide",
    "modifications",
    "charge",
    "spectrum_title",
    "spectrum_file",
    "proteins",
    "observed_retention_time",
    "database",
    "psm_score",
    "q-value",
    "PEP",
]


class IonbotReader(ReaderBase):
    def __init__(
        self,
        filename: str | Path,
        *args,
        **kwargs,
    ) -> None:

        super().__init__(filename, *args, **kwargs)
        self.filename = filename

    def __iter__(self) -> Iterable[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        with open(self.filename, "rt") as open_file:
            reader = csv.DictReader(open_file, delimiter=",")
            for row in reader:
                psm = self._get_peptide_spectrum_match(row)
                yield psm

    def read_file(self) -> PSMList:
        """Read full Peptide Record PSM file into a PSMList object."""
        psm_list = []
        with open(self.filename) as ionbot_in:
            reader = csv.DictReader(ionbot_in, delimiter=",")
            for row in reader:
                psm_list.append(self._get_peptide_spectrum_match(row))
        return PSMList(psm_list=psm_list)

    def _get_peptide_spectrum_match(self, psm_dict: dict[str, str | float]) -> PSM:

        return PSM(
            peptidoform=self._parse_peptidoform(
                psm_dict["database_peptide"],
                psm_dict["modifications"],
                psm_dict["charge"],
            ),
            spectrum_id=psm_dict["spectrum_title"],
            run=psm_dict["spectrum_file"],
            is_decoy=psm_dict["database"] == "D",
            score=float(psm_dict["psm_score"]),
            # precursor_mz=float(psm_dict["m/z"]),
            retention_time=float(psm_dict["observed_retention_time"]),
            protein_list=psm_dict["proteins"].split(
                "|"
            ),  # what is the ionbot separator?
            source="Ionbot",
            qvalue=float(psm_dict["q-value"]),
            pep=float(psm_dict["PEP"]),
            provenance_data=({"Ionbot_filename": str(self.filename)}),
            metadata={
                col: str(psm_dict[col])
                for col in psm_dict.keys()
                if col not in REQUIRED_COLUMNS
            },
        )

    @staticmethod
    def _parse_peptidoform(peptide, modifications, charge):
        peptide = peptide = [""] + list(peptide) + [""]
        pattern = re.compile(r"^(?P<U>\[\S*?\])?(?P<mod>.*?)(?P<AA>\[\S*?\])?$")

        for position, label in zip(
            modifications.split("|")[::2], modifications.split("|")[1::2]
        ):
            mod_match = pattern.search(label)

            if mod_match.group("U"):
                parsed_label = "U:" + mod_match.group("U")[1:-1]

            else:
                parsed_label = mod_match.group("mod")

            # if (mod_match.group("AA")) and (
            #     mod_match.group("AA")[1:-1] != peptide[int(position)]
            # ):
            #     print(mod_match.group("AA")[1:-1], peptide[int(position)])

            peptide[int(position)] += f"[{parsed_label}]"

        peptide[0] = peptide[0] + "-" if peptide[0] else ""
        peptide[-1] = "-" + peptide[-1] if peptide[-1] else ""
        proforma_seq = "".join(peptide)

        # Add charge state
        if charge:
            proforma_seq += f"/{charge}"

        return proforma_seq


class InvalidPeprecError(PSMUtilsIOException):
    """Invalid Peptide Record file."""

    pass


class InvalidIonbotModificationError(InvalidPeprecError):
    """Invalid Peptide Record modification."""

    pass
