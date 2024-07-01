"""
Interface with ionbot PSM files.

Currently only supports the ionbot.first.csv files.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Dict, Iterable, Union

from psm_utils.io._base_classes import ReaderBase
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList
from psm_utils.io._utils import set_csv_field_size_limit

set_csv_field_size_limit()

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
        """
        Reader for ``ionbot.first.csv`` PSM files.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to PSM file.

        Examples
        --------

        IonbotReader supports iteration:

        >>> from psm_utils.io.ionbot import IonbotReader
        >>> for psm in IonbotReader("ionbot.first.csv"):
        ...     print(psm.peptidoform.proforma)
        ACDEK
        AC[Carbamidomethyl]DEFGR
        [Acetyl]-AC[Carbamidomethyl]DEFGHIK

        Or a full file can be read at once into a :py:class:`psm_utils.psm_list.PSMList`
        object:

        >>> ionbot_reader = IonbotReader("ionbot.first.csv")
        >>> psm_list = ionbot_reader.read_file()

        """
        super().__init__(filename, *args, **kwargs)
        self.filename = filename

    def __iter__(self) -> Iterable[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        with open(self.filename, "rt") as open_file:
            reader = csv.DictReader(open_file, delimiter=",")
            for row in reader:
                yield self._get_peptide_spectrum_match(row)

    def read_file(self) -> PSMList:
        """Read full PSM file into a PSMList object."""
        return PSMList(psm_list=[psm for psm in self])

    def _get_peptide_spectrum_match(self, psm_dict: Dict[str, str | float]) -> PSM:
        return PSM(
            peptidoform=self._parse_peptidoform(
                psm_dict["matched_peptide"],
                psm_dict["modifications"],
                psm_dict["charge"],
            ),
            spectrum_id=psm_dict["spectrum_title"],
            run=psm_dict["spectrum_file"],
            is_decoy=(
                True
                if psm_dict["database"] == "D"
                else False if psm_dict["database"] == "T" else None
            ),
            score=float(psm_dict["psm_score"]),
            precursor_mz=float(psm_dict["m/z"]),
            retention_time=float(psm_dict["observed_retention_time"]),
            protein_list=psm_dict["proteins"].split("||"),
            source="ionbot",
            qvalue=float(psm_dict["q-value"]),
            pep=float(psm_dict["PEP"]),
            provenance_data=({"ionbot_filename": str(self.filename)}),
            metadata={
                col: str(psm_dict[col]) for col in psm_dict.keys() if col not in REQUIRED_COLUMNS
            },
        )

    @staticmethod
    def _parse_peptidoform(
        peptide: str, modifications: str, charge: Union[str, int]
    ) -> Peptidoform:
        """Parse peptide, modifications, and charge to Peptidoform."""
        # Split peptide into list of amino acids with termini
        peptide = peptide = [""] + list(peptide) + [""]

        # Add modifications
        pattern = re.compile(r"^(?P<U>\[\S*?\])?(?P<mod>.*?)(?P<AA>\[\S*?\])?$")
        for position, label in zip(modifications.split("|")[::2], modifications.split("|")[1::2]):
            mod_match = pattern.search(label)
            if mod_match.group("U"):
                parsed_label = "U:" + mod_match.group("U")[1:-1]
            else:
                parsed_label = mod_match.group("mod")
            peptide[int(position)] += f"[{parsed_label}]"

        # Add terminal modifications
        peptide[0] = peptide[0] + "-" if peptide[0] else ""
        peptide[-1] = "-" + peptide[-1] if peptide[-1] else ""
        proforma_seq = "".join(peptide)

        # Add charge state
        proforma_seq += f"/{charge}"

        return Peptidoform(proforma_seq)


class InvalidIonbotModificationError(PSMUtilsIOException):
    """Invalid Peptide Record modification."""

    pass
