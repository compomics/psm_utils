"""Interface to MaxQuant msms.txt PSM files."""

from __future__ import annotations

import csv
import logging
import re
from itertools import compress
from pathlib import Path
from typing import Union

import numpy as np

from psm_utils.exceptions import PSMUtilsException
from psm_utils.io._base_classes import ReaderBase
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PeptideSpectrumMatch
from psm_utils.psm_list import PSMList

logger = logging.getLogger(__name__)

MSMS_DEFAULT_COLUMNS = {
    "Raw file",
    "Scan number",
    "Charge",
    "Length",
    "Sequence",
    "Modified sequence",
    "Proteins",
    "Missed cleavages",
    "Mass",
    "Mass error [Da]",
    "Mass error [ppm]",
    "Reverse",
    "Retention time",
    "PEP",
    "Score",
    "Delta score",
    "Localization prob",
    "Matches",
    "Intensities",
    "Mass Deviations [Da]",
    "Mass Deviations [ppm]",
    "Intensity coverage",
    "id",
}


class MSMSReader(ReaderBase):
    """Reader for MaxQuant msms.txt PSM files."""

    def __init__(
        self,
        filename: Union[str, Path],
    ) -> None:
        """
        Reader for MaxQuant msms.txt PSM files.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to PSM file.
        decoy_prefix: str, optional
            Protein name prefix used to denote decoy protein entries. Default:
            ``"DECOY_"``.

        Examples
        --------
        :py:class:`MSMSReader` supports iteration:

        >>> from psm_utils.io.maxquant import MSMSReader
        >>> for psm in MSMSReader("msms.txt"):
        ...     print(psm.peptide.proforma)
        WFEELSK
        NDVPLVGGK
        GANLGEMTNAGIPVPPGFC[+57.022]VTAEAYK
        ...

        Or a full file can be read at once into a :py:class:`psm_utils.psm_list.PSMList`
        object:

        >>> reader = MSMSReader("msms.txt")
        >>> psm_list = reader.read_file()

        """

        super().__init__(filename)

        self._rename_mapping = None
        self._mass_error_unit = None

        self._validate_msms()

    def __iter__(self):
        """Iterate over file and return PSMs one-by-one"""
        with open(self.filename) as msms_in:
            reader = csv.DictReader(msms_in, delimiter="\t")
            for psm_dict in reader:
                psm = self._get_peptide_spectrum_match(psm_dict)
                yield psm

    def __next__(self) -> PeptideSpectrumMatch:
        return super().__next__()

    def read_file(self) -> PSMList:
        """Read full MaxQuant msms.txt PSM file into a PSMList object."""
        psm_list = []
        with open(self.filename) as msms_in:
            reader = csv.DictReader(msms_in, delimiter="\t")
            for psm_dict in reader:
                psm_list.append(self._get_peptide_spectrum_match(psm_dict))
        return PSMList(psm_list)

    def _validate_msms(self) -> None:
        with open(self.filename, "r") as msms_file:
            msms_reader = csv.DictReader(msms_file, delimiter="\t")
            self._evaluate_columns(msms_reader.fieldnames)
            self._set_mass_error_unit(msms_reader.fieldnames)
            self._rename_mapping = self._fix_column_case(msms_reader.fieldnames)

    @staticmethod
    def _evaluate_columns(columns) -> bool:
        """Case insensitive column evaluation msms file."""
        columns = list(map(lambda col: col.lower(), columns))
        column_check = [
            True if col.lower() in columns else False for col in MSMS_DEFAULT_COLUMNS
        ]
        if all(column_check):
            return True
        else:
            raise MSMSParsingError(
                f"Missing columns: {list(compress(columns, column_check))}"
            )

    @staticmethod
    def _fix_column_case(columns: list[str]) -> dict[str, str]:
        """
        Create mapping for column names with the correct case.

        Using `_evaluate_columns`, we can load required columns in a case-insensitive
        manner. As a result, the column name case must be fixed for downstream usage.
        """
        case_mapping = {col.lower(): col for col in MSMS_DEFAULT_COLUMNS}
        required_col = list(map(lambda col: col.lower(), MSMS_DEFAULT_COLUMNS))
        rename_mapping = {
            case_mapping[col.lower()]: col
            for col in columns
            if col.lower() in required_col
        }
        return rename_mapping

    def _set_mass_error_unit(self, columns) -> None:
        """Get mass error unit from DataFrame columns."""
        if "Mass error [Da]" in columns:
            self._mass_error_unit = "Da"
        elif "Mass error [ppm]" in columns:
            self._mass_error_unit = "ppm"
        else:
            raise NotImplementedError(f"MSMS.txt mass error unit not supported.")

    def _get_peptide_spectrum_match(
        self, psm_dict: dict[str, Union[str, float]]
    ) -> PeptideSpectrumMatch:
        """Return a PeptideSpectrumMatch object from MaxQuant msms.txt PSM file."""

        psm = PeptideSpectrumMatch(
            peptide=self._parse_peptidoform(
                psm_dict["Modified sequence"], psm_dict["Charge"]
            ),
            spectrum_id=psm_dict["Scan number"],
            run=psm_dict["Raw file"],
            is_decoy=psm_dict["Reverse"] == "+",
            score=float(psm_dict["Score"]),
            precursor_mz=float(psm_dict["m/z"]),
            retention_time=float(psm_dict["Retention time"]),
            protein_list=psm_dict["Proteins"].split(";"),
            source="msms",
            provenance_data=({"msms_filename": str(self.filename)}),
            metadata={
                "Delta Score": psm_dict["Delta score"],
                "Localization prob": psm_dict["Localization prob"],
                "PepLen": psm_dict["Length"],
                "Precursor Intensity": psm_dict["Precursor Intensity"],
                "dM": psm_dict[f"Mass error [{self._mass_error_unit}]"],
                "Matches": psm_dict["Matches"],
                "Intensities": psm_dict["Intensities"],
                "Intensity coverage": psm_dict["Intensity coverage"],
                f"Mass Deviations [{self._mass_error_unit}]": psm_dict[
                    self._rename_mapping[
                        f"Mass Deviations [{self._mass_error_unit}]"
                    ]
                ],
            },
        )
        return psm

    @staticmethod
    def _parse_peptidoform(modified_seq: str, charge: int) -> Peptidoform:
        """Parse modified sequence to :py:class:`psm_utils.peptidoform.Peptidoform`."""

        # pattern to match open and closed round brackets
        pattern = re.compile(r"\(((?:[^)(]+|\((?:[^)(]+|\([^)(]*\))*\))*)\)")
        modified_seq = modified_seq.strip("_")
        modified_seq_len = len(modified_seq)

        for match in pattern.finditer(modified_seq):
            # take match object string, remove leading and trailing bracket
            # and escape remaining brackets
            se_mod_string = match[1].replace("(", "\(").replace(")", "\)")

            # if N-term mod
            if match.start() == 0:
                modified_seq = re.sub(
                    f"\({se_mod_string}\)", f"[{match[1]}]-", modified_seq
                )

            # if C-term mod
            elif match.end() == modified_seq_len:
                modified_seq = re.sub(
                    f"\({se_mod_string}\)", f"-[{match[1]}]", modified_seq
                )

            # if modification on amino acid
            else:
                modified_seq = re.sub(
                    f"\({se_mod_string}\)", f"[{match[1]}]", modified_seq
                )

        modified_seq += f"/{charge}"

        return Peptidoform(modified_seq)


class MSMSParsingError(PSMUtilsException):
    """Error while parsing MaxQuant msms.txt PSM file."""
