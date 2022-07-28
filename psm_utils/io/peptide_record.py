"""Interface with Peptide Record PSM files."""

from __future__ import annotations

import csv
from collections import namedtuple
from pathlib import Path
from typing import Optional, Union

import pandas as pd

from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.psm import PeptideSpectrumMatch
from psm_utils.psm_list import PSMList


class PeptideRecord:
    """Peptide Record (PEPREC) file."""

    required_columns = [
        "peptide",
        "modifications"
    ]
    optional_columns = ["spec_id","charge","observed_retention_time", "label", "score"]

    def __init__(
        self,
        filename: Union[str, Path],
        required_columns: list[str] = None,
        optional_columns: list[str] = None,
    ) -> None:
        """
        Peptide Record (PEPREC) file.

        Helper class for handling PeptideRecord class. Upon initialization, the header
        is validated, separator inferred, and presence of required columns checked.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to PSM file.
        required_columns: list[str]
            Override default columns.
        optional_columns: list[str]
            Override default columns.

        """
        self.filename = filename
        if required_columns:
            self.required_columns = required_columns
        if optional_columns:
            self.optional_columns = optional_columns

        self._validate_header(filename)
        self.separator = self._infer_separator(filename)
        self._validate_required_columns()

    @staticmethod
    def peprec_to_proforma(
        peptide: str, modifications: str, label_mapping: Optional[dict[str, str]] = None
    ) -> str:
        """
        Convert PeptideRecord (PEPREC) peptide and modifications to ProForma string.

        Parameters
        ----------
        peptide: str
            Stripped peptide sequence.
        modifications: str
            Peptide modifications in PEPREC format (e.g., ``4|Oxidation``)
        label_mapping: dict[str, str], optional
            Mapping to replace modifications labels from :py:obj:`key` to
            :py:obj:`value`.

        Returns
        -------
        proforma_seq: str
            ProForma sequence

        Raises
        ------
        InvalidPeprecModificationError
            If a PEPREC modification cannot be parsed.

        """
        if not label_mapping:
            label_mapping = {}

        # List of peptide sequence with added terminal positions
        peptide = [""] + list(peptide) + [""]

        # Add modification labels
        for position, old_label in zip(
            modifications.split("|")[::2], modifications.split("|")[1::2]
        ):
            # Apply label_mapping
            if old_label in label_mapping:
                new_label = label_mapping[old_label]
            else:
                new_label = old_label

            try:
                peptide[int(position)] += f"[{new_label}]"
            except ValueError:
                raise InvalidPeprecModificationError(
                    f"Could not parse PEPREC modification `{modifications}`."
                )

        # Add dashes between residues and termini
        peptide[0] = peptide[0] + "-" if peptide[0] else ""
        peptide[-1] = "-" + peptide[-1] if peptide[-1] else ""

        proforma_seq = "".join(peptide)
        return proforma_seq

    @staticmethod
    def _validate_header(filename: Union[str, Path]):
        """Validate PEPREC header."""
        with open(filename, "rt") as f:
            line = f.readline()

    @staticmethod
    def _infer_separator(filename: Union[str, Path]) -> str:
        """Infer separator in PEPREC file."""
        with open(filename, "rt") as f:
            line = f.readline()
            if "\t" in line:
                separator = "\t"
            else:
                separator = ","
        return separator

    def _validate_required_columns(self) -> None:
        """Raise InvalidPeprecError if not all required columns are present."""
        with open(self.filename, "rt") as f:
            reader = csv.reader(f, delimiter=self.separator)
            columns = next(reader)
            for rc in self.required_columns:
                if rc not in columns:
                    raise InvalidPeprecError(f"Required column missing: `{rc}`")


class PeptideRecordReader(ReaderBase):
    def __init__(
        self,
        filename: Union[str, Path],
        modification_definitions: list[dict[str, str]] = None,
    ) -> None:
        """
        Reader for Peptide Record (PEPREC) PSM files.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to PSM file.
        modification_definitions: list[dict[str, str]]
            List of modification definition dictionaries. See
            `Parsing of modification labels`_ for more information.
        """

        super().__init__(filename, modification_definitions)
        self._peprec = PeptideRecord(self.filename)

        # TODO
        self.label_mapping = {}

        columns = self._peprec.required_columns + self._peprec.optional_columns
        self.PeprecEntry = namedtuple(
            "PeprecEntry", columns, defaults=[None for _ in columns]
        )

    def __next__(self) -> PeptideSpectrumMatch:
        return super().__next__()

    def __iter__(self):
        return super().__iter__()

    def read_file(self) -> PSMList:
        """Read full MaxQuant msms.txt PSM file into a PSMList object."""

        psm_list = []

        with open(self.filename) as peprec_in:
            reader = csv.DictReader(peprec_in, delimiter=self._peprec.separator)
            for row in reader:
                entry = self.PeprecEntry(**row)

                # Parse sequence and modifications
                proforma = PeptideRecord.peprec_to_proforma(
                    entry.peptide, entry.modifications, label_mapping=self.label_mapping
                )

                # Parse decoy label
                if entry.label:
                    is_decoy_map = {"-1": True, "1": False}
                    try:
                        is_decoy = is_decoy_map[entry.label]
                    except ValueError:
                        InvalidPeprecError(
                            f"Could not parse value for `label` {entry.label}. Should be `1` or `-1`."
                        )
                else:
                    is_decoy = None

                psm_list.append(PeptideSpectrumMatch(
                    peptide=proforma,
                    spectrum_id=entry.spec_id,
                    precursor_charge=entry.charge,
                    is_decoy=is_decoy,
                    retention_time=entry.observed_retention_time,
                    score=entry.score,
                    source="PeptideRecord",
                    provenance_data={"peprec_filename": self.filename}
                ))

        return PSMList(psm_list)


class PeptideRecordWriter(WriterBase):
    # TODO Implement
    pass


class InvalidPeprecError(Exception):
    """Invalid PEPREC file."""

    pass


class PeptideRecordOld:
    """Peptide record (PEPREC)."""

    def __init__(
        self,
        path: Optional[str] = None,
        context: str = "default",
        extra_required_columns: Optional[list[str]] = None,
    ):
        """
        Peptide record (PEPREC).

        Parameters
        ----------
        path: str, None
            Path to PEPREC file. If not None, file is read on instance creation.
        context: str
            Context of PEPREC. Is used for determining the required columns. Can be any
            of the following: {default, ms2rescore}.
        extra_required_columns: list[str]
            Extra required columns to be validated.

        Attributes
        ----------
        df: pandas.DataFrame
            DataFrame containing peptide record content.

        Methods
        -------
        from_csv(path, **kwargs)
            Read PEPREC from CSV.
        to_csv(path, **kwargs)
            Save PEPREC to CSV, if Path is None, overwrite existing PEPREC file.

        """
        self.path = path
        self.context = context
        self.extra_required_columns = extra_required_columns
        self.df = None

        if self.path:
            self.read_csv(self.path)

    def __repr__(self):
        """Get string representation of PeptideRecord."""
        return self.df.__repr__()

    @property
    def df(self) -> Optional[pd.DataFrame]:
        """Get DataFrame with PeptideRecord."""
        return self._df

    @df.setter
    def df(self, value: Optional[pd.DataFrame]):
        """Set DataFrame with PeptideRecord."""
        if isinstance(value, pd.DataFrame):
            self._validate_column_names(value)
        self._df = value

    def read_csv(self, path: str, **kwargs):
        """Read PEPREC from CSV."""
        self._validate_header(path)
        sep = self._infer_separator(path)
        self.df = pd.read_csv(path, sep=sep, index_col=None, **kwargs)
        self.df["modifications"] = self.df["modifications"].fillna("-")

    def to_csv(self, path: Optional[str] = None, **kwargs):
        """Save PEPREC to CSV, if Path is None, overwrite existing PEPREC file."""
        if not path and self.path:
            path = self.path
        elif not path and not self.path:
            raise ValueError("No known path to write PEPREC file.")
        self.df.to_csv(path, sep=" ", index=False, float_format="%.3f", **kwargs)

    def reorder_columns(self):
        """Reorder columns with canonical PeptideRecord columns first."""
        original_columns = self.df.columns.to_list()
        for col in self._required_columns:
            if col in original_columns:
                original_columns.remove(col)
            else:
                self._required_columns.remove(col)
        reordered_columns = self._required_columns + original_columns
        self.df = self.df[reordered_columns]

    @property
    def _required_columns(self):
        """Get required columns depending on PeptideRecord context."""
        required_columns = [
            "spec_id",
            "peptide",
            "modifications",
            "charge",
        ]
        if self.context == "ms2rescore":
            required_columns.extend(
                [
                    "psm_score",
                    "observed_retention_time",
                ]
            )
        elif self.context == "retention_time":
            required_columns.remove("modifications")
        return required_columns

    def _validate_header(self, path: str):
        """Validate PEPREC header."""
        with open(path, "rt") as f:
            line = f.readline()
            if line[:7] != "spec_id":
                raise InvalidPeprecError("PEPREC header should start with `spec_id`.")

    def _validate_column_names(self, df: pd.DataFrame):
        """Validate header of PEPREC file."""
        for col in self._required_columns:
            if col not in df.columns:
                raise InvalidPeprecError(
                    f"Required column `{col}` missing from header."
                )

    @classmethod
    def from_dataframe(cls, dataframe: pd.DataFrame):
        """Create PeptideRecord instance with data from a pandas DataFrame."""
        peprec = cls()
        peprec.df = dataframe
        peprec.reorder_columns()
        return peprec

    @classmethod
    def from_csv(cls, *args, **kwargs):
        """Create PeptideRecord instance with data from CSV file."""
        peprec = cls()
        peprec.read_csv(*args, **kwargs)
        return peprec

    @staticmethod
    def _infer_separator(path: str) -> str:
        """Infer separator in PEPREC file."""
        with open(path, "rt") as f:
            line = f.readline()
            separator = line[7]
        return separator

    def validate_decoy_presence(self):
        "check whether decoys are present in peprec file"

        if sum(self.df["Label"] == -1) == 0:
            raise InvalidPeprecError(
                "No decoys present, make sure to provide decoys along targets"
            )
        else:
            pass


class InvalidPeprecError(PSMUtilsIOException):
    """Invalid Peptide Record file."""

    pass


class InvalidPeprecModificationError(InvalidPeprecError):
    """Invalid Peptide Record modification."""

    pass
