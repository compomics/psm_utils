"""
Interface with Peptide Record PSM files.

Peptide Record (or PEPREC) is a legacy PSM file type developed at CompOmics as input
format for `MS²PIP <https://github.com/compomics/ms2pip_c/>`_. It is a simple and
flexible delimited text file where each row represents a single PSM. Required columns
are:

- ``peptide``: Simple, stripped peptide sequence (e.g., ``ACDE``).
- ``modifications``: Amino acid modifications in a custom format (see below).

Depending on the use case, more columns can be required or optional:

- ``spec_id``: Spectrum identifier; usually the identifier used in the spectrum file.
- ``charge``: Peptide precursor charge.
- ``observed_retention_time``: Observed retention time.
- ``predicted_retention_time``: Predicted retention time.
- ``label``: Target/decoy: ``1`` for target PSMs, ``-1`` for decoy PSMs.
- ``score``: Primary search engine score (e.g., the score used for q-value calculation).

Peptide modifications are denoted as a pipe-separated list of pipe-separated
*location → label* pairs for each modification. The location is an integer counted
starting at 1 for the first amino acid. ``0`` is reserved for N-terminal modifications
and ``-1`` for C-terminal modifications. Unmodified peptides can be marked with a hyphen
(``-``). For example:

==================================== ==============================================================================
 PEPREC modification(s)               Explanation
==================================== ==============================================================================
 ``-``                               Unmodified
 ``1|Oxidation``                     ``Oxidation`` on the first amino acid
 ``1|Oxidation|5|Carbamidomethyl``   ``Oxidation`` on the first amino acid and ``Carbamidomethyl`` on the fifth
 ``0|Acetylation``                   ``Acetylation`` on the N-terminus
==================================== ==============================================================================


Full PEPREC example:

.. code-block::

    spec_id,modifications,peptide,charge
    peptide1,-,ACDEK,2
    peptide2,2|Carbamidomethyl,ACDEFGR,3
    peptide3,0|Acetyl|2|Carbamidomethyl,ACDEFGHIK,2

"""

from __future__ import annotations

import csv
from collections import namedtuple
from pathlib import Path
from typing import Iterable, NamedTuple, Optional, Union

import pandas as pd

from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.psm import PeptideSpectrumMatch
from psm_utils.psm_list import PSMList


class PeptideRecord:

    required_columns = ["peptide", "modifications"]
    optional_columns = [
        "spec_id",
        "charge",
        "observed_retention_time",
        "predicted_retention_time",
        "label",
        "score",
    ]

    def __init__(
        self,
        filename: Union[str, Path],
        required_columns: list[str] = None,
        optional_columns: list[str] = None,
    ) -> None:
        """
        Helper class for handling PeptideRecord files.

        Upon initialization, the separator inferred and presence of required columns
        is checked.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to PSM file.
        required_columns: list[str]
            Override default columns.
        optional_columns: list[str]
            Override default columns.

        Attributes
        ----------
        separator: str
            Separator (delimiter) used in PeptideRecord file.

        Raises
        ------
        InvalidPeprecError
            If PeptideRecord separator cannot be inferred from header.

        """
        self.filename = filename
        self.separator = None

        if required_columns:
            self.required_columns = required_columns
        else:
            self.required_columns = self.required_columns.copy()  # Copy from class
        if optional_columns:
            self.optional_columns = optional_columns
        else:
            self.optional_columns = self.optional_columns.copy()  # Copy from class

        self._infer_separator()
        self._validate_required_columns()

    def __repr__(self) -> str:
        return f"PeptideRecord('{self.filename}')"

    def _infer_separator(self) -> None:
        """Infer separator used in PeptideRecord file."""
        with open(self.filename, "rt") as f:
            line = f.readline().strip()
            for sep in ["\t", ",", ";", " "]:
                cols = line.split(sep)
                if all(rc in cols for rc in self.required_columns):
                    self.separator = sep
                    break
            else:
                raise InvalidPeprecError(
                    "Could not infer separator. Please validate the PeptideRecord "
                    "header and/or the `required_columns` setting."
                )

    def _validate_required_columns(self) -> None:
        """Raise InvalidPeprecError if not all required columns are present."""
        with open(self.filename, "rt") as f:
            reader = csv.reader(f, delimiter=self.separator)
            columns = next(reader)
            for rc in self.required_columns:
                if rc not in columns:
                    raise InvalidPeprecError(f"Required column missing: `{rc}`")

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


class PeptideRecordReader(ReaderBase):
    def __init__(
        self,
        filename: Union[str, Path],
        modification_definitions: list[dict[str, str]] = None,
    ) -> None:
        """
        Reader for Peptide Record PSM files.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to PSM file.
        modification_definitions: list[dict[str, str]]
            List of modification definition dictionaries. See
            `Parsing of modification labels`_ for more information.

        Examples
        --------

        PeptideRecordReader supports iteration:

        >>> from psm_utils.io.peptide_record import PeptideRecordReader
        >>> for psm in PeptideRecordReader("peprec.txt"):
        ...     print(psm.peptide.proforma)
        ACDEK
        AC[Carbamidomethyl]DEFGR
        [Acetyl]-AC[Carbamidomethyl]DEFGHIK

        Or a full file can be read at once into a :py:class:`psm_utils.psm_list.PSMList`
        object:

        >>> peprec_reader = PeptideRecordReader("peprec.txt")
        >>> psm_list = peprec_reader.read_file()

        """

        super().__init__(filename, modification_definitions)
        self._peprec = PeptideRecord(self.filename)

        # Define named tuple for single PeptideRecord entries, based on
        # configured columns
        columns = self._peprec.required_columns + self._peprec.optional_columns
        self.PeprecEntry = namedtuple(
            "PeprecEntry", columns, defaults=[None for _ in columns]
        )

    def __iter__(self) -> Iterable[PeptideSpectrumMatch]:
        """Iterate over file and return PSMs one-by-one."""
        with open(self.filename, "rt") as open_file:
            reader = csv.DictReader(open_file, delimiter=self._peprec.separator)
            for row in reader:
                entry = self.PeprecEntry(**row)
                psm = self._peprec_entry_to_psm(entry)
                yield psm

    def read_file(self) -> PSMList:
        """Read full PeptideRecord PSM file into a PSMList object."""
        psm_list = []
        with open(self.filename) as peprec_in:
            reader = csv.DictReader(peprec_in, delimiter=self._peprec.separator)
            for row in reader:
                entry = self.PeprecEntry(**row)
                psm_list.append(self._peprec_entry_to_psm(entry))
        return PSMList(psm_list)

    def _peprec_entry_to_psm(self, entry: NamedTuple) -> PeptideSpectrumMatch:
        """Parse single PeptideRecord entry to a `PeptideSpectrumMatch`."""
        # Parse sequence and modifications
        proforma = PeptideRecord.peprec_to_proforma(
            entry.peptide, entry.modifications, label_mapping=self._label_mapping
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

        return PeptideSpectrumMatch(
            peptide=proforma,
            spectrum_id=entry.spec_id,
            precursor_charge=entry.charge,
            is_decoy=is_decoy,
            retention_time=entry.observed_retention_time,
            score=entry.score,
            source="PeptideRecord",
            provenance_data={"peprec_filename": self.filename},
        )


class PeptideRecordWriter(WriterBase):
    # TODO Implement
    pass


class InvalidPeprecError(Exception):
    """Invalid PEPREC file."""

    pass


class _PeptideRecordOld:
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
