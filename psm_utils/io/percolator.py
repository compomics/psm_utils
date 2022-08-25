"""
Reader and writers for Percolator Tab ``PIN``/``POUT`` PSM files.

The tab-delimited input and output format for Percolator are defined on the
`Percolator GitHub Wiki pages <https://github.com/percolator/percolator/wiki/Interface>`_.

Notes
-----

* While :py:class:`PercolatorTabReader` supports reading the peptide notation with
  preceding and following amino acids (e.g. ``R.ACDEK.F``),
  :py:class:`PercolatorTabWriter` simply writes peptides in Proforma format, without
  preceding and following amino acids.


"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Iterable, Optional, Union

from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PeptideSpectrumMatch
from psm_utils.psm_list import PSMList


class PercolatorTabReader(ReaderBase):
    def __init__(
        self,
        filename: Union[str, Path],
        score_column=None,
        retention_time_column=None,
        mz_column=None,
    ) -> None:
        """
        Reader for Percolator Tab PIN/POUT PSM file.

        As the score, retention time, and precursor m/z are often embedded as feature
        columns, but not with a fixed column name, their respective column names need to
        be provided as parameters to the class. If not provided, these properties will
        not be added to the resulting :py:class:`~psm_utils.psm.PeptideSpectrumMatch`.
        Nevertheless, they will still be added to its
        :py:attr:`~psm_utils.psm.PeptideSpectrumMatch.rescoring_features` property
        dictionary, along with the other features.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to PSM file.
        score_column: str, optional
            Name of the column that holds the primary PSM score.
        retention_time_column: str, optional
            Name of the column that holds the retention time.
        mz_column: str, optional
            Name of the column that holds the precursor m/z.

        """
        super().__init__(filename)
        self.score_column = score_column
        self.rt_column = retention_time_column
        self.mz_column = mz_column

        # Constant properties
        self._protein_separator = "|||"
        self.non_feature_columns = [
            "SpecId",
            "PSMId",
            "Label",
            "ScanNr",
            "Peptide",
            "Proteins",
        ]

        # Properties derived from header
        self.fieldnames = self._read_header(self.filename)
        self.id_column = self.fieldnames[0]
        self.charge_column, self.charge_onehot_columns = self._infer_charge_columns(
            self.fieldnames
        )

        # Validate column names from parameters
        for col in [self.score_column, self.rt_column, self.mz_column]:
            if col and col not in self.fieldnames:
                raise ValueError(
                    f"Column `{col}` not found in header of Percolator Tab file "
                    f"`{self.filename}`."
                )

    def __iter__(self) -> Iterable[PeptideSpectrumMatch]:
        """Iterate over file and return PSMs one-by-one."""
        with _PercolatorTabIO(
            self.filename, "rt", protein_separator=self._protein_separator
        ) as open_file:
            reader = csv.DictReader(open_file, delimiter="\t")
            for entry in reader:
                if entry[self.id_column] == "DefaultDirection":
                    pass
                else:
                    psm = self._parse_entry(entry)
                    yield psm

    def read_file(self) -> PSMList:
        """Read full PSM file into a PSMList object."""
        psm_list = []
        for psm in self.__iter__():
            psm_list.append(psm)
        return PSMList(psm_list)

    @staticmethod
    def _read_header(filename):
        with open(filename, "rt") as f:
            fieldnames = f.readline().strip().split("\t")
        return fieldnames

    @staticmethod
    def _infer_charge_columns(fieldnames):
        """Infer columns that hold the precursor charge from the header fieldnames."""
        # Infer single charge column
        charge_column = None
        for col in ["charge", "Charge"]:
            if col in fieldnames:
                charge_column = col

        # Infer one-hot encoded charge columns
        charge_onehot_columns = {}
        for col in fieldnames:
            match = re.fullmatch("(charge|Charge)([0-9]+)", col)
            if match:
                # charge state -> column name mapping
                charge_onehot_columns[int(match[2])] = match[0]

        return charge_column, charge_onehot_columns

    @staticmethod
    def _parse_peptidoform(percolator_peptide, charge):
        """Parse Percolator TSV peptide notation to Peptidoform."""
        # Remove leading and trailing amino acids
        match = re.match("^[A-Z-]\.(.+)\.[A-Z-]$", percolator_peptide)
        peptidoform = match[1] if match else percolator_peptide
        if charge:
            peptidoform += f"/{charge}"
        return Peptidoform(peptidoform)

    def _parse_charge(self, entry) -> Union[int, None]:
        """Parse charge state from single or one-hot encoded charge state."""
        if self.charge_column:
            return int(entry["charge"])
        elif self.charge_onehot_columns:
            for charge_state, column_name in self.charge_onehot_columns.items():
                if entry[column_name] == "1":
                    return charge_state
        else:
            return None

    def _parse_entry(self, entry):
        """Parse Percolator TSV entry to PeptideSpectrumMatch."""
        label = entry["Label"]
        is_decoy = True if label == "-1" else False if label == "1" else None
        rescoring_features = {
            k: v for k, v in entry.items() if k not in self.non_feature_columns
        }
        charge = self._parse_charge(entry)
        peptide = self._parse_peptidoform(entry["Peptide"], charge)
        psm = PeptideSpectrumMatch(
            peptide=peptide,
            spectrum_id=entry[self.id_column],
            is_decoy=is_decoy,
            score=float(entry[self.score_column]) if self.score_column else None,
            precursor_mz=float(entry[self.mz_column]) if self.mz_column else None,
            retention_time=float(entry[self.rt_column]) if self.rt_column else None,
            protein_list=entry["Proteins"].split(self._protein_separator),
            source="percolator",
            provenance_data={"filename": str(self.filename)},
            rescoring_features=rescoring_features,
        )
        return psm


class PercolatorTabWriter(WriterBase):
    def __init__(
        self, filename: Union[str, Path], feature_names: Optional[list[str]] = None
    ) -> None:
        """
        Writer for Percolator TSV "PIN" and "POUT" PSM files.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to PSM file.
        feature_names: list[str], optional
            List of feature names to extract from PSMs and write to file. List values
            should correspond to keys in the
            :py:class:`~psm_utils.psm.PeptideSpectrumMatch.rescoring_features` property.
            If :py:const:`None`, no rescoring features will be written to the file. If appending to
            an existing file, the existing header will be used to determine the feature
            names.

        """
        super().__init__(filename)
        self.feature_names = list(feature_names) if feature_names else []

        self._columns = (
            ["PSMId", "Label", "ScanNr"] + self.feature_names + ["Peptide", "Proteins"]
        )
        self._open_file = None
        self._writer = None
        self._protein_separator = "|||"

    def __enter__(self) -> PercolatorTabWriter:
        """Either open existing file in append mode or new file in write mode."""
        file_existed = self.filename.is_file()
        mode = "at" if file_existed else "wt"
        self._open_file = _PercolatorTabIO(
            self.filename, mode, newline="", protein_separator=self._protein_separator
        )
        if file_existed:
            with open(self.filename, "rt") as open_file:
                for line in open_file:
                    fieldnames = line.strip().split("\t")
                    break
                else:
                    raise ValueError(
                        f"File {self.filename} exists, but is not a valid Percolator Tab file."
                    )
        else:
            fieldnames = self._columns
        self._writer = csv.DictWriter(
            self._open_file,
            fieldnames=fieldnames,
            extrasaction="ignore",
            delimiter="\t",
        )
        if not file_existed:
            self._writer.writeheader()
        return self

    def __exit__(self, *args, **kwargs) -> None:
        self._open_file.close()
        self._open_file = None
        self._writer = None

    def write_psm(self, psm: PeptideSpectrumMatch):
        """Write a single PSM to the PSM file."""
        try:
            self._writer.writerow(self._psm_to_entry(psm))
        except AttributeError:
            raise PSMUtilsIOException(
                "`write_psm` method can only be called if `PercolatorTabWriter` is "
                "opened in context (i.e., using the `with` statement)."
            )

    def write_file(self, psm_list: PSMList):
        """Write an entire PSMList to the PSM file."""
        with open(self.filename, "wt", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=self._columns, delimiter="\t")
            writer.writeheader()
            for psm in psm_list:
                writer.writerow(self._psm_to_entry(psm))

    def _psm_to_entry(self, psm: PeptideSpectrumMatch):
        """Parse PeptideSpectrumMatch to Percolator Tab entry."""
        entry = {
            "PSMId": psm.spectrum_id,
            "Label": None if psm.is_decoy is None else -1 if psm.is_decoy else 1,
            "ScanNr": None,  # TODO
            "Peptide": psm.peptide.proforma,
            "Proteins": self._protein_separator.join(psm.protein_list),
        }
        try:
            entry.update(psm.rescoring_features)
        except TypeError:
            pass
        return entry


class _PercolatorTabIO:
    def __init__(self, *args, protein_separator="|||", **kwargs) -> None:
        """File reader and writer for Percolator Tab files with fixed Proteins tab."""
        self._open_file = open(*args, **kwargs)
        self.protein_separator = protein_separator

    def __enter__(self, *args, **kwargs) -> _PercolatorTabIO:
        return self

    def __exit__(self, *args, **kwargs):
        self.close()

    def __iter__(self):
        """Yield lines in file with Proteins tab replaced by separator."""
        number_of_columns = None
        for i, line in enumerate(self._open_file):
            if i == 0:
                number_of_columns = len(line.split("\t"))
                yield line
            elif i == 1 & line.startswith("DefaultDirection"):
                yield line
            else:
                r = line.strip().split("\t")
                row_columns = r[: number_of_columns - 1]
                row_proteins = r[number_of_columns - 1 :]
                row_columns.append(self.protein_separator.join(row_proteins))
                line = "\t".join(row_columns) + "\n"
                yield line

    def close(self):
        self._open_file.close()

    def write(self, __s: str):
        """Write line to file with Proteins separator replaced by tab."""
        __s = __s.replace(self.protein_separator, "\t")
        self._open_file.write(__s)
