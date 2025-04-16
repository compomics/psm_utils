"""
Reader and writers for Percolator Tab ``PIN``/``POUT`` PSM files.

The tab-delimited input and output format for Percolator are defined on the
`Percolator GitHub Wiki pages <https://github.com/percolator/percolator/wiki/Interface>`_.

Notes
-----

* While :py:class:`PercolatorTabReader` supports reading the peptide notation with
  preceding and following amino acids (e.g. ``R.ACDEK.F``), these amino acids are not stored and
  are not written by :py:class:`PercolatorTabWriter`.

"""

from __future__ import annotations

import csv
import gzip
import logging
import re
from pathlib import Path
from typing import Iterable, List, Optional, Tuple, Union

from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io._utils import set_csv_field_size_limit
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

LOGGER = logging.getLogger(__name__)
set_csv_field_size_limit()


class PercolatorTabReader(ReaderBase):
    def __init__(
        self,
        filename: str | Path,
        score_column=None,
        retention_time_column=None,
        mz_column=None,
        *args,
        **kwargs,
    ) -> None:
        """
        Reader for Percolator Tab PIN/POUT PSM file.

        As the score, retention time, and precursor m/z are often embedded as feature
        columns, but not with a fixed column name, their respective column names need to
        be provided as parameters to the class. If not provided, these properties will
        not be added to the resulting :py:class:`~psm_utils.psm.PSM`.
        Nevertheless, they will still be added to its
        :py:attr:`~psm_utils.psm.PSM.rescoring_features` property
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
        super().__init__(filename, *args, **kwargs)
        self.score_column = score_column
        self.rt_column = retention_time_column
        self.mz_column = mz_column

        # Constant properties
        self._protein_separator = "|||"
        self.non_feature_columns = [
            "specid",
            "psmid",
            "label",
            "scannr",
            "peptide",
            "proteins",
            "q-value",
            "posterior_error_prob",
            "proteinids",
        ]

        # Properties derived from header
        self.fieldnames = self._read_header(self.filename)
        self.id_column = self.fieldnames[0]
        self.charge_column, self.charge_onehot_columns = self._infer_charge_columns(
            self.fieldnames
        )

        # Check if `score` in header if score_column was not defined (default in POUT)
        if not score_column and "score" in self.fieldnames:
            self.score_column = "score"

        # Validate column names from parameters
        for col in [self.score_column, self.rt_column, self.mz_column]:
            if col and col.lower() not in self.fieldnames:
                raise PercolatorIOException(
                    f"Column `{col}` not found in header of Percolator Tab file "
                    f"`{self.filename}`."
                )

    def __iter__(self) -> Iterable[PSM]:
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

    @staticmethod
    def _read_header(filename):
        if str(filename).endswith(".gz"):
            with gzip.open(filename, "rt") as f:
                fieldnames = f.readline().strip().lower().split("\t")
        else:
            with open(filename, "rt") as f:
                fieldnames = f.readline().strip().lower().split("\t")
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
        # Remove leading and trailing amino acids (e.g., R.PEPTIDE.S -> PEPTIDE)
        match = re.match(r"^(?:[A-Z-])?\.(.+)\.(?:[A-Z-])?$", percolator_peptide)
        peptidoform = match[1] if match else percolator_peptide
        # Handle Comet's n-terminal modification format: n[42.0106]PEPTIDE
        peptidoform = re.sub(r'^n\[([+-]?[\w\.]*?)\]', r'[\1]-', peptidoform)
        
        # Ensure positive values inside square brackets have a '+' sign
        peptidoform = re.sub(r'\[(\d+[\.]*\d*)]', r'[+\1]', peptidoform)
        
        if charge:
            peptidoform += f"/{charge}"
        return Peptidoform(peptidoform)

    def _parse_charge(self, entry) -> int | None:
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
        """Parse Percolator TSV entry to PSM."""
        label = entry["label"] if "label" in entry else None
        is_decoy = True if label == "-1" else False if label == "1" else None
        rescoring_features = {
            k: str(v) for k, v in entry.items() if k not in self.non_feature_columns
        }
        charge = self._parse_charge(entry)
        peptidoform = self._parse_peptidoform(entry["peptide"], charge)
        protein_list = (
            entry["proteins"].split(self._protein_separator)
            if "proteins" in entry
            else entry["proteinids"].split(self._protein_separator)
            if "proteinids" in entry
            else None
        )
        psm = PSM(
            peptidoform=peptidoform,
            spectrum_id=entry[self.id_column],
            is_decoy=is_decoy,
            score=float(entry[self.score_column.lower()]) if self.score_column else None,
            qvalue=entry["q-value"] if "q-value" in entry else None,
            pep=entry["posterior_error_prob"] if "posterior_error_prob" in entry else None,
            precursor_mz=float(entry[self.mz_column.lower()]) if self.mz_column else None,
            retention_time=float(entry[self.rt_column.lower()]) if self.rt_column else None,
            protein_list=protein_list,
            source="percolator",
            provenance_data={"filename": str(self.filename)},
            rescoring_features=rescoring_features,
        )
        return psm


class PercolatorTabWriter(WriterBase):
    def __init__(
        self,
        filename: str | Path,
        style: Optional[str] = None,
        feature_names: Optional[list[str]] = None,
        add_basic_features: bool = False,
        *args,
        **kwargs,
    ) -> None:
        """
        Writer for Percolator TSV "PIN" and "POUT" PSM files.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to PSM file.
        style: str, optional
            Percolator Tab style. One of {``pin``, ``pout``}. If ``pin``, the columns
            ``SpecId``, ``Label``, ``ScanNr``, ``ChargeN``, ``PSMScore``, ``Peptide``, and
            ``Proteins`` are written alongside the requested feature names
            (see ``feature_names``). If ``pout``, the columns ``PSMId``, ``Label``, ``score``,
            ``q-value``, ``posterior_error_prob``, ``peptide``, and ``proteinIds`` are written.
            By default, the style is inferred from the file name extension.
        feature_names: list[str], optional
            List of feature names to extract from PSMs and write to file. List values
            should correspond to keys in the
            :py:class:`~psm_utils.psm.PSM.rescoring_features` property.
            If :py:const:`None`, no rescoring features will be written to the file. If appending to
            an existing file, the existing header will be used to determine the feature
            names. Only has effect with ``pin`` style.
        add_basic_features: bool, optional
            If :py:const:`True`, add ``PSMScore`` and ``ChargeN`` features to the file. Only has
            effect with ``pin`` style. Default is :py:const:`False`.

        """
        super().__init__(filename, *args, **kwargs)
        self.feature_names = list(feature_names) if feature_names else []
        self.add_basic_features = add_basic_features

        if not style:
            suffix = self.filename.suffix.lower()
            if suffix == ".pin":
                self.style = "pin"
            elif suffix == ".pout":
                self.style = "pout"
            else:
                raise PercolatorIOException(
                    f"Could not infer Percolator Tab style from file extension `{suffix}`. "
                    "Please provide the `style` parameter."
                )
        else:
            self.style = style

        if self.style == "pin":
            basic_features = ["PSMScore", "ChargeN"] if add_basic_features else []
            self._columns = (
                ["SpecId", "Label", "ScanNr"]
                + basic_features
                + self.feature_names
                + ["Peptide", "Proteins"]
            )
        elif self.style == "pout":
            self._columns = [
                "PSMId",
                "Label",
                "score",
                "q-value",
                "posterior_error_prob",
                "peptide",
                "proteinIds",
            ]
        else:
            raise ValueError("Invalid Percolator Tab style. Should be one of {`pin`, `pout`}.")

        self._open_file = None
        self._writer = None
        self._protein_separator = "|||"
        self._current_scannr = 0

    def __enter__(self) -> PercolatorTabWriter:
        """Either open existing file in append mode or new file in write mode."""
        file_existed = self.filename.is_file()
        mode = "at" if file_existed else "wt"
        self._open_file = _PercolatorTabIO(
            self.filename, mode, newline="", protein_separator=self._protein_separator
        )
        if file_existed:
            fieldnames, self._current_scannr = self._parse_existing_file(self.filename)
        else:
            fieldnames = self._columns
            self._current_scannr = -1
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
        self._current_scannr = None

    def write_psm(self, psm: PSM):
        """Write a single PSM to the PSM file."""
        entry = self._psm_to_entry(psm)
        self._current_scannr += 1
        entry["ScanNr"] = self._current_scannr
        try:
            self._writer.writerow(entry)
        except AttributeError as e:
            raise PSMUtilsIOException(
                f"`write_psm` method can only be called if `{self.__class__.__qualname__}`"
                "is opened in context (i.e., using the `with` statement)."
            ) from e
        else:
            self._current_scannr = entry["ScanNr"]

    def write_file(self, psm_list: PSMList):
        """Write an entire PSMList to the PSM file."""
        with _PercolatorTabIO(
            self.filename, "wt", newline="", protein_separator=self._protein_separator
        ) as f:
            writer = csv.DictWriter(
                f, fieldnames=self._columns, delimiter="\t", extrasaction="ignore"
            )
            writer.writeheader()
            for psm in psm_list:
                writer.writerow(self._psm_to_entry(psm))

    def _psm_to_entry(self, psm: PSM):
        """Parse PSM to Percolator Tab entry."""
        if self.style == "pin":
            entry = {
                "SpecId": psm.spectrum_id,
                "Label": None if psm.is_decoy is None else -1 if psm.is_decoy else 1,
                "Peptide": "." + re.sub(r"/\d+$", "", psm.peptidoform.proforma) + ".",
                "Proteins": self._protein_separator.join(psm.protein_list)
                if psm.protein_list
                else "PEP_" + psm.peptidoform.proforma,
            }
            if self.add_basic_features:
                entry.update({"ChargeN": psm.peptidoform.precursor_charge, "PSMScore": psm.score})
            try:
                entry.update(psm.rescoring_features)
            except TypeError:
                pass
        elif self.style == "pout":
            entry = {
                "PSMId": psm.spectrum_id,
                "Label": None if psm.is_decoy is None else -1 if psm.is_decoy else 1,
                "score": psm.score,
                "q-value": psm.qvalue,
                "posterior_error_prob": psm.pep,
                "peptide": psm.peptidoform.proforma,
                "proteinIds": self._protein_separator.join(psm.protein_list)
                if psm.protein_list
                else None,
            }
        return entry

    @staticmethod
    def _parse_existing_file(
        filename: Union[str, Path], style: str
    ) -> Tuple[List[str], Optional[int]]:
        """Parse existing Percolator Tab file to determine fieldnames and last ScanNr."""
        # Get fieldnames
        if str(filename).endswith(".gz"):
            open_func = gzip.open
        else:
            open_func = open
            
        with open_func(filename, "rt") as open_file:
            for line in open_file:
                fieldnames = line.strip().split("\t")
                break
            else:
                raise PercolatorIOException(
                    f"Existing file {filename} is not a valid Percolator Tab file."
                )
        if not _fieldnames_are_valid(fieldnames, style):
            raise PercolatorIOException(
                f"Existing file {filename} is not a valid Percolator Tab file of style {style}."
            )

        # Get last ScanNr
        last_scannr = None
        with open_func(filename, "rt") as open_file:
            # Read last line
            open_file.seek(0)
            last_line = None
            for line in open_file:
                if line.strip():
                    last_line = line
        if last_line:
            # Parse last line
            last_line_items = {k: v for k, v in zip(fieldnames, last_line.strip().split("\t"))}
            try:
                last_scannr = int(last_line_items["ScanNr"])
            except (KeyError, ValueError):
                pass

        if last_scannr is None:
            last_scannr = -1
            LOGGER.warning(
                f"Could not determine last ScanNr from file {filename}. Starting ScanNr from 0."
            )

        return fieldnames, last_scannr


class _PercolatorTabIO:
    def __init__(self, *args, protein_separator="|||", **kwargs) -> None:
        """File reader and writer for Percolator Tab files with fixed Proteins tab."""
        filename = args[0]
        if str(filename).endswith(".gz"):
            self._open_file = gzip.open(*args, **kwargs)
        else:
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
                yield line.lower()
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


def _fieldnames_are_valid(fieldnames: List[str], style: str) -> bool:
    """Check if fieldnames are valid for Percolator Tab style."""
    if style == "pin":
        required_columns = ["SpecId", "Label", "ScanNr"]
    elif style == "pout":
        required_columns = ["PSMId", "score", "q-value", "posterior_error_prob"]
    else:
        raise ValueError("Invalid Percolator Tab style. Should be one of {`pin`, `pout`}.")
    return all(col in fieldnames for col in required_columns)


def join_pout_files(
    target_filename: str | Path,
    decoy_filename: str | Path,
    output_filename: str | Path,
):
    """
    Join target and decoy Percolator Out (POUT) files into single PercolatorTab file.

    Parameters
    ----------
    target_filename: str, Path
    decoy_filename: str, Path
    output_filename: str, Path

    """
    target_reader = PercolatorTabReader(target_filename, score_column="score")
    decoy_reader = PercolatorTabReader(decoy_filename, score_column="score")
    with PercolatorTabWriter(output_filename, style="pout") as writer:
        for psm in target_reader:
            psm.is_decoy = False
            writer.write_psm(psm)
        for psm in decoy_reader:
            psm.is_decoy = True
            writer.write_psm(psm)


class PercolatorIOException(PSMUtilsIOException):
    """Exception for Percolator Tab file I/O errors."""

    pass
