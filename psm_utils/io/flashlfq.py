"""
Reader and writer for the FlashLFQ generic TSV format.

See the `FlashLFQ documentation <https://github.com/smith-chem-wisc/FlashLFQ/wiki/Identification-Input-Formats>`_
for more information on the format.

Notes
-----
- The FlashLFQ format does not contain the actual spectrum identifier. When reading a FlashLFQ
  file, the spectrum identifier is set to the row number in the file.
- The FlashLFQ format does not contain the precursor m/z, but the theoretical monoisotopic mass.
  This value is not read into the PSM object, but can be calculated from the peptidoform.
- To read from a FlashLFQ file, the ``Full Sequence`` column is expected to contain a ProForma v2
  compatible peptidoform notation.

"""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Optional, Union

import numpy as np

from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io._utils import set_csv_field_size_limit
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

set_csv_field_size_limit()

LOGGER = logging.getLogger(__name__)


class FlashLFQReader(ReaderBase):
    """Reader for FlashLFQ TSV format."""

    required_columns = ["Full Sequence", "Precursor Charge"]

    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""
        with open(self.filename, "rt") as open_file:
            reader = csv.DictReader(open_file, delimiter="\t")
            if not all(col in reader.fieldnames for col in self.required_columns):
                raise PSMUtilsIOException(
                    f"FlashLFQ TSV file must contain the following columns: {self.required_columns}"
                )
            for i, row in enumerate(reader):
                yield self._parse_entry(row, spectrum_id=str(i))

    def _parse_entry(self, entry: dict, spectrum_id) -> PSM:
        """Parse single FlashLFQ TSV entry to :py:class:`~psm_utils.psm.PSM`."""
        # Replace empty strings with None
        entry = {k: v if v else None for k, v in entry.items()}

        # Parse entry
        return PSM(
            peptidoform=f"{entry['Full Sequence']}/{entry['Precursor Charge']}",
            spectrum_id=spectrum_id,
            run=entry.get("File Name"),
            retention_time=entry.get("Scan Retention Time"),
            protein_list=self._parse_protein_list(entry.get("Protein Accessions")),
        )

    @staticmethod
    def _parse_protein_list(protein_accessions: Optional[str]) -> list[str]:
        """Parse protein list string to list of protein accessions."""
        if not protein_accessions:
            return []
        elif ";" in protein_accessions:  # Docs define separator as semicolon
            return protein_accessions.split(";")
        elif "|" in protein_accessions:  # Example file uses pipe
            return protein_accessions.split("|")
        else:
            return [protein_accessions]  # Single protein


class FlashLFQWriter(WriterBase):
    """Reader for FlashLFQ TSV format."""

    def __init__(
        self,
        filename: Union[str, Path],
        *args,
        fdr_threshold: float = 0.01,
        only_targets: bool = True,
        **kwargs,
    ):
        """
        Reader for psm_utils TSV format.

        Parameters
        ----------
        filename
            Path to PSM file.
        fdr_threshold
            FDR threshold for filtering PSMs.
        only_targets
            If True, only target PSMs are written to file. If False, both target and decoy PSMs
            are written.

        """
        super().__init__(filename, *args, **kwargs)

        self.fdr_threshold = fdr_threshold
        self.only_targets = only_targets

        self._open_file = None
        self._writer = None
        self.fieldnames = None

    def __enter__(self) -> FlashLFQWriter:
        if Path(self.filename).is_file():
            # Get fieldnames from existing file
            with open(self.filename, "rt") as open_file:
                # Get fieldnames
                self.fieldnames = open_file.readline().strip().split("\t")
            mode = "at"
        else:
            # Set default fieldnames
            self.fieldnames = [
                "File Name",
                "Base Sequence",
                "Full Sequence",
                "Peptide Monoisotope Mass",
                "Scan Retention Time",
                "Precursor Charge",
                "Protein Accessions",
            ]
            mode = "wt"

        # Open file and writer
        self._open_file = open(self.filename, mode, newline="")
        self._writer = csv.DictWriter(
            self._open_file,
            fieldnames=self.fieldnames,
            extrasaction="ignore",
            delimiter="\t",
        )

        if mode == "wt":
            self._writer.writeheader()

        return self

    def __exit__(self, *args, **kwargs) -> None:
        self._open_file.close()
        self._open_file = None
        self._writer = None

    def write_psm(self, psm: PSM):
        """
        Write a single PSM to new or existing PSM file.

        Parameters
        ----------
        psm
            PSM object to write.

        """
        if psm.qvalue and psm.qvalue > self.fdr_threshold:
            return
        if self.only_targets and psm.is_decoy:
            return

        entry = self._psm_to_entry(psm)
        try:
            self._writer.writerow(entry)
        except AttributeError as e:
            raise PSMUtilsIOException(
                f"`write_psm` method can only be called if `{self.__class__.__qualname__}`"
                "is opened in context (i.e., using the `with` statement)."
            ) from e

    def write_file(self, psm_list: PSMList):
        """
        Write an entire PSMList to a new PSM file.

        Parameters
        ----------
        psm_list
            PSMList object to write to file.

        """
        # Filter out decoys
        if self.only_targets:
            # Accept both None and False
            target_mask = np.array([not psm.is_decoy for psm in psm_list])
            LOGGER.debug(f"Skipping {~target_mask.sum()} decoy PSMs for FlashLFQ file.")
        else:
            target_mask = np.ones(len(psm_list), dtype=bool)

        # Filter out PSMs above FDR threshold
        if any(psm.qvalue is None for psm in psm_list):
            LOGGER.warning(
                "Not all PSMs have a q-value. Skipping FDR filtering for FlashLFQ file."
            )
            fdr_mask = np.ones(len(psm_list), dtype=bool)
        else:
            fdr_mask = psm_list["qvalue"] <= self.fdr_threshold
        filtered_by_fdr = (~fdr_mask & target_mask).sum()
        LOGGER.debug(f"Skipping {filtered_by_fdr} PSMs above FDR threshold for FlashLFQ file.")

        filtered_psm_list = psm_list[target_mask & fdr_mask]

        with open(self.filename, "wt", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=self.fieldnames, delimiter="\t", extrasaction="ignore"
            )
            writer.writeheader()
            for psm in filtered_psm_list:
                writer.writerow(self._psm_to_entry(psm))

    @staticmethod
    def _psm_to_entry(psm: PSM) -> dict:
        """Convert :py:class:`~psm_utils.psm.PSM` to FlashLFQ TSV entry."""
        return {
            "File Name": psm.run,
            "Base Sequence": psm.peptidoform.sequence,
            "Full Sequence": psm.peptidoform.modified_sequence,
            "Peptide Monoisotope Mass": psm.peptidoform.theoretical_mass,
            "Scan Retention Time": psm.retention_time,
            "Precursor Charge": psm.peptidoform.precursor_charge,
            "Protein Accessions": ";".join(psm.protein_list),
        }
