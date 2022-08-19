"""
Reader and writer for a simple, lossless psm_utils TSV format.

Most PSM file formats will introduce a loss of some information when reading,
writing, or converting with :py:mod:`psm_utils.io` due to differences between file
formats. In contrast, :py:class:`psm_utils.io.psm_list.PSMList` objects can be written
to — or read from — this simple TSV format without any information loss (with exception
of the free-form :py:attr:`spectrum` attribute).

The format follows basic TSV rules, using tab as delimiter, and supports quoting when
a field contains the delimiter. Peptidoforms are written in the `HUPO-PSI ProForma 2.0
<https://psidev.info/proforma>`_ notation.

Required and optional columns equate to the required and optional attributes of
:py:class:`psm_utils.psm.PeptideSpectrumMatch`. Dictionary items in
:py:attr:`provenance_data`, :py:attr:`metadata`, and :py:attr:`rescoring_features`
are flattened to separate columns, each with their column names prefixed with
``provenance:``, ``meta:``, and ``rescoring:``, respectively.


**Examples**

.. code-block::
    :caption: Minimal :py:mod:`psm_utils` TSV file

    peptide	spectrum_id
    RNVIDKVAK/2	1
    KHLEQHPK/2	2
    ...

.. code-block::
    :caption: Recommended :py:mod:`psm_utils` TSV file, compatible with `HUPO-PSI Universal Spectrum Identifier <https://www.psidev.info/usi>`_

    peptide	spectrum_id	run	collection
    VLHPLEGAVVIIFK/2	17555	Adult_Frontalcortex_bRP_Elite_85_f09	PXD000561
    ...

.. code-block::
    :caption: Full :py:mod:`psm_utils` TSV file, converted from a Percolator Tab file

    peptide	spectrum_id	run	collection	spectrum	is_decoy	score	precursor_mz	retention_time	protein_list	source	provenance:filename	rescoring:ExpMass	rescoring:CalcMass	rescoring:hyperscore	rescoring:deltaScore	rescoring:frac_ion_b	rescoring:frac_ion_y	rescoring:Mass	rescoring:dM	rescoring:absdM	rescoring:PepLen	rescoring:Charge2	rescoring:Charge3	rescoring:Charge4	rescoring:enzN	rescoring:enzC	rescoring:enzInt
    RNVIDKVAK/2	_3_2_1				False	20.3	1042.64		['DECOY_sp|Q8U0H4_REVERSED|RTCB_PYRFU-tRNA-splicing-ligase-RtcB-OS=Pyrococcus-furiosus...']	percolator	pyro.t.xml.pin	1042.64	1042.64	20.3	6.6	0.444444	0.333333	1042.64	0.0003	0.0003	9	1	0	0	1	0	1
    KHLEQHPK/2	_4_2_1				False	26.5	1016.56		['sp|Q8TZD9|RS15_PYRFU-30S-ribosomal-protein-S15-OS=Pyrococcus-furiosus-(strain-ATCC...']	percolator	pyro.t.xml.pin	1016.56	1016.56	26.5	18.5	0.375	0.75	1016.56	0.001	0.001	8	1	0	0	1	0	0
    ...


"""
from __future__ import annotations

import ast
import csv
import dataclasses
from pathlib import Path
from typing import Optional, Union

from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.psm import PeptideSpectrumMatch
from psm_utils.psm_list import PSMList


class TSVReader(ReaderBase):
    """Reader for psm_utils TSV format."""

    def __init__(self, filename: Union[str, Path]) -> None:
        super().__init__(filename)

    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""
        with open(self.filename, "rt") as open_file:
            reader = csv.DictReader(open_file, delimiter="\t")
            for row in reader:
                yield PeptideSpectrumMatch(**self._parse_entry(row))

    def read_file(self) -> PSMList:
        """Read full PSM file into a PSMList object."""
        return PSMList([psm for psm in self.__iter__()])

    @staticmethod
    def _parse_entry(entry: dict):
        """Parse single TSV entry to :py:class:`psm_utils.psm.PeptideSpectrumMatch`."""
        # Replace empty strings with None
        entry = {k: v if v else None for k, v in entry.items()}

        # Parse protein list
        entry["protein_list"] = ast.literal_eval(entry["protein_list"])

        # Extract dict properties
        parsed_entry = {}
        provenance_data = {}
        metadata = {}
        rescoring_features = {}
        for k, v in entry.items():
            if k.startswith("provenance:"):
                provenance_data[k[11:]] = v
            elif k.startswith("meta:"):
                metadata[k[5:]] = v
            elif k.startswith("rescoring:"):
                rescoring_features[k[10:]] = v
            else:
                parsed_entry[k] = v

        parsed_entry.update(
            {
                "provenance_data": provenance_data,
                "metadata": metadata,
                "rescoring_features": rescoring_features,
            }
        )

        return parsed_entry


class TSVWriter(WriterBase):
    """Reader for psm_utils TSV format."""

    def __init__(self, filename: Union[str, Path], example_psm: Optional[PeptideSpectrumMatch] = None):
        """
        Reader for psm_utils TSV format.

        Parameters
        ----------
        filename: str, Pathlib.Path
            Path to PSM file.
        example_psm: psm_utils.psm.PeptideSpectrumMatch, optional
            Example PSM, required to extract the column names when writing to a new
            file. Should contain all fields that are to be written to the PSM file,
            i.e., all items in the :py:attr:`provenance_data`, :py:attr:`metadata`, and
            :py:attr:`rescoring_features` attributes. In other words, items that are
            not present in the example PSM will not be written to the file, even though
            they are present in other PSMs passed to :py:meth:`write_psm` or
            :py:meth:`write_file`.
        """
        super().__init__(filename)

        self._open_file = None
        self._writer = None

        if example_psm:
            self.fieldnames = self._psm_to_entry(example_psm).keys()
        else:
            self.fieldnames = None

    def __enter__(self) -> TSVWriter:
        if Path(self.filename).is_file():
            with open(self.filename, "rt") as open_file:
                # Get fieldnames
                self.fieldnames = open_file.readline().strip().split("\t")
                # Check if newline is needed
                for line in open_file:
                    pass
                ends_on_newline = line[-1] in ["\n", "\r"]


            self._open_file = open(self.filename, "at", newline="")
            if not ends_on_newline:
                self._open_file.write()


            self._writer = csv.DictWriter(
                self._open_file,
                fieldnames=self.fieldnames,
                extrasaction="ignore",
                delimiter="\t",
            )
        else:
            if not self.fieldnames:
                raise ValueError("`example_psm` required when writing to new file.")
            self._open_file = open(self.filename, "wt", newline="")
            self._writer = csv.DictWriter(
                self._open_file,
                fieldnames=self.fieldnames,
                extrasaction="ignore",
                delimiter="\t",
            )
            self._writer.writeheader()
        return self

    def __exit__(self, *args, **kwargs) -> None:
        self._open_file.close()
        self._open_file = None
        self._writer = None

    def write_psm(self, psm: PeptideSpectrumMatch):
        """
        Write a single PSM to new or existing PSM file.

        Parameters
        ----------
        psm: PeptideSpectrumMatch
            PeptideSpectrumMatch object to write.

        """
        if not self._open_file:
            raise PSMUtilsIOException(
                f"`write_psm` method can only be called if `{self.__class__.__qualname__}`"
                "is opened in context (i.e., using the `with` statement)."
            )
        self._writer.writerow(self._psm_to_entry(psm))

    def write_file(self, psm_list: PSMList):
        """
        Write an entire PSMList to a new PSM file.

        Parameters
        ----------
        psm_list: PSMList
            PSMList object to write to file.

        """
        if not self.fieldnames:
            raise ValueError("`example_psm` required when writing to new file.")
        with open(self.filename, "wt", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=self.fieldnames, delimiter="\t")
            writer.writeheader()
            for psm in psm_list:
                writer.writerow(self._psm_to_entry(psm))

    @staticmethod
    def _psm_to_entry(psm: PeptideSpectrumMatch) -> dict:
        entry = dataclasses.asdict(psm)

        # Convert Peptidoform to proforma sequence
        entry["peptide"] = entry["peptide"].proforma

        # Drop spectrum
        del entry["spectrum"]

        # Flatten dictionary items
        if entry["provenance_data"]:
            entry.update(
                {"provenance:" + k: v for k, v in entry["provenance_data"].items()}
            )
        if entry["metadata"]:
            entry.update({"meta:" + k: v for k, v in entry["metadata"].items()})
        if entry["rescoring_features"]:
            entry.update(
                {"rescoring:" + k: v for k, v in entry["rescoring_features"].items()}
            )
        del entry["provenance_data"]
        del entry["metadata"]
        del entry["rescoring_features"]

        return entry
