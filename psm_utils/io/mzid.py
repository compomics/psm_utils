from __future__ import annotations

import logging
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Union

from pyteomics import mzid

from psm_utils.exceptions import PSMUtilsException
from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PeptideSpectrumMatch
from psm_utils.psm_list import PSMList

logger = logging.getLogger(__name__)


class MzidReader(ReaderBase):
    """Reader for MaxQuant msms.txt PSM files."""

    def __init__(self, filename: Union[str, Path]) -> None:
        super().__init__(filename)
        """
        Reader for MZID Record PSM files.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to PSM file.

        Examples
        --------

        MzidReader supports iteration:

        >>> from psm_utils.io.mzid import MzidReader
        >>> for psm in MzidReader("peptides_1_1_0.mzid"):
        ...     print(psm.peptide.proforma)
        ACDEK
        AC[Carbamidomethyl]DEFGR
        [Acetyl]-AC[Carbamidomethyl]DEFGHIK

        Or a full file can be read at once into a :py:class:`psm_utils.psm_list.PSMList`
        object:

        >>> mzid_reader = MzidReader("peprec.txt")
        >>> psm_list = mzid_reader.read_file()

        """
        self.filename = filename
        self.source = self._infer_source()
        self.searchengine_key_dict = self._get_searchengine_specific_keys()

    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""
        with mzid.read(self.filename) as reader:

            for spectrum in reader:
                spectrum_title = spectrum[self.searchengine_key_dict["spectrum_key"]]
                rawfile = self._get_rawfile_name(spectrum["location"])
                for psm_dict in spectrum["SpectrumIdentificationItem"]:
                    psm = self._get_peptide_spectrum_match(
                        spectrum_title, rawfile, psm_dict
                    )
                    yield psm

    def __next__(self) -> PeptideSpectrumMatch:
        return super().__next__()

    def read_file(self) -> PSMList:
        """Read full mzid file to PSM list object"""

        with mzid.read(self.filename) as reader:
            # Go over each Spectrum that has a match
            psm_list = []
            # Go over each PSM
            for spectrum in reader:
                spectrum_title = spectrum[self.searchengine_key_dict["spectrum_key"]]
                rawfile = self._get_rawfile_name(spectrum["location"])
                for psm in spectrum["SpectrumIdentificationItem"]:
                    psm_list.append(
                        self._get_peptide_spectrum_match(spectrum_title, rawfile, psm)
                    )

        return PSMList(psm_list)

    def _get_peptide_spectrum_match(
        self,
        spectrum_title: str,
        rawfile: str,
        SpectrumIdentificationItem: dict[str, Union[str, float, list]],
    ) -> PeptideSpectrumMatch:
        """Parse single mzid Record entry to a `PeptideSpectrumMatch`."""

        try:
            peptide = self._parse_peptidoform(
                SpectrumIdentificationItem["PeptideSequence"],
                SpectrumIdentificationItem["Modification"],
            )
        except KeyError:
            peptide = Peptidoform(SpectrumIdentificationItem["PeptideSequence"])

        isdecoy, protein_list = self._parse_PeptideEvidenceRef(
            SpectrumIdentificationItem["PeptideEvidenceRef"]
        )
        # retention time is often missing from mzid
        try:
            rt = float(SpectrumIdentificationItem[self.searchengine_key_dict["rt_key"]])
        except KeyError:
            rt = float("nan")

        psm = PeptideSpectrumMatch(
            peptide=peptide,
            spectrum_id=str(spectrum_title),
            run=str(rawfile),
            is_decoy=isdecoy,
            score=float(
                SpectrumIdentificationItem[self.searchengine_key_dict["score_key"]]
            ),
            precursor_charge=int(SpectrumIdentificationItem["chargeState"]),
            # Or use calculatedMassToCharge?
            precursor_mz=float(SpectrumIdentificationItem["experimentalMassToCharge"]),
            # Retention time often not included in mzid output
            retention_time=rt,
            protein_list=protein_list,
            source=self.source,
            provenance_data=({f"{self.source}_filename": self.filename}),
            metadata=self._get_searchengine_specific_metadata(
                SpectrumIdentificationItem
            ),
        )
        return psm

    @staticmethod
    def _parse_peptidoform(seq: str, modification_list: list[dict]):
        """Return a peptido form"""
        # TODO: not able to differentiate between c-terminal mod and mod on last aa

        modification_lengths = 0
        for mod in modification_list:
            if mod["location"] == 0:
                proforma_mod = f"[{mod['name']}]-"
                seq = proforma_mod + seq

            else:
                proforma_mod = f"[{mod['name']}]"
                seq = (
                    seq[0 : mod["location"] + modification_lengths]
                    + proforma_mod
                    + seq[mod["location"] + modification_lengths :]
                )

            modification_lengths += len(proforma_mod)

        return Peptidoform(seq)

    def _infer_source(self):
        """Get the source of the mzid file"""

        mzid_xml = ET.parse(self.filename)
        root = mzid_xml.getroot()
        name_space = self._get_xml_namespace(root.tag)

        return root.find(f".//{name_space}AnalysisSoftware").attrib["name"]

    @staticmethod
    def _get_xml_namespace(root_tag):
        """Get the namespace of the xml root"""

        m = re.match(r"\{.*\}", root_tag)
        return m.group(0) if m else ""

    @staticmethod
    def _parse_PeptideEvidenceRef(PeptideEvidenceList: list[dict]):
        """Parse peptide evidence References of PSM"""

        isdecoy = PeptideEvidenceList[0]["isDecoy"]

        protein_list = [
            d["accession"] for d in PeptideEvidenceList if "accession" in d.keys()
        ]

        return isdecoy, protein_list

    def _get_searchengine_specific_keys(self):
        """Get searchengine specific score"""
        # TODO works for PEAKS pro?
        if "PEAKS" in self.source:
            return {
                "score_key": "PEAKS:peptideScore",
                "rt_key": "retention time",
                "spectrum_key": "spectrumID",
            }

        elif self.source == "MS-GF+":
            return {
                "score_key": "MS-GF:RawScore",
                "rt_key": "scan start time",
                "spectrum_key": "spectrum title",
            }

    def _get_searchengine_specific_metadata(self, SpectrumIdentificationItem):
        """Get searchengine specific psm metadata"""

        metadata = {
            "calculatedMassToCharge": SpectrumIdentificationItem[
                "calculatedMassToCharge"
            ],
            "rank": SpectrumIdentificationItem["rank"],
            "Length": len(SpectrumIdentificationItem["PeptideSequence"]),
        }
        if self.source == "MS-GF+":
            metadata.update(
                {
                    "MS-GF:DeNovoScore": SpectrumIdentificationItem[
                        "MS-GF:DeNovoScore"
                    ],
                    "MS-GF:EValue": SpectrumIdentificationItem["MS-GF:EValue"],
                    "MS-GF:DeNovoScore": SpectrumIdentificationItem[
                        "MS-GF:DeNovoScore"
                    ],
                    "IsotopeError": SpectrumIdentificationItem["IsotopeError"],
                }
            )
        return metadata

    @staticmethod
    def _get_rawfile_name(file_location: str) -> str:
        """Get rawfile name out of mzid file location or filename"""

        filename = file_location.rsplit("/", 1)[1]
        return filename.rsplit(".", 1)[0]


class MzidWriter(WriterBase):
    """Writer for MaxQuant msms.txt files (not implemented)."""

    def __init__(self, filename):
        raise NotImplementedError()

    def write_psm(self, psm: PeptideSpectrumMatch):
        """Write a single PSM to the MaxQuant msms.txt PSM file."""
        raise NotImplementedError()

    def write_file(self, psm_list: PSMList):
        """Write entire PSMList to the MaxQuant msms.txt PSM file."""
        raise NotImplementedError()
