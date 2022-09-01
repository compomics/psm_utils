from __future__ import annotations

import logging
import re
import xml.etree.ElementTree as ET
from multiprocessing.sharedctypes import Value
from pathlib import Path
from tkinter import scrolledtext
from typing import Union

from pyteomics import mzid

from psm_utils.exceptions import PSMUtilsException
from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PeptideSpectrumMatch
from psm_utils.psm_list import PSMList

logger = logging.getLogger(__name__)

STANDARD_SEARCHENGINE_SCORES = [
    "Amanda:AmandaScore",
    "Andromeda:score",
    "Byonic:Score",
    "Comet:Xcorr",
    "DeBunker:score",
    "IdentityE Score",
    "KSDP score",
    "MS-GF:RawScore",
    "MSFit:Mowse score",
    "MSPathFinder:RawScore" "MSPepSearch:score",
    "Mascot:score",
    "MetaMorpheus:score",
    "OMSSA:evalue",
    "OpenPepXL:score",
    "PEAKS:peptideScore",
    "PeptideShaker PSM score",
    "Phenyx:Pepzscore",
    "ProLuCID:xcorr",
    "ProSight:specral C-score",
    "Profound:z value",
    "ProteinProspector:score",
    "ProteinScape:SequestMetaScore",
    "ProteomeDiscoverer:Delta Score",
    "SEQUEST:xcorr",
    "SIM-XL score ",
    "SQID:score ",
    "Sonar:Score",
    "SpectrumMill:Score",
    "TopMG:spectral E-Value",
    "X!Tandem:hyperscore",
    "ZCore:probScore:",
    "Percolator:scrolledtext",
    "xi:score",
]


class MzidReader(ReaderBase):
    def __init__(self, filename: Union[str, Path]) -> None:

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

        >>> mzid_reader = MzidReader("peptides_1_1_0.mzid")
        >>> psm_list = mzid_reader.read_file()

        """
        super().__init__(filename)
        self.source = self._infer_source()
        self.searchengine_key_dict = self._get_searchengine_specific_keys()

    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""

        with mzid.read(str(self.filename.absolute())) as reader:

            for spectrum in reader:
                spectrum_title = spectrum[self.searchengine_key_dict["spectrum_key"]]
                rawfile = self._get_rawfile_name(spectrum["location"])

                for psm_dict in spectrum["SpectrumIdentificationItem"]:
                    if not self.searchengine_key_dict["score"]:
                        self.searchengine_key_dict["score"] = self._infer_score(
                            spectrum["SpectrumIdentificationItem"].keys()
                        )
                    psm = self._get_peptide_spectrum_match(
                        spectrum_title, rawfile, psm_dict
                    )
                    yield psm

    def read_file(self) -> PSMList:
        """Read full mzid file to PSM list object."""

        with mzid.read(str(self.filename.absolute())) as reader:
            # Go over each Spectrum that has a match
            psm_list = []
            # Go over each PSM
            for spectrum in reader:
                spectrum_title = spectrum[self.searchengine_key_dict["spectrum_key"]]
                rawfile = self._get_rawfile_name(spectrum["location"])

                for psm in spectrum["SpectrumIdentificationItem"]:
                    if not self.searchengine_key_dict["score"]:
                        self.searchengine_key_dict["score"] = self._infer_score(
                            spectrum["SpectrumIdentificationItem"].keys()
                        )
                    psm_list.append(
                        self._get_peptide_spectrum_match(spectrum_title, rawfile, psm)
                    )

        return PSMList(psm_list)

    def _get_peptide_spectrum_match(
        self,
        spectrum_title: str,
        rawfile: str,
        spectrum_identification_item: dict[str, Union[str, float, list]],
    ) -> PeptideSpectrumMatch:
        """Parse single mzid Record entry to a `PeptideSpectrumMatch`."""

        try:
            modifications = spectrum_identification_item["Modification"]
        except KeyError:
            modifications = []
        sequence = spectrum_identification_item["PeptideSequence"]
        peptide = self._parse_peptidoform(sequence, modifications)

        isdecoy, protein_list = self._parse_peptide_evidence_ref(
            spectrum_identification_item["PeptideEvidenceRef"]
        )
        # retention time is often missing from mzid
        try:
            rt = float(
                spectrum_identification_item[self.searchengine_key_dict["rt_key"]]
            )
        except KeyError:
            rt = float("nan")

        psm = PeptideSpectrumMatch(
            peptide=peptide,
            spectrum_id=str(spectrum_title),
            run=str(rawfile),
            is_decoy=isdecoy,
            score=float(
                spectrum_identification_item[self.searchengine_key_dict["score_key"]]
            ),
            precursor_charge=int(spectrum_identification_item["chargeState"]),
            precursor_mz=float(
                spectrum_identification_item["experimentalMassToCharge"]
            ),
            retention_time=rt,
            protein_list=protein_list,
            source=self.source,
            provenance_data=({f"mzid_filename": self.filename}),
            metadata=self._get_searchengine_specific_metadata(
                spectrum_identification_item
            ),
        )
        return psm

    @staticmethod
    def _parse_peptidoform(seq: str, modification_list: list[dict]):
        """Parse mzid sequence and modifications to Peptidoform."""

        peptide = [""] + list(seq) + [""]

        # Add modification labels
        for mod in modification_list:
            peptide[int(mod["location"])] += f"[{mod['name']}]"

        # Add dashes between residues and termnin, and join sequence
        peptide[0] = peptide[0] + "-" if peptide[0] else ""
        peptide[-1] = "-" + peptide[-1] if peptide[-1] else ""
        proforma_seq = "".join(peptide)

        return Peptidoform(proforma_seq)

    def _infer_source(self):
        """Get the source of the mzid file."""

        mzid_xml = ET.parse(self.filename)
        root = mzid_xml.getroot()
        name_space = self._get_xml_namespace(root.tag)

        return root.find(f".//{name_space}AnalysisSoftware").attrib["name"]

    @staticmethod
    def _get_xml_namespace(root_tag):
        """Get the namespace of the xml root."""

        m = re.match(r"\{.*\}", root_tag)
        return m.group(0) if m else ""

    @staticmethod
    def _parse_peptide_evidence_ref(peptide_evidence_list: list[dict]):
        """Parse peptide evidence References of PSM."""

        isdecoy = peptide_evidence_list[0]["isDecoy"]

        protein_list = [
            d["accession"] for d in peptide_evidence_list if "accession" in d.keys()
        ]

        return isdecoy, protein_list

    def _get_searchengine_specific_keys(self):
        """Get searchengine specific keys."""

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

        # if source not known return standard keys and none as score key
        return {
            "score_key": None,
            "rt_key": "retention time",
            "spectrum_key": "spectrumID",
        }

    def _get_searchengine_specific_metadata(self, SpectrumIdentificationItem):
        """Get searchengine specific psm metadata."""

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
    def _infer_score(identification_keys: list):
        """Infer the score used when source is not unknown"""

        score = None
        for score in STANDARD_SEARCHENGINE_SCORES:
            if score in identification_keys:
                score_key = score
                break
        if not score:
            raise UnknownMzidScore("No known score metric found in Mzid file")

        return score_key

    @staticmethod
    def _get_rawfile_name(file_location: str) -> str:
        """Get rawfile name out of mzid file location or filename."""

        return Path(file_location).stem


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


class UnknownMzidScore(PSMUtilsException):
    """Exception while handling or parsing peptide modifications."""

    pass
