from __future__ import annotations

import logging
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Optional, Union

from pyteomics import mzid, mass
from psims.mzid import MzIdentMLWriter
from rich.progress import Progress


from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io.exceptions import PSMUtilsIOException
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
    def __init__(self, filename: Union[str, Path], *args, **kwargs) -> None:

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
        super().__init__(filename, *args, **kwargs)
        self._source = self._infer_source()
        self._searchengine_key_dict = self._get_searchengine_specific_keys()

    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""
        with mzid.read(str(self.filename)) as reader:
            for spectrum in reader:
                spectrum_title = spectrum[self._searchengine_key_dict["spectrum_key"]]
                raw_file = Path(spectrum["location"]).stem
                for entry in spectrum["SpectrumIdentificationItem"]:
                    psm = self._get_peptide_spectrum_match(
                        spectrum_title, raw_file, entry
                    )
                    yield psm

    def read_file(self) -> PSMList:
        """Read full mzid file to PSM list object."""
        return PSMList(psm_list=[psm for psm in self.__iter__()])

    @staticmethod
    def _get_xml_namespace(root_tag):
        """Get the namespace of the xml root."""
        m = re.match(r"\{.*\}", root_tag)
        return m.group(0) if m else ""

    def _infer_source(self):
        """Get the source of the mzid file."""
        mzid_xml = ET.parse(self.filename)
        root = mzid_xml.getroot()
        name_space = self._get_xml_namespace(root.tag)
        return root.find(f".//{name_space}AnalysisSoftware").attrib["name"]

    def _get_searchengine_specific_keys(self):
        """Get searchengine specific keys."""
        if "PEAKS" in self._source:
            return {
                "score_key": "PEAKS:peptideScore",
                "rt_key": "retention time",
                "spectrum_key": "spectrumID",
            }
        elif self._source == "MS-GF+":
            return {
                "score_key": "MS-GF:RawScore",
                "rt_key": "scan start time",
                "spectrum_key": "spectrum title",
            }
        # if source is not known return standard keys and infer score key
        return {
            "score_key": self._infer_score_name(),
            "rt_key": "retention time",
            "spectrum_key": "spectrumID",
        }

    def _infer_score_name(self) -> str:
        """Infer the score from the known list of PSM scores."""
        with mzid.read(str(self.filename)) as reader:
            for spectrum in reader:
                sii_keys = spectrum["SpectrumIdentificationItem"][0].keys()
                break
        for score in STANDARD_SEARCHENGINE_SCORES:
            if score in sii_keys:
                return score
        else:
            raise UnknownMzidScore("No known score metric found in mzIdentML file.")

    @staticmethod
    def _parse_peptidoform(seq: str, modification_list: list[dict], charge: int):
        """Parse mzid sequence and modifications to Peptidoform."""
        peptide = [""] + list(seq) + [""]

        # Add modification labels
        for mod in modification_list:
            peptide[int(mod["location"])] += f"[{mod['name']}]"

        # Add dashes between residues and termini, and join sequence
        peptide[0] = peptide[0] + "-" if peptide[0] else ""
        peptide[-1] = "-" + peptide[-1] if peptide[-1] else ""
        proforma_seq = "".join(peptide)

        # Add charge state
        proforma_seq += f"/{charge}"

        return Peptidoform(proforma_seq)

    @staticmethod
    def _parse_peptide_evidence_ref(peptide_evidence_list: list[dict]):
        """Parse peptide evidence References of PSM."""
        isdecoy = peptide_evidence_list[0]["isDecoy"]
        protein_list = [
            d["accession"] for d in peptide_evidence_list if "accession" in d.keys()
        ]
        return isdecoy, protein_list

    def _get_searchengine_specific_metadata(self, spectrum_identification_item):
        """Get searchengine specific psm metadata."""
        sii = spectrum_identification_item
        metadata = {
            "calculatedMassToCharge": sii["calculatedMassToCharge"],
            "rank": sii["rank"],
            "length": len(sii["PeptideSequence"]),
        }
        if self._source == "MS-GF+":
            metadata.update(
                {
                    "MS-GF:DeNovoScore": sii["MS-GF:DeNovoScore"],
                    "MS-GF:EValue": sii["MS-GF:EValue"],
                    "MS-GF:DeNovoScore": sii["MS-GF:DeNovoScore"],
                    "IsotopeError": sii["IsotopeError"],
                }
            )
        return metadata

    def _get_peptide_spectrum_match(
        self,
        spectrum_title: str,
        raw_file: str,
        spectrum_identification_item: dict[str, Union[str, float, list]],
    ) -> PeptideSpectrumMatch:
        """Parse single mzid entry to :py:class:`~psm_utils.peptidoform.Peptidoform`."""
        sii = spectrum_identification_item

        try:
            modifications = sii["Modification"]
        except KeyError:
            modifications = []
        sequence = sii["PeptideSequence"]
        peptide = self._parse_peptidoform(sequence, modifications, sii["chargeState"])
        is_decoy, protein_list = self._parse_peptide_evidence_ref(
            sii["PeptideEvidenceRef"]
        )
        try:
            rt = float(sii[self._searchengine_key_dict["rt_key"]])
        except KeyError:
            rt = float("nan")

        psm = PeptideSpectrumMatch(
            peptide=peptide,
            spectrum_id=spectrum_title,
            run=raw_file,
            is_decoy=is_decoy,
            score=sii[self._searchengine_key_dict["score_key"]],
            precursor_mz=sii["experimentalMassToCharge"],
            retention_time=rt,
            protein_list=protein_list,
            source=self._source,
            provenance_data={"mzid_filename": str(self.filename)},
            metadata=self._get_searchengine_specific_metadata(sii),
        )
        return psm


class MzidWriter(WriterBase):
    """Writer for MzIdentMl files."""

    def __init__(
        self,
        filename: Union[str, Path],
        example_psm: Optional[PeptideSpectrumMatch] = None,
        *args,
        **kwargs,
    ):
        """
        Writer for psm_utils MZID format.

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
        super().__init__(filename, *args, **kwargs)

        self._writer = None

    def __enter__(self) -> MzidWriter:
        pass

    def __exit__(self, *args, **kwargs) -> None:
        pass

    def write_psm(self, psm: PeptideSpectrumMatch):
        """Write a single PSM to a Peptidehit object."""
        raise NotImplemented("Mzid writer does not support write psm")

    def write_file(self, psm_list: PSMList, disable_progressbar: bool = True):
        """Write entire PSMList to mzid file."""
        file = open(self.filename, "wb")
        writer = MzIdentMLWriter(file, close=True)
        with Progress(disable=disable_progressbar) as progress:

            with writer:
                writer.controlled_vocabularies()
                writer.provenance(
                    software={
                        "name": "psm_utils",
                        "uri": "https://github.com/compomics/psm_utils",
                        "version": "0.0.1",  # function to add automatically
                    }
                )
                writer.register("SpectraData", 1)
                writer.register("SearchDatabase", 1)
                writer.register("SpectrumIdentificationList", 1)
                writer.register("SpectrumIdentificationProtocol", 1)

                proteins = set()
                peptide_ids = set()
                peptide_evidence_ids = set()

                proteins = {
                    prot
                    for prot_list in list(psm_list["protein_list"])
                    for prot in prot_list
                }

                spec_id_dict = psm_list.get_psm_dict()
                task1 = progress.add_task(
                    "[cyan]Writing Proteins to mzid", total=len(proteins)
                )
                task2 = progress.add_task(
                    "[cyan]Writing Peptides/PeptideEvidences to mzid",
                    total=len(psm_list),
                )
                task3 = progress.add_task(
                    "[cyan]Writing SpectrumIdentificationResults to mzid",
                    total=len(psm_list),
                )

                with writer.sequence_collection():

                    for prot in proteins:
                        writer.write_db_sequence(prot, None, id=prot, params=[])
                        progress.update(task1, advance=1)
                    for psm in psm_list:
                        peptide = psm["peptide"]
                        if peptide not in peptide_ids:
                            writer.write_peptide(**self._create_peptide_object(peptide))
                            peptide_ids.add(peptide)

                        for protein in psm["protein_list"]:
                            peptide_evidence_id = (
                                f"PeptideEvidence_{peptide.proforma}_{protein}"
                            )
                            if peptide_evidence_id not in peptide_evidence_ids:
                                peptide_evidence_ids.add(peptide_evidence_id)
                                writer.write_peptide_evidence(
                                    peptide_id="Peptide_" + peptide.proforma,
                                    db_sequence_id=protein,
                                    id=peptide_evidence_id,
                                    start_position=None,
                                    end_position=None,
                                    is_decoy=psm["is_decoy"],
                                )
                        progress.update(task2, advance=1)
                with writer.analysis_collection():
                    writer.SpectrumIdentification([1], [1]).write(writer)

                with writer.analysis_protocol_collection():
                    writer.spectrum_identification_protocol()  # build without?

                with writer.data_collection():
                    spectra_data, spectra_data_id_dict = self._transform_spectra_data(
                        spec_id_dict=spec_id_dict
                    )
                    writer.inputs(
                        source_files=[],
                        # search_databases=transform_search_database(), # if fasta file is given we can parse here and add protein information
                        spectra_data=spectra_data,
                    )
                with writer.analysis_data():
                    with writer.spectrum_identification_list(id=1):
                        for collection in spec_id_dict.keys():
                            for run in spec_id_dict[collection].keys():
                                spectra_data_id = spectra_data_id_dict[
                                    "/".join(filter(None, [collection, run]))
                                ]
                                for spec_id in spec_id_dict[collection][run].keys():
                                    identified_psms = spec_id_dict[collection][run][
                                        spec_id
                                    ]
                                    writer.write_spectrum_identification_result(
                                        **self._transform_spectrum_identification_result(
                                            spec_id, identified_psms, spectra_data_id
                                        )
                                    )
                                    progress.update(
                                        task3,
                                        advance=len(
                                            spec_id_dict[collection][run][spec_id]
                                        ),
                                    )
        try:
            writer.close()
        except OSError:
            pass

    @staticmethod
    def _create_peptide_object(peptidoform):
        """Create mzid peptide object from petidoform"""

        parsed_sequence = peptidoform.parsed_sequence
        properties = peptidoform.properties
        peptide = peptidoform.sequence
        modifications = []

        # makes no effort to handle terminal modifications or snps at this time
        for loc, (aa, mod) in enumerate(parsed_sequence, start=1):
            if mod:
                mzid_mod = {
                    "location": loc,
                    "name": mod[0].name,  # what with multiple mods on one aa
                    "monoisotopic_mass_delta": mod[0].mass,
                }
                modifications.append(mzid_mod)

        if properties["n_term"]:
            modifications.append(
                {
                    "location": 0,
                    "name": mod[0].name,  # what with multiple mods on one aa
                    "monoisotopic_mass_delta": mod[0].mass,
                }
            )

        if properties["c_term"]:
            modifications.append(
                {
                    "location": len(peptide) + 1,
                    "name": mod[0].name,  # what with multiple mods on one aa
                    "monoisotopic_mass_delta": mod[0].mass,
                }
            )

        peptide_object = {
            "id": "Peptide_" + peptidoform.proforma,
            "modifications": modifications,
            "peptide_sequence": peptide,
        }
        return peptide_object

    def _transform_search_database(self):
        """Create mzid databse object"""

        # TODO: Create this and link with protein object when fasta file is provided
        return {
            "file_format": "fasta format",
            "name": "",
            "id": 1,
            "location": "",
            "params": [],
        }

    @staticmethod
    def _transform_spectra_data(spec_id_dict: dict):
        """Get all the unique spectra data from psm list spectrum id dict"""

        collection_run_id_dict = {}
        spectra_data = []
        i = 1
        for collection in spec_id_dict.keys():
            for run in spec_id_dict[collection].keys():
                collection_run_id = "/".join(filter(None, [collection, run]))

                if collection_run_id not in collection_run_id_dict.keys():
                    collection_run_id_dict[collection_run_id] = i

                    spectra_data_object = {
                        "id": i,
                        "location": collection_run_id,
                        "spectrum_id_format": "multiple peak list nativeID format",
                        # 'file_format': #TODO can we infer this?
                    }

                spectra_data.append(spectra_data_object)
        return spectra_data, collection_run_id_dict

    @staticmethod
    def _transform_spectrum_identification_item(candidate_psm, i=0):
        """Create mzid spectrum identification item for each candidate psm"""

        peptide = candidate_psm["peptide"].proforma
        items = []
        for prot_acc in candidate_psm["protein_list"]:
            theo_mass = candidate_psm["peptide"].theoretical_mass
            charge = candidate_psm["peptide"].precursor_charge
            theo_mz = (theo_mass + (mass.nist_mass["H"][1][0] * charge)) / charge
            try:
                exper_mz = candidate_psm["precursor_mz"]
                # mass_error = abs(theo_mz - exper_mz)
            except:
                exper_mz = None
                # mass_error = None

            candidate_psm_dict = {
                "charge_state": charge,
                "peptide_id": f"Peptide_" + peptide,
                "peptide_evidence_id": f"PeptideEvidence_{peptide}_{prot_acc}",
                "id": f"SII_{candidate_psm['spectrum_id']}_{peptide}_{prot_acc}",
                "score": {
                    "score": candidate_psm["score"]
                },  # add which search engine score cv param
                "experimental_mass_to_charge": exper_mz,
                "calculated_mass_to_charge": theo_mz,
                "rank": candidate_psm["rank"],
                "params": [
                    {x: candidate_psm["metadata"][x]}
                    for x in candidate_psm["metadata"].keys()
                ],
            }
            items.append(candidate_psm_dict)
        return items

    def _transform_spectrum_identification_result(
        self, spec_id, identified_psms, spectra_data_id
    ):
        """Create mzid spectrum identification result object for each psm that match the same spectrum"""

        spectrum_id_result = {
            "id": f"SIR_{spec_id}",
            "spectrum_id": spec_id,
            "spectra_data_id": spectra_data_id,
            "params": [
                {
                    "name": "scan start time",
                    "value": identified_psms[0]["retention_time"],
                }
            ],  # add function to determine the unit
        }
        identifications = []
        for i, candidate in enumerate(identified_psms):
            identifications.extend(
                self._transform_spectrum_identification_item(candidate, i)
            )
        spectrum_id_result["identifications"] = identifications
        return spectrum_id_result


class UnknownMzidScore(PSMUtilsIOException):
    """Exception while handling or parsing peptide modifications."""

    pass
