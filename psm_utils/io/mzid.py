"""
Reader and writers for the HUPO-PSI mzIdentML format.

See `psidev.info/mzidentml <https://psidev.info/mzidentml>`_ for more info on the
format.

"""

from __future__ import annotations

import logging
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Optional, Union

from psims.mzid import MzIdentMLWriter
from pyteomics import mzid, proforma
from rich.progress import Progress

from psm_utils import __version__
from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

logger = logging.getLogger(__name__)

# Excerpt from MS:1001143 items (PSM-level search engine specific statistic)
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
    "Percolator:score",
    "xi:score",
]


class MzidReader(ReaderBase):
    def __init__(self, filename: Union[str, Path], *args, **kwargs) -> None:

        """
        Reader for mzIdentML PSM files.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to PSM file.

        Examples
        --------

        MzidReader supports iteration:

        >>> from psm_utils.io.mzid import MzidReader
        >>> for psm in MzidReader("peptides_1_1_0.mzid"):
        ...     print(psm.peptidoform.proforma)
        ACDEK
        AC[Carbamidomethyl]DEFGR
        [Acetyl]-AC[Carbamidomethyl]DEFGHIK

        Or a full file can be read at once into a :py:class:`psm_utils.psm_list.PSMList`
        object:

        >>> mzid_reader = MzidReader("peptides_1_1_0.mzid")
        >>> psm_list = mzid_reader.read_file()

        """
        super().__init__(filename, *args, **kwargs)
        self._non_metadata_keys = None
        self._score_key = None
        self._rt_key = None

        self._source = self._infer_source()

    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""
        with mzid.read(str(self.filename)) as reader:
            for spectrum in reader:
                spectrum_id = spectrum["spectrumID"]
                raw_file = Path(spectrum["location"]).stem
                for entry in spectrum["SpectrumIdentificationItem"]:
                    if not self._non_metadata_keys:
                        self._get_non_metadata_keys(entry.keys())
                    try:
                        psm = self._get_peptide_spectrum_match(
                            spectrum_id, raw_file, entry, spectrum["spectrum title"]
                        )
                    except KeyError:
                        psm = self._get_peptide_spectrum_match(
                            spectrum_id, raw_file, entry
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

    @staticmethod
    def _infer_score_name(keys) -> str:
        """Infer the score from the known list of PSM scores."""

        for score in STANDARD_SEARCHENGINE_SCORES:
            if score in keys:
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

    def _get_peptide_spectrum_match(
        self,
        spectrum_id: str,
        raw_file: str,
        spectrum_identification_item: dict[str, Union[str, float, list]],
        spectrum_title: Optional[str] = None,
    ) -> PSM:
        """Parse single mzid entry to :py:class:`~psm_utils.peptidoform.Peptidoform`."""
        sii = spectrum_identification_item

        try:
            modifications = sii["Modification"]
        except KeyError:
            modifications = []
        sequence = sii["PeptideSequence"]
        peptidoform = self._parse_peptidoform(
            sequence, modifications, sii["chargeState"]
        )
        is_decoy, protein_list = self._parse_peptide_evidence_ref(
            sii["PeptideEvidenceRef"]
        )
        try:
            precursor_mz = sii["experimentalMassToCharge"]
        except KeyError:
            precursor_mz = None

        if self._rt_key:
            rt = sii[self._rt_key]
        else:
            rt = float("nan")

        metadata = {
            col: str(sii[col])
            for col in sii.keys()
            if col not in self._non_metadata_keys
        }
        if spectrum_title:
            metadata["spectrum title"] = spectrum_title

        psm = PSM(
            peptidoform=peptidoform,
            spectrum_id=spectrum_id,
            run=raw_file,
            is_decoy=is_decoy,
            score=sii[self._score_key],
            precursor_mz=precursor_mz,
            retention_time=rt,
            protein_list=protein_list,
            rank=sii["rank"],
            source=self._source,
            provenance_data={"mzid_filename": str(self.filename)},
            metadata=metadata,
        )
        return psm

    def _get_non_metadata_keys(self, keys: list):
        """Gather all the keys that should not be written to metadata"""

        # All keys required to create PSM object
        default_keys = [
            "chargeState",
            "rank",
            "PeptideSequence",
            "experimentalMassToCharge",
            "PeptideEvidenceRef",
            "Modification",
        ]
        # Get the score key and add to default keys
        self._score_key = self._infer_score_name(keys)
        default_keys.append(self._score_key)

        # Get retention time key
        for rt_key in ["retention time", "scan start time"]:
            if rt_key in keys:
                self._rt_key = rt_key
                default_keys.append(rt_key)
                break

        # Keys that are not necessary for metadata
        self._non_metadata_keys = ["ContactRole", "passThreshold"]
        self._non_metadata_keys.extend(default_keys)


class MzidWriter(WriterBase):
    """Writer for MzIdentMl files."""

    def __init__(
        self,
        filename: Union[str, Path],
        show_progressbar: bool = False,
        *args,
        **kwargs,
    ):
        """
        Writer for mzIdentML PSM files.

        Parameters
        ----------
        filename: str, Pathlib.Path
            Path to PSM file.
        show_progressbar: bool, optional
            Show progress bar for conversion process. (default: False)

        Notes
        -----
        Unlike other psm_utils.io writer classes, :py:class:`MzidWriter` does not
        support writing a single PSM to a file with the :py:meth:`write_psm` method.
        Only writing a full PSMList to a file at once with the :py:meth:`write_file`
        method is currently supported.

        """
        super().__init__(filename, *args, **kwargs)
        self.show_progressbar = show_progressbar
        self._writer = None

    def __enter__(self) -> MzidWriter:
        return self

    def __exit__(self, *args, **kwargs) -> None:
        pass

    def write_psm(self, psm: PSM):
        """
        Write a single PSM to the PSM file.

        This method is currently not supported (see Notes).

        Raises
        ------
        NotImplementedError
            MzidWriter currently does not support write_psm.
        """
        raise NotImplementedError("MzidWriter currently does not support write_psm.")

    def write_file(self, psm_list: PSMList):
        """Write entire PSMList to mzid file."""
        file = open(self.filename, "wb")
        with Progress(disable=(not self.show_progressbar)) as progress:
            with MzIdentMLWriter(file, close=True) as writer:
                writer.controlled_vocabularies()
                writer.provenance(
                    software={
                        "name": "psm_utils",
                        "uri": "https://github.com/compomics/psm_utils",
                        "version": __version__,
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
                    "[cyan]Writing Peptide and PeptideEvidence items",
                    total=len(psm_list),
                )
                task3 = progress.add_task(
                    "[cyan]Writing SpectrumIdentificationResults",
                    total=len(psm_list),
                )

                with writer.sequence_collection():
                    for prot in proteins:
                        writer.write_db_sequence(prot, None, id=prot, params=[])
                        progress.update(task1, advance=1)
                    for psm in psm_list:
                        peptide = psm["peptidoform"]
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

    @staticmethod
    def _create_peptide_object(peptidoform):
        """Create mzid peptide object from Peptidoform."""

        def parse_modifications(modifications: list[proforma.TagBase], location: int):
            modification_list = []
            if modifications:
                for mod in modifications:
                    try:
                        modification_list.append(
                            {
                                "location": location,
                                "name": mod.name,
                                "monoisotopic_mass_delta": mod.mass,
                            }
                        )
                    except AttributeError:
                        modification_list.append(
                            {
                                "location": location,
                                "monoisotopic_mass_delta": mod.mass,
                            }
                        )
            return modification_list

        # Parse modifications
        modifications = []
        for loc, (aa, mods) in enumerate(peptidoform.parsed_sequence, start=1):
            modifications.extend(parse_modifications(mods, loc))
        modifications.extend(parse_modifications(peptidoform.properties["n_term"], 0))
        modifications.extend(
            parse_modifications(
                peptidoform.properties["c_term"], len(peptidoform.sequence) + 1
            )
        )

        peptide_object = {
            "id": "Peptide_" + peptidoform.proforma,
            "peptide_sequence": peptidoform.sequence,
            "modifications": modifications,
        }

        return peptide_object

    def _transform_search_database(self):
        """Create mzid database object."""
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
        """Get all the unique spectra data from PSMList spectrum id dict."""
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
        """Create SpectrumIdentificationItem for each candidate PSM."""
        peptide = candidate_psm["peptidoform"].proforma
        if candidate_psm["metadata"]:
            params = [{k: v} for k, v in candidate_psm["metadata"].items()]
        else:
            params = []
        candidate_psm_dict = {
            "charge_state": candidate_psm["peptidoform"].precursor_charge,
            "peptide_id": f"Peptide_" + peptide,
            # TODO: add which search engine score cv param
            "score": {"score": candidate_psm["score"]},
            "experimental_mass_to_charge": candidate_psm["precursor_mz"],
            "calculated_mass_to_charge": candidate_psm["peptidoform"].theoretical_mz,
            "rank": candidate_psm["rank"],
            "params": params,
        }
        items = []
        for prot_acc in candidate_psm["protein_list"]:
            protein_specific_items = {
                "peptide_evidence_id": f"PeptideEvidence_{peptide}_{prot_acc}",
                "id": f"SII_{candidate_psm['spectrum_id']}_{peptide}_{prot_acc}",
            }
            items.append(dict(candidate_psm_dict, **protein_specific_items))
        return items

    def _transform_spectrum_identification_result(
        self, spec_id, identified_psms, spectra_data_id
    ):
        """Create mzid SpectrumIdentificationResult object for PSMs that match the same spectrum."""
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
