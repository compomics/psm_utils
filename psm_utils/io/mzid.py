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
from typing import Union

from psims.mzid import MzIdentMLWriter
from pyteomics import mzid, proforma
from rich.progress import Progress

from psm_utils import __version__
from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

from lxml import etree
import copy

logger = logging.getLogger(__name__)

# Excerpt from MS:1001143 items (PSM-level search engine specific statistic)
# Not all child terms are used, as not all statistics are direct scores.
# Items are sorted by priority (if more scores are present, the first found one is used)
STANDARD_SEARCHENGINE_SCORES = [
    "peptideshaker psm score",
    "amanda:amandascore",
    "andromeda:score",
    "byonic:score",
    "comet:xcorr",
    "debunker:score",
    "identitye score",
    "ksdp score",
    "ms-gf:specevalue",
    "ms-gf:evalue",
    "ms-gf:rawscore",
    "ms-gf:denovoscore",
    "msfit:mowse score",
    "mspathfinder:rawscore",
    "mspepsearch:score",
    "mascot:score",
    "metamorpheus:score",
    "omssa:evalue",
    "openpepxl:score",
    "peaks:peptidescore",
    "phenyx:pepzscore",
    "prolucid:xcorr",
    "prosight:specral c-score",
    "profound:z value",
    "proteinprospector:score",
    "proteinscape:sequestmetascore",
    "proteomediscoverer:delta score",
    "proteome discoverer delta score",
    "sequest:xcorr",
    "sim-xl score ",
    "sqid:score ",
    "sonar:score",
    "spectrummill:score",
    "topmg:spectral e-value",
    "x!tandem:hyperscore",
    "zcore:probscore:",
    "percolator:score",
    "xi:score",
    "search engine specific score",
]
Q_VALUE_TERMS = [
    "PSM-level q-value",
    "MS-GF:QValue",
    "MSPathFinder:QValue",
    "ProSight:spectral Q-value",
]
PEP_TERMS = [
    "PSM-level local FDR",
    "Andromeda:PEP",
    "Byonic:PEP",
    "MaxQuant-DIA PEP",
    "MS-GF:PEP",
    "OpenMS:ConsensusID PEP",
    "percolator:PEP",
]


class MzidReader(ReaderBase):
    def __init__(self, filename: str | Path, *args, score_key: str = None, **kwargs) -> None:
        """
        Reader for mzIdentML PSM files.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to PSM file.
        score_key: str, optional
            Name of the score metric to use as PSM score. If not provided, the score metric is
            inferred from the file if one of the child parameters of ``MS:1001143`` is present.

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

        Notes
        -----
        - :py:class:`MzidReader` looks for the ``retention time`` or ``scan start time`` cvParams
          in both SpectrumIdentificationResult and SpectrumIdentificationItem levels. Note that
          according to the mzIdentML specification document (v1.1.1) neither cvParams are expected
          to be present at either levels.
        - For the :py:attr:`PSM.spectrum_id` property, the ``spectrum title`` cvParam is preferred
          over the ``spectrumID`` attribute, as these titles always match the titles in the  peak
          list files. ``spectrumID`` is then saved in ``PSM.metadata["mzid_spectrum_id"]``.
          If ``spectrum title`` is absent, ``spectrumID`` is saved to :py:attr:`PSM.spectrum_id`.

        """
        super().__init__(filename, *args, **kwargs)
        self._non_metadata_keys = ["ContactRole", "passThreshold"]
        self._score_key = score_key
        self._rt_key = None
        self._spectrum_rt_key = None
        self._qvalue_key = None
        self._pep_key = None
        self._im_key = None

        self._source = self._infer_source()

    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""
        with mzid.read(str(self.filename)) as reader:
            # Parse spectrum metadata
            self._get_toplevel_non_metadata_keys(reader[0].keys())
            # Parse PSM non-metadata keys, rt key and score key
            self._get_non_metadata_keys(reader[0]["SpectrumIdentificationItem"][0].keys())

            for spectrum in reader:
                # Parse spectrum metadata
                spectrum_id = spectrum["spectrumID"]
                spectrum_title = (
                    spectrum["spectrum title"] if "spectrum title" in spectrum else None
                )
                run = Path(spectrum["location"]).stem if "location" in spectrum else None
                rt = float(spectrum[self._spectrum_rt_key]) if self._spectrum_rt_key else None
                ion_mobility = float(spectrum[self._im_key]) if self._im_key else None

                # Parse PSMs from spectrum
                for entry in spectrum["SpectrumIdentificationItem"]:
                    yield self._get_peptide_spectrum_match(
                        spectrum_id, spectrum_title, run, rt, ion_mobility, entry
                    )

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
        try:
            return root.find(f".//{name_space}AnalysisSoftware").attrib["name"]
        except KeyError:
            return None

    @staticmethod
    def _parse_peptidoform(seq: str, modification_list: list[dict], charge: Union[int, None]):
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
        if charge:
            proforma_seq += f"/{charge}"

        return Peptidoform(proforma_seq)

    @staticmethod
    def _parse_peptide_evidence_ref(peptide_evidence_list: list[dict]):
        """
        Parse PeptideEvidence list of PSM.

        Notes
        -----
        If multiple PeptideEvidence entries are associated with the PSM, the PSM is only considered
        a decoy entry if ALL PeptideEvidence entries are decoy entries. If a target PeptideEvidence
        entry is present, it should get priority over decoy entries. In theory, no overlap between
        target and decoy peptide sequence should be present in the search space, although this
        might not have been filtered for by the search engine.

        """
        isdecoy = all(
            [entry["isDecoy"] if "isDecoy" in entry else None for entry in peptide_evidence_list]
        )
        protein_list = [d["accession"] for d in peptide_evidence_list if "accession" in d.keys()]
        return isdecoy, protein_list

    def _get_peptide_spectrum_match(
        self,
        spectrum_id: str,
        spectrum_title: Union[str, None],
        run: Union[str, None],
        rt: Union[float, None],
        ion_mobility: Union[float, None],
        spectrum_identification_item: dict[str, str | float | list],
    ) -> PSM:
        """Parse single mzid entry to :py:class:`~psm_utils.peptidoform.Peptidoform`."""
        sii = spectrum_identification_item
        try:
            modifications = sii["Modification"]
        except KeyError:
            modifications = []
        sequence = sii["PeptideSequence"]
        charge = sii["chargeState"] if "chargeState" in sii else None
        peptidoform = self._parse_peptidoform(sequence, modifications, charge)
        is_decoy, protein_list = self._parse_peptide_evidence_ref(sii["PeptideEvidenceRef"])
        try:
            precursor_mz = sii["experimentalMassToCharge"]
        except KeyError:
            precursor_mz = None

        # Override spectrum-level RT if present at PSM level
        if self._rt_key:
            rt = float(sii[self._rt_key])

        metadata = {col: str(sii[col]) for col in sii.keys() if col not in self._non_metadata_keys}

        # Prefer spectrum title (identical to titles in MGF), fall back to spectrumID
        if spectrum_title:
            metadata["mzid_spectrum_id"] = spectrum_id
            psm_spectrum_id = spectrum_title
        else:
            psm_spectrum_id = spectrum_id

        try:
            score = sii[self._score_key]
        except KeyError:
            score = None
        psm = PSM(
            peptidoform=peptidoform,
            spectrum_id=psm_spectrum_id,
            run=run,
            is_decoy=is_decoy,
            score=score,
            qvalue=sii[self._qvalue_key] if self._qvalue_key else None,
            pep=sii[self._pep_key] if self._pep_key else None,
            precursor_mz=precursor_mz,
            retention_time=rt,
            ion_mobility=ion_mobility,
            protein_list=protein_list,
            rank=sii["rank"] if "rank" in sii else None,
            source=self._source,
            provenance_data={"mzid_filename": str(self.filename)},
            metadata=metadata,
        )
        return psm

    def _get_non_metadata_keys(self, keys: list):
        """Gather all the keys at PSM-level that should not be written to metadata."""
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
        if not self._score_key:
            self._score_key = self._infer_score_name(keys)
        if self._score_key:
            default_keys.append(self._score_key)
        else:
            logger.warning(
                "No known score metric found in mzIdentML file. Scores will be set to None."
            )

        # Get the q-value key and add to default keys
        self._qvalue_key = self._infer_qvalue_name(keys)
        if self._qvalue_key:
            default_keys.append(self._qvalue_key)

        # Get the PEP key and add to default keys
        self._pep_key = self._infer_pep_name(keys)
        if self._pep_key:
            default_keys.append(self._pep_key)

        # Get retention time key
        for rt_key in ["retention time", "scan start time"]:
            if rt_key in keys:
                self._rt_key = rt_key
                default_keys.append(rt_key)
                break

        # Keys that are not necessary for metadata
        self._non_metadata_keys.extend(default_keys)

    def _get_toplevel_non_metadata_keys(self, keys: list):
        """Gather all keys at spectrum-level that should not be written to metadata."""
        # Check if RT is encoded in spectrum metadata
        for key in ["retention time", "scan start time"]:
            if key in keys:
                self._spectrum_rt_key = key
                self._non_metadata_keys.append(key)
                break

        # Check if ion mobility is encoded in spectrum metadata
        for im_key in ["inverse reduced ion mobility"]:
            if im_key in keys:
                self._im_key = im_key
                self._non_metadata_keys.append(im_key)
                break

    @staticmethod
    def _infer_score_name(keys) -> str:
        """Infer the score from the list of known PSM scores."""
        lower_keys = {key.lower(): key for key in keys}
        for score in STANDARD_SEARCHENGINE_SCORES:
            if score in lower_keys:
                return lower_keys[score]

    @staticmethod
    def _infer_qvalue_name(keys) -> Union[str, None]:
        """Infer the q-value term from the list of known terms."""
        for qvalue in Q_VALUE_TERMS:
            if qvalue in keys:
                return qvalue
        else:
            return None

    @staticmethod
    def _infer_pep_name(keys) -> Union[str, None]:
        """Infer the PEP term from the list of known terms."""
        for pep in PEP_TERMS:
            if pep in keys:
                return pep
        else:
            return None


class MzidQuickReader(ReaderBase):
    def __init__(self, filename: str | Path, *args, score_key: str = None, **kwargs) -> None:
        """
        Quick, and not totally complete reader for mzIdentML PSM files.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to PSM file.
        score_key: str, optional
            Name of the score metric to use as PSM score. If not provided, the score metric is
            inferred from the file if one of the child parameters of ``MS:1001143`` is present.

        Examples
        --------

        MzidQuickReader supports iteration like MzidReader.

        Notes
        -----
        - :py:class:`MzidReader` looks for the ``retention time`` or ``scan start time`` cvParams
          in both SpectrumIdentificationResult and SpectrumIdentificationItem levels. Note that
          according to the mzIdentML specification document (v1.1.1) neither cvParams are expected
          to be present at either levels.
        - For the :py:attr:`PSM.spectrum_id` property, the ``spectrum title`` cvParam is preferred
          over the ``spectrumID`` attribute, as these titles always match the titles in the  peak
          list files. ``spectrumID`` is then saved in ``PSM.metadata["mzid_spectrum_id"]``.
          If ``spectrum title`` is absent, ``spectrumID`` is saved to :py:attr:`PSM.spectrum_id`.

        """
        super().__init__(filename, *args, **kwargs)
        self._non_metadata_keys = ["ContactRole", "passThreshold"]
        self._score_key = score_key
        self._rt_key = None
        self._spectrum_rt_key = None
        self._qvalue_key = None
        self._pep_key = None
        self._im_key = None

        # some helper-dictionaries
        self.peptides_dict = {}
        self.peptide_evidences_dict = {}
        self.db_sequences_dict = {}
        self.search_dbs_dict = {}
        self.spectra_data_dict = {}

        self._preparse_references()

    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""
        first_sir = True

        for event, element in etree.iterparse(str(self.filename), events=("end", ), tag=("{*}SpectrumIdentificationResult")):
            spectrum = self._parse_sir(element)

            if first_sir:
                # Parse spectrum metadata
                self._get_toplevel_non_metadata_keys(spectrum.keys())

                # Parse PSM non-metadata keys, rt key and score key
                self._get_non_metadata_keys(spectrum["SpectrumIdentificationItem"][0].keys())
                first_sir = False
            
            spectrum_id = spectrum["spectrumID"]
            spectrum_title = (
                spectrum["spectrum title"] if "spectrum title" in spectrum else None
            )
            run = Path(spectrum["location"]).stem if "location" in spectrum else None
            rt = float(spectrum[self._spectrum_rt_key]) if self._spectrum_rt_key else None
            ion_mobility = float(spectrum[self._im_key]) if self._im_key else None
            
            # Parse PSMs from spectrum
            for entry in spectrum["SpectrumIdentificationItem"]:
                yield self._get_peptide_spectrum_match(
                    spectrum_id, spectrum_title, run, rt, ion_mobility, entry
                )

    def _preparse_references(self) -> None:
        """pre-parses all information relevant for references"""
        for _, element in etree.iterparse(str(self.filename), events=("end", ), tag=("{*}Peptide", "{*}PeptideEvidence", "{*}DBSequence", "{*}SearchDatabase", "{*}SpectraData", "{*}AnalysisSoftware")):
            tag = element.tag.rpartition("}")[2]

            if (tag == "Peptide"):
                self.peptides_dict |= MzidQuickReader._parse_peptide(element)
            elif (tag == "PeptideEvidence"):
                self.peptide_evidences_dict |= MzidQuickReader._parse_peptideevidence(element)
            elif (tag == "DBSequence"):
                self.db_sequences_dict |= MzidQuickReader._parse_dbsequence(element)
            elif (tag == "SearchDatabase"):
                self.search_dbs_dict |= MzidQuickReader._parse_searchdb(element)
            elif (tag == "SpectraData"):
                self.spectra_data_dict |= MzidQuickReader._parse_spectradata(element)
            elif (tag == "AnalysisSoftware"):
                self._source = MzidQuickReader._parse_analysissoftware(element)
            
            element.clear()

    @staticmethod
    def _parse_peptide(peptide_element: etree.Element) -> dict[str, dict]:
        pep_id = None
        sequence = None
        modifications = []

        # parse the Peptide's attributes
        for idx, item in peptide_element.items():
            if (idx == "id"):
                pep_id = item
        
        for event, item in etree.iterwalk(peptide_element, events=("start", "end"), tag=("{*}PeptideSequence", "{*}Modification")):
            if event == "start":
                # strip the namespace
                tag = item.tag.rpartition("}")[2]

                if (tag == "PeptideSequence"):
                    sequence = item.text
                elif (tag == "Modification"):
                    modifications.append(MzidQuickReader._parse_modification(item))
            else:
                item.clear()
        
        return {pep_id: {"PeptideSequence": sequence, "Modification" : modifications}}
    
    @staticmethod
    def _parse_modification(modification: etree.Element) -> dict[str, str]:
        # parse the Modification's attributes
        params = MzidQuickReader._parse_elements_attributes(modification)

        for event, item in etree.iterwalk(modification, events=("start", "end")):
            if event == "start":
                tag = item.tag.rpartition("}")[2]

                if (tag == "cvParam"):
                    params_attributes = MzidQuickReader._parse_elements_attributes(item)
                    params["name"] = params_attributes["name"]
            else:
                item.clear()
        
        return params
    
    @staticmethod
    def _parse_elements_attributes(param: etree.Element):
        attributes = {}
        for idx, item in param.items():
            attributes[idx] = item
        
        return attributes

    @staticmethod
    def _parse_peptideevidence(pepevidence_element: etree.Element) -> dict[str, dict]:
        # parse the PeptideEvidence's attributes
        attributes = MzidQuickReader._parse_elements_attributes(pepevidence_element)
        pep_ev_id = attributes["id"]
        del attributes["id"]

        # transform some types
        if "end" in attributes.keys():
            attributes["end"] = int(attributes["end"])
        if "isDecoy" in attributes.keys():
            attributes["isDecoy"] = MzidQuickReader._text_to_boolean(attributes["isDecoy"])
        if "start" in attributes.keys():
            attributes["start"] = int(attributes["start"])
        
        # there could be cvParams or userParams, but they don't have any dedicated meaninfg (yet)
        return {pep_ev_id : attributes}

    @staticmethod
    def _text_to_boolean(text: str) -> bool:
        return text.lower() in ("yes", "true", "t", "1")

    @staticmethod
    def _parse_dbsequence(dbseq_element: etree.Element) -> dict[str, dict]:
        dbseq_id = None

        # parse the DBSequences's attributes
        attributes = MzidQuickReader._parse_elements_attributes(dbseq_element)
        dbseq_id = attributes["id"]
        del attributes["id"]

        # transform some types
        if "length" in attributes.keys():
            attributes["length"] = int(attributes["length"])

        # get cvParams and userParams (mapping: name -> value)
        for event, item in etree.iterwalk(dbseq_element, events=("start", "end",), tag=("{*}cvParam", "{*}userParam")):
            if event == "start":
                param_name, param_val = MzidQuickReader._parse_param_name_and_value(item)
                if param_name is not None:
                    attributes[param_name] = param_val

            else:
                item.clear()

        # there is also Seq but ignore this for psm_utils
        return {dbseq_id: attributes}

    @staticmethod
    def _parse_searchdb(searchdb_element: etree.Element) -> dict[str, dict]:
        db_id = None

        # parse the SearchDB's attributes
        attributes = MzidQuickReader._parse_elements_attributes(searchdb_element)
        db_id = attributes["id"]
        del attributes["id"]

        # transform some types
        if "numDatabaseSequences" in attributes.keys():
            attributes["numDatabaseSequences"] = int(attributes["numDatabaseSequences"])
        if "numResidues" in attributes.keys():
            attributes["numResidues"] = int(attributes["numResidues"])
        
        for event, item in etree.iterwalk(searchdb_element, events=("start", "end"), tag=("{*}FileFormat", "{*}DatabaseName")):
            if event == "start":
                # strip the namespace
                tag = item.tag.rpartition("}")[2]

                # just take the name of the first userParam or cvParam in the FileFormat or DatabaseName
                for _, ff_item in etree.iterwalk(item, events=("end",), tag=("{*}cvParam", "{*}userParam")):
                    attributes[tag] = MzidQuickReader._parse_elements_attributes(ff_item)["name"]
                
                # there could also be cvParams, but ignore them for now
            else:            
                item.clear()
        
        return {db_id: attributes}

    @staticmethod
    def _parse_spectradata(spectradata_element: etree.Element) -> dict[str, dict]:
        specdata_id = None

        # parse the SearchDB's attributes
        attributes = MzidQuickReader._parse_elements_attributes(spectradata_element)
        specdata_id = attributes["id"]
        del attributes["id"]

        for event, item in etree.iterwalk(spectradata_element, events=("start", "end",), tag=("{*}FileFormat", "{*}SpectrumIDFormat")):
            if event == "start":
                # strip the namespace
                tag = item.tag.rpartition("}")[2]

                # just take the name of the first userParam or cvParam in the FileFormat or DatabaseName
                for _, item in etree.iterwalk(item, events=("end",), tag=("{*}cvParam", "{*}userParam")):
                    attributes[tag] = MzidQuickReader._parse_elements_attributes(item)["name"]
                
                # there could also be cvParams, but ignore them for now
            else:
                item.clear()
        
        return {specdata_id: attributes}

    def _parse_sir(self, sir_element: etree.Element) -> dict:
        # parse the SearchDB's attributes
        attributes = MzidQuickReader._parse_elements_attributes(sir_element)
        attributes["SpectrumIdentificationItem"] = []

        if "spectraData_ref" in attributes.keys() and attributes["spectraData_ref"] in self.spectra_data_dict.keys():
            spectra_data = self.spectra_data_dict[attributes["spectraData_ref"]]
            attributes |= spectra_data
            del attributes["spectraData_ref"]
        
        for event, item in etree.iterwalk(sir_element, events=("start", "end"), tag=("{*}SpectrumIdentificationItem", "{*}cvParam", "{*}userParam")):
            if event == "start":
                # strip the namespace
                tag = item.tag.rpartition("}")[2]

                if (tag == "SpectrumIdentificationItem"):
                    attributes["SpectrumIdentificationItem"].append(self._parse_sii(item))
                elif (tag == "cvParam") or (tag == "userParam"):
                    param_name, param_val = MzidQuickReader._parse_param_name_and_value(item)
                    if param_name is not None:
                        attributes[param_name] = param_val
            else:
                item.clear()
        
        return attributes
    
    def _parse_sii(self, sii_element: etree.Element) -> dict:
        # parse the SearchDB's attributes
        attributes = MzidQuickReader._parse_elements_attributes(sii_element)
        del attributes["id"]

        # transform some types
        if "calculatedMassToCharge" in attributes.keys():
            attributes["calculatedMassToCharge"] = float(attributes["calculatedMassToCharge"])
        if "calculatedPI" in attributes.keys():
            attributes["calculatedPI"] = float(attributes["calculatedPI"])
        if "chargeState" in attributes.keys():
            attributes["chargeState"] = int(attributes["chargeState"])
        if "experimentalMassToCharge" in attributes.keys():
            attributes["experimentalMassToCharge"] = float(attributes["experimentalMassToCharge"])
        if "passThreshold" in attributes.keys():
            attributes["passThreshold"] = MzidQuickReader._text_to_boolean(attributes["passThreshold"])
        if "rank" in attributes.keys():
            attributes["rank"] = int(attributes["rank"])
        
        # get the peptide information
        peptide_data = self.peptides_dict[attributes["peptide_ref"]]
        attributes |= peptide_data
        del attributes["peptide_ref"]
        
        attributes["PeptideEvidenceRef"] = []

        for event, item in etree.iterwalk(sii_element, events=("start", "end"), tag=("{*}PeptideEvidenceRef", "{*}cvParam", "{*}userParam")):
            if event == "start":
                # strip the namespace
                tag = item.tag.rpartition("}")[2]

                if (tag == "PeptideEvidenceRef"):
                    pep_evidence_data = self._parse_peptide_evidence_ref(item)
                    attributes["PeptideEvidenceRef"].append(pep_evidence_data)
                elif (tag == "cvParam") or (tag == "userParam"):
                    param_name, param_val = MzidQuickReader._parse_param_name_and_value(item)
                    if param_name is not None:
                        attributes[param_name] = param_val
            else:
                item.clear()
        
        return attributes
    
    @staticmethod
    def _parse_param_name_and_value(param_item: etree.Element) -> tuple[str, str | float]:
        param_name = None
        param_val = None
        param_attrs = MzidQuickReader._parse_elements_attributes(param_item)
        if "name" in param_attrs.keys() and "value" in param_attrs.keys():
            param_name = param_attrs["name"]
            try:
                param_val = float(param_attrs["value"])
            except ValueError:
                param_val = str(param_attrs["value"])
        
        return param_name, param_val

    def _parse_peptide_evidence_ref(self, pepevidenceref_item: etree.Element) -> dict[str, dict]:
        pep_evidence_attrs = MzidQuickReader._parse_elements_attributes(pepevidenceref_item)
        pep_evidence_data = copy.deepcopy(self.peptide_evidences_dict[pep_evidence_attrs["peptideEvidence_ref"]])

        # add actual information from DBSequence
        db_sequence_data = self.db_sequences_dict[pep_evidence_data["dBSequence_ref"]]
        pep_evidence_data |= db_sequence_data
        del pep_evidence_data["dBSequence_ref"]

        search_db_data = self.search_dbs_dict[pep_evidence_data["searchDatabase_ref"]]
        pep_evidence_data |= search_db_data
        del pep_evidence_data["searchDatabase_ref"]
        
        peptide_data = self.peptides_dict[pep_evidence_data["peptide_ref"]]
        pep_evidence_data |= peptide_data
        del pep_evidence_data["peptide_ref"]

        return pep_evidence_data


    @staticmethod
    def _parse_analysissoftware(spectradata_element: etree.Element) -> str:
        software_name = None
        for event, item in etree.iterwalk(spectradata_element, events=("start", "end",), tag=("{*}SoftwareName")):
            if event == "start":
                # just take the name of the first userParam or cvParam in the SoftwareName
                for _, sub_item in etree.iterwalk(item, events=("end",), tag=("{*}cvParam", "{*}userParam")):
                    software_name = MzidQuickReader._parse_elements_attributes(sub_item)["name"]
                
                # there can be other tags, not needed for now
            else:
                item.clear()
        
        return software_name

    @staticmethod
    def _parse_peptidoform(seq: str, modification_list: list[dict], charge: Union[int, None]):
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
        if charge:
            proforma_seq += f"/{charge}"

        return Peptidoform(proforma_seq)

    @staticmethod
    def _get_accessions_from_peptide_evidence_ref(peptide_evidence_list: list[dict]):
        """
        Parse PeptideEvidence list of PSM.

        Notes
        -----
        If multiple PeptideEvidence entries are associated with the PSM, the PSM is only considered
        a decoy entry if ALL PeptideEvidence entries are decoy entries. If a target PeptideEvidence
        entry is present, it should get priority over decoy entries. In theory, no overlap between
        target and decoy peptide sequence should be present in the search space, although this
        might not have been filtered for by the search engine.

        """
        isdecoy = all(
            [entry["isDecoy"] if "isDecoy" in entry else None for entry in peptide_evidence_list]
        )
        protein_list = [d["accession"] for d in peptide_evidence_list if "accession" in d.keys()]
        return isdecoy, protein_list

    def _get_peptide_spectrum_match(
        self,
        spectrum_id: str,
        spectrum_title: Union[str, None],
        run: Union[str, None],
        rt: Union[float, None],
        ion_mobility: Union[float, None],
        spectrum_identification_item: dict[str, str | float | list],
    ) -> PSM:
        """Parse single mzid entry to :py:class:`~psm_utils.peptidoform.Peptidoform`."""
        sii = spectrum_identification_item
        try:
            modifications = sii["Modification"]
        except KeyError:
            modifications = []
        sequence = sii["PeptideSequence"]
        charge = sii["chargeState"] if "chargeState" in sii else None
        peptidoform = self._parse_peptidoform(sequence, modifications, charge)
        is_decoy, protein_list = self._get_accessions_from_peptide_evidence_ref(sii["PeptideEvidenceRef"])
        try:
            precursor_mz = sii["experimentalMassToCharge"]
        except KeyError:
            precursor_mz = None
        
        # Override spectrum-level RT if present at PSM level
        if self._rt_key:
            rt = float(sii[self._rt_key])

        metadata = {col: str(sii[col]) for col in sii.keys() if col not in self._non_metadata_keys}

        # Prefer spectrum title (identical to titles in MGF), fall back to spectrumID
        if spectrum_title:
            metadata["mzid_spectrum_id"] = spectrum_id
            psm_spectrum_id = spectrum_title
        else:
            psm_spectrum_id = spectrum_id

        try:
            score = sii[self._score_key]
        except KeyError:
            score = None
        psm = PSM(
            peptidoform=peptidoform,
            spectrum_id=psm_spectrum_id,
            run=run,
            is_decoy=is_decoy,
            score=score,
            qvalue=sii[self._qvalue_key] if self._qvalue_key else None,
            pep=sii[self._pep_key] if self._pep_key else None,
            precursor_mz=precursor_mz,
            retention_time=rt,
            ion_mobility=ion_mobility,
            protein_list=protein_list,
            rank=sii["rank"] if "rank" in sii else None,
            source=self._source,
            provenance_data={"mzid_filename": str(self.filename)},
            metadata=metadata,
        )
        return psm

    def _get_non_metadata_keys(self, keys: list):
        """Gather all the keys at PSM-level that should not be written to metadata."""
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
        if not self._score_key:
            self._score_key = self._infer_score_name(keys)
        if self._score_key:
            default_keys.append(self._score_key)
        else:
            logger.warning(
                "No known score metric found in mzIdentML file. Scores will be set to None."
            )

        # Get the q-value key and add to default keys
        self._qvalue_key = self._infer_qvalue_name(keys)
        if self._qvalue_key:
            default_keys.append(self._qvalue_key)

        # Get the PEP key and add to default keys
        self._pep_key = self._infer_pep_name(keys)
        if self._pep_key:
            default_keys.append(self._pep_key)

        # Get retention time key
        for rt_key in ["retention time", "scan start time"]:
            if rt_key in keys:
                self._rt_key = rt_key
                default_keys.append(rt_key)
                break

        # Keys that are not necessary for metadata
        self._non_metadata_keys.extend(default_keys)

    def _get_toplevel_non_metadata_keys(self, keys: list):
        """Gather all keys at spectrum-level that should not be written to metadata."""
        # Check if RT is encoded in spectrum metadata
        for key in ["retention time", "scan start time"]:
            if key in keys:
                self._spectrum_rt_key = key
                self._non_metadata_keys.append(key)
                break

        # Check if ion mobility is encoded in spectrum metadata
        for im_key in ["inverse reduced ion mobility"]:
            if im_key in keys:
                self._im_key = im_key
                self._non_metadata_keys.append(im_key)
                break

    @staticmethod
    def _infer_score_name(keys) -> str:
        """Infer the score from the list of known PSM scores."""
        lower_keys = {key.lower(): key for key in keys}
        for score in STANDARD_SEARCHENGINE_SCORES:
            if score in lower_keys:
                return lower_keys[score]

    @staticmethod
    def _infer_qvalue_name(keys) -> Union[str, None]:
        """Infer the q-value term from the list of known terms."""
        for qvalue in Q_VALUE_TERMS:
            if qvalue in keys:
                return qvalue
        else:
            return None

    @staticmethod
    def _infer_pep_name(keys) -> Union[str, None]:
        """Infer the PEP term from the list of known terms."""
        for pep in PEP_TERMS:
            if pep in keys:
                return pep
        else:
            return None


class MzidWriter(WriterBase):
    """Writer for MzIdentMl files."""

    def __init__(
        self,
        filename: str | Path,
        *args,
        show_progressbar: bool = False,
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
        - Unlike other psm_utils.io writer classes, :py:class:`MzidWriter` does not support writing
          a single PSM to a file with the :py:meth:`write_psm` method. Only writing a full PSMList
          to a file at once with the :py:meth:`write_file` method is currently supported.
        - While not required according to the mzIdentML specification document (v1.1.1), the
          retention time is written as cvParam ``retention time`` to the
          SpectrumIdentificationItem element. As the actual unit is not known in psm_utils,
          the unit is written as seconds.
        - As the actual PSM score type is not known in psm_utils, the score is written as cvParam
          ``MS:1001153`` to the SpectrumIdentificationItem element.

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
        with Progress(disable=not self.show_progressbar) as progress:
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
                    if prot_list
                    for prot in prot_list
                }

                spec_id_dict = psm_list.get_psm_dict()
                task1 = progress.add_task("[cyan]Writing Proteins to mzid", total=len(proteins))
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

                        if psm["protein_list"]:
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
                        # # if fasta file is given, we can parse here and add protein information
                        # search_databases=transform_search_database(),
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
                                    identified_psms = spec_id_dict[collection][run][spec_id]
                                    writer.write_spectrum_identification_result(
                                        **self._transform_spectrum_identification_result(
                                            spec_id, identified_psms, spectra_data_id
                                        )
                                    )
                                    progress.update(
                                        task3,
                                        advance=len(spec_id_dict[collection][run][spec_id]),
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
            parse_modifications(peptidoform.properties["c_term"], len(peptidoform.sequence) + 1)
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
    def _transform_spectrum_identification_item(candidate_psm):
        """Create SpectrumIdentificationItem for each candidate PSM."""
        peptide = candidate_psm["peptidoform"].proforma
        if candidate_psm["metadata"]:
            params = [{k: v} for k, v in candidate_psm["metadata"].items()]
        else:
            params = []

        for key, label, unit in [
            ("retention_time", "retention time", "second"),
            ("qvalue", "PSM-level q-value", None),
            ("pep", "PSM-level local FDR", None),
        ]:
            if candidate_psm[key]:
                if unit:
                    params.append({"name": label, "value": candidate_psm[key], "unit_name": unit})
                else:
                    params.append({"name": label, "value": candidate_psm[key]})

        candidate_psm_dict = {
            "charge_state": candidate_psm["peptidoform"].precursor_charge,
            "peptide_id": f"Peptide_{peptide}",
            "score": {"search engine specific score": candidate_psm["score"]},
            "experimental_mass_to_charge": candidate_psm["precursor_mz"],
            "calculated_mass_to_charge": candidate_psm["peptidoform"].theoretical_mz,
            "rank": candidate_psm["rank"],
            "params": params,
        }
        items = []
        if candidate_psm["protein_list"]:
            for prot_acc in candidate_psm["protein_list"]:
                protein_specific_items = {
                    "peptide_evidence_id": f"PeptideEvidence_{peptide}_{prot_acc}",
                    "id": f"SII_{candidate_psm['spectrum_id']}_{peptide}_{prot_acc}",
                }
                items.append(dict(candidate_psm_dict, **protein_specific_items))
        return items

    def _transform_spectrum_identification_result(self, spec_id, identified_psms, spectra_data_id):
        """Create mzid SpectrumIdentificationResult object for PSMs that match the same spectrum."""
        spectrum_id_result = {
            "id": f"SIR_{spec_id}",
            "spectrum_id": spec_id,
            "spectra_data_id": spectra_data_id,
        }
        identifications = []
        for candidate in identified_psms:
            identifications.extend(self._transform_spectrum_identification_item(candidate))
        spectrum_id_result["identifications"] = identifications
        return spectrum_id_result


class UnknownMzidScore(PSMUtilsIOException):
    """No known score metric found in mzIdentML file."""
