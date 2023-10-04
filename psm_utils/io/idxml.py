"""
Interface with OpenMS idXML PSM files.


Notes
-----

* idXML supports multiple peptide hits (identifications) per spectrum. Each peptide hit
  is parsed as an individual :py:class:`~psm_utils.psm.PSM` object.

"""


from __future__ import annotations

from abc import abstractmethod
import re

try:
    import pyopenms as oms
except ImportError:
    from pyteomics.openms import idxml

from psm_utils.io._base_classes import ReaderBase
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

# Patterns to match open and closed round/square brackets
MOD_PATTERN = re.compile(r"\(((?:[^)(]+|\((?:[^)(]+|\([^)(]*\))*\))*)\)")
MOD_PATTERN_NTERM = re.compile(r"^\.\[((?:[^][]+|\[(?:[^][]+|\[[^][]*\])*\])*)\]")
MOD_PATTERN_CTERM = re.compile(r"\.\[((?:[^][]+|\[(?:[^][]+|\[[^][]*\])*\])*)\]$")

class AIdXMLReader(ReaderBase):
    @abstractmethod
    def __iter__(self):
        raise NotImplementedError()
    
    @abstractmethod
    def _parse_psm(self):
        raise NotImplementedError()
    
    @staticmethod
    def _parse_peptidoform(sequence: str, charge: int):
        """
        Parse idXML peptide to :py:class:`~psm_utils.peptidoform.Peptidoform`.
        Notes
        -----
        Implemented according to the documentation on
        `github.com/OpenMS/OpenMS <https://github.com/OpenMS/OpenMS/blob/8cb90/src/openms/include/OpenMS/CHEMISTRY/AASequence.h>`_
        . The differentiation between square- and round bracket notation is removed
        after parsing.
        """
        sequence = MOD_PATTERN.sub(r"[\1]", sequence)
        if sequence[:2] == ".[":
            sequence = MOD_PATTERN_NTERM.sub(r"[\1]-", sequence)
        if sequence[-1] == "]":
            sequence = MOD_PATTERN_CTERM.sub(r"-[\1]", sequence)
        sequence = sequence.strip(".")
        sequence += f"/{charge}"
        return sequence


class IdXMLReader(AIdXMLReader):
    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""
        with idxml.read(str(self.filename)) as reader:
            for entry in reader:
                for peptide_hit in entry["PeptideHit"]:
                    yield self._parse_psm(entry, peptide_hit)

    def _parse_psm(self, entry: dict, peptide_hit: dict) -> PSM:
        """Parse idXML PSM to :py:class:`~psm_utils.psm.PSM`."""
        return PSM(
            peptidoform=self._parse_peptidoform(peptide_hit["sequence"], peptide_hit["charge"]),
            spectrum_id=entry["spectrum_reference"],
            is_decoy= peptide_hit["target_decoy"] == "decoy",
            score=peptide_hit["score"],
            precursor_mz=entry["MZ"],
            retention_time=entry["RT"],
            protein_list=[protein["accession"] for protein in peptide_hit["protein"]],
            source="idXML",
            metadata={
                "idxml:score_type": str(entry["score_type"]),
                "idxml:higher_score_better": str(entry["higher_score_better"]),
                "idxml:significance_threshold": str(entry["significance_threshold"]),
            },
        )


class OmsIdXMLReader(AIdXMLReader):
    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""
        # TODO: This is a bit hacky, is there a cleaner way to do this?
        try:
            protein_ids, peptide_ids = [], []
            oms.IdXMLFile().load(str(self.filename), protein_ids, peptide_ids)
        except NameError:
            raise ModuleNotFoundError("Install pyopenms to use OmsIdXMLReader.")
        for peptide_id in peptide_ids:
            for peptide_hit in peptide_id.getHits():
                yield self._parse_psm(peptide_id, peptide_hit)

    def _parse_psm(self, peptide_id: oms.PeptideIdentification, peptide_hit: oms.PeptideHit) -> PSM:
        """Parse idXML PSM to :py:class:`~psm_utils.psm.PSM` using pyopenms"""
        return PSM(
            peptidoform = self._parse_peptidoform(peptide_hit.getSequence().toString(), peptide_hit.getCharge()),
            spectrum_id = peptide_id.getMetaValue("spectrum_reference"),
            is_decoy = peptide_hit.getMetaValue("target_decoy") == "decoy",
            score = peptide_hit.getScore(),
            precursor_mz = peptide_id.getMZ(),
            retention_time = peptide_id.getRT(),
            # TODO: ion mobility is currently not parsed to idXML so this needs to be read out from mzml
            protein_list=[accession.decode() for accession in peptide_hit.extractProteinAccessionsSet()],
            rank = peptide_hit.getRank() + 1, # 0-based to 1-based
            source="idXML",
            # TODO: Store all peptide hit related metadata?
            metadata = {
                "idxml:score_type": str(peptide_id.getScoreType()),
                "idxml:higher_score_better": str(peptide_id.isHigherScoreBetter()),
                "idxml:significance_threshold": str(peptide_id.getSignificanceThreshold())
            }
            # TODO: provenance_data
            # TODO: rescoring_features (Is there a fixed nomenclature for using these?)
        )