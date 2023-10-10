"""
Interface with OpenMS idXML PSM files.


Notes
-----

* idXML supports multiple peptide hits (identifications) per spectrum. Each peptide hit
  is parsed as an individual :py:class:`~psm_utils.psm.PSM` object.

"""


from __future__ import annotations

from abc import abstractmethod
from typing import Union, Iterable, Tuple, List
import re
from pathlib import Path

import pyopenms as oms

from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList
from psm_utils.exceptions import PSMUtilsException

# Patterns to match open and closed round/square brackets
MOD_PATTERN = re.compile(r"\(((?:[^)(]+|\((?:[^)(]+|\([^)(]*\))*\))*)\)")
MOD_PATTERN_NTERM = re.compile(r"^\.\[((?:[^][]+|\[(?:[^][]+|\[[^][]*\])*\])*)\]")
MOD_PATTERN_CTERM = re.compile(r"\.\[((?:[^][]+|\[(?:[^][]+|\[[^][]*\])*\])*)\]$")

# Extracted from the OpenMS PSMFeatureExtractor, which adds and manipulates features that will be given to percolator
# https://github.com/OpenMS/OpenMS/blob/342f6524e76a2bab3dcb428ba2f4aa2d6bfe8483/src/topp/PSMFeatureExtractor.cpp
RESCORING_FEATURE_LIST = [
    "isotope_error"
    "MS:1002049",  # MSGFPlus unchanged RawScore
    "MS:1002050",  # MSGFPlus unchanged DeNovoScore
    "MSGF:ScoreRatio",
    "MSGF:Energy",
    "MSGF:lnEValue",
    "MSGF:lnExplainedIonCurrentRatio",
    "MSGF:lnNTermIonCurrentRatio",
    "MSGF:lnCTermIonCurrentRatio",
    "MSGF:lnMS2IonCurrent",
    "MSGF:MeanErrorTop7",
    "MSGF:sqMeanErrorTop7",
    "MSGF:StdevErrorTop7",
    "XTANDEM:hyperscore",
    "XTANDEM:deltascore",
    "MS:1001330",  # expect_score
    "hyperscore",  # MSfragger
    "nextscore",  # MSfragger
    "COMET:deltaCn",  # recalculated deltaCn = (current_XCorr - 2nd_best_XCorr) / max(current_XCorr, 1)
    "COMET:deltaLCn",  # deltaLCn = (current_XCorr - worst_XCorr) / max(current_XCorr, 1)
    "COMET:lnExpect",  # log(E-value)
    "MS:1002252",  # unchanged XCorr
    "MS:1002255",  # unchanged Sp = number of candidate peptides
    "COMET:lnNumSP",  # log(number of candidate peptides)
    "COMET:lnRankSP",  # log(rank based on Sp score)
    "COMET:IonFrac",  # matched_ions / total_ions
    "MS:1001171",  # unchanged mScore
    "MASCOT:delta_score",  # delta score based on mScore
    "CONCAT:lnEvalue",
    "CONCAT:deltaLnEvalue",
    "SAGE:ln(-poisson)"
    "SAGE:ln(delta_best)",
    "SAGE:ln(delta_next)",
    "SAGE:ln(matched_intensity_pct)",
    "SAGE:longest_b",
    "SAGE:longest_y",
    "SAGE:longest_y_pct",
    "SAGE:matched_peaks",
    "SAGE:scored_candidates"
]

class IdXMLReader(ReaderBase):
    def __init__(self, filename: Union[Path, str], *args, **kwargs) -> None:
        """
        Reader for idXML files.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to idXML file.

        Examples
        --------
        >>> from psm_utils.io import IdXMLReader
        >>> reader = IdXMLReader("example.idXML")
        >>> psm_list = [psm for psm in reader]
        """
        super().__init__(filename, *args, **kwargs)
        self.filename = Path(filename)
        self.rescoring_features = None

    def __iter__(self) -> Iterable[PSM]:
        """
        Iterate over file and return PSMs one-by-one.
        """
        protein_ids, peptide_ids = self._parse_idxml()
        self.rescoring_features = self._get_rescoring_features(peptide_ids[0].getHits()[0])
        for peptide_id in peptide_ids:
            for peptide_hit in peptide_id.getHits():
                yield self._parse_psm(protein_ids, peptide_id, peptide_hit)

    def _parse_idxml(self) -> Tuple(oms.ProteinIdentification, oms.PeptideIdentification):
        """
        Parse idXML using pyopenms and do some sanity checks to make sure the file is not empty.
        """
        protein_ids, peptide_ids = [], []
        oms.IdXMLFile().load(str(self.filename), protein_ids, peptide_ids)

        if len(protein_ids) == 0:
            raise IdXMLReaderEmptyListException(f"File {self.filename} contains no proteins. Nothing to parse.")
        elif len(peptide_ids) == 0:
            raise IdXMLReaderEmptyListException(f"File {self.filename} contains no PeptideIdentifications. Nothing to parse.")
        elif len(peptide_ids[0].getHits()) == 0:
            raise IdXMLReaderEmptyListException(f"File {self.filename} contains no PeptideHits. Nothing to parse.")
        else:
            return protein_ids, peptide_ids

    def _parse_peptidoform(self, sequence: str, charge: int) -> str:
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


    def _parse_psm(self, protein_ids: oms.ProteinIdentification, peptide_id: oms.PeptideIdentification, peptide_hit: oms.PeptideHit) -> PSM:
        """
        Parse idXML :py:class:`~pyopenms.PeptideHit` to :py:class:`~psm_utils.psm.PSM`.
        Use additional information from :py:class:`~pyopenms.ProteinIdentification` and :py:class:`~pyopenms.PeptideIdentification`
        to annotate parameters of the :py:class:`~psm_utils.psm.PSM` object.
        """
        return PSM(
            peptidoform = self._parse_peptidoform(peptide_hit.getSequence().toString(), peptide_hit.getCharge()),
            spectrum_id = peptide_id.getMetaValue("spectrum_reference"),
            run = self._get_run(protein_ids, peptide_id),
            is_decoy = self._is_decoy(peptide_hit),
            score = peptide_hit.getScore(),
            precursor_mz = peptide_id.getMZ(),
            retention_time = peptide_id.getRT(),
            # NOTE: In future update of OpenMS we will be able to read out ion_mobility from idXML file
            protein_list=[accession.decode() for accession in peptide_hit.extractProteinAccessionsSet()],
            rank = peptide_hit.getRank() + 1, # 0-based to 1-based
            source="idXML",
            # This is needed to calculate a qvalue before rescoring the PSMList
            metadata = {
                "idxml:score_type": str(peptide_id.getScoreType()),
                "idxml:higher_score_better": str(peptide_id.isHigherScoreBetter()),
                "idxml:significance_threshold": str(peptide_id.getSignificanceThreshold())
            },
            rescoring_features = {key: float(peptide_hit.getMetaValue(key)) for key in self.rescoring_features \
                                  if self._is_float(peptide_hit.getMetaValue(key))}
        )

    def _get_run(self, protein_ids: oms.ProteinIdentification, peptide_id: oms.PeptideIdentification) -> str:
        """
        Get run name from idXML using pyopenms
        If idXML contains merge index, use it to annotate the run name without extension
        """
        if peptide_id.metaValueExists("id_merge_index"):
            return Path(protein_ids[0].getMetaValue("spectra_data")[peptide_id.getMetaValue("id_merge_index")].decode()).stem
        else:
            return Path(protein_ids[0].getMetaValue("spectra_data")[0].decode()).stem

    def _get_rescoring_features(self, peptide_hit: oms.PeptideHit) -> List[str]:
        """
        Get list of rescoring features from idXML
        """
        # Fill the key list with all the keys from the PeptideHit in the pyopenms way
        keys = []
        peptide_hit.getKeys(keys)
        # Convert from byte to str and remove the keys, which values are not in the fixed list of features
        keys = [key.decode() for key in keys if  key.decode() in RESCORING_FEATURE_LIST]

        return keys

    def _is_float(self, element: any) -> bool:
        """
        Check if element is float.
        """
        if element is None: 
            return False
        try:
            float(element)
            return True
        except ValueError:
            return False
        
    def _is_decoy(self, peptide_hit: oms.PeptideHit) -> bool:
        """
        Check if PSM is target or decoy
        """
        if peptide_hit.metaValueExists("target_decoy"):
            return peptide_hit.getMetaValue("target_decoy") == "decoy"
        else:
            return None

class IdXMLReaderEmptyListException(PSMUtilsException):
    """Exception in psm_utils.io.IdXMLReader"""
