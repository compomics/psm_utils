"""Reader for Proteome Discoverer MSF PSM files."""

import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Union

import pyteomics.proforma as proforma
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import psm_utils.io._pd_msf_tables as msf
from psm_utils import PSM, Peptidoform
from psm_utils.io._base_classes import ReaderBase

logger = logging.getLogger(__name__)

COMPATIBLE_VERSIONS = [79]


class MSFReader(ReaderBase):
    """Reader for Proteome Discoverer MSF files."""

    def __init__(
        self,
        filename: Union[str, Path],
        *args,
        **kwargs,
    ) -> None:
        """
        Reader for Proteome Discoverer MSF file.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to MSF file.

        """
        super().__init__(filename, *args, **kwargs)

        self._engine = create_engine(f"sqlite:///{self.filename.as_posix()}")
        self._session = sessionmaker(bind=self._engine)()

        self._check_version()

    def __len__(self):
        """Return number of PSMs in file."""
        return sum(
            self._session.query(peptide).count() for peptide in [msf.Peptide, msf.PeptideDecoy]
        )

    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""
        for is_decoy in [False, True]:
            modifications = self._get_modifications(is_decoy)
            terminal_modifications = self._get_terminal_modifications(is_decoy)
            protein_entries = self._get_protein_entries(is_decoy)
            main_score = self._get_main_score(is_decoy)
            secondary_scores = self._get_secondary_scores(is_decoy)

            for entry in self._iter_peptides(is_decoy):
                peptide_id = entry.PeptideDecoy.PeptideID if is_decoy else entry.Peptide.PeptideID
                yield self._parse_entry(
                    entry,
                    modifications[peptide_id],
                    terminal_modifications[peptide_id],
                    protein_entries[peptide_id],
                    main_score[peptide_id],
                    secondary_scores[peptide_id],
                    is_decoy,
                )

    def _check_version(self):
        """Check if MSF file version is compatible."""
        version = self._session.query(msf.SchemaInfo.Version).first()[0]
        if version not in COMPATIBLE_VERSIONS:
            logger.warning(
                f"MSF file version {version} might not be compatible with this reader. "
                f"Checked versions are: {COMPATIBLE_VERSIONS}."
            )

    def _iter_peptides(self, is_decoy: bool):
        """Iterate over peptides in MSF file."""
        Peptide = msf.PeptideDecoy if is_decoy else msf.Peptide
        for entry in (
            self._session.query(Peptide, msf.SpectrumHeader, msf.MassPeak, msf.FileInfo)
            .select_from(Peptide)
            .join(msf.SpectrumHeader, Peptide.SpectrumID == msf.SpectrumHeader.SpectrumID)
            .join(msf.MassPeak, msf.MassPeak.MassPeakID == msf.SpectrumHeader.MassPeakID)
            .join(msf.FileInfo, msf.FileInfo.FileID == msf.MassPeak.FileID)
        ):
            yield entry

    def _get_modifications(self, is_decoy: bool) -> Dict[int, Tuple[int, int]]:
        """Get all modifications per peptide ID."""
        PeptidesAminoAcidModification = (
            msf.PeptidesAminoAcidModificationsDecoy
            if is_decoy
            else msf.PeptidesAminoAcidModification
        )
        query = (
            self._session.query(
                PeptidesAminoAcidModification.PeptideID,
                PeptidesAminoAcidModification.Position,
                msf.AminoAcidModification.UnimodAccession,
            )
            .select_from(PeptidesAminoAcidModification)
            .join(
                msf.AminoAcidModification,
                PeptidesAminoAcidModification.AminoAcidModificationID
                == msf.AminoAcidModification.AminoAcidModificationID,
            )
        )
        modifications_by_peptide = defaultdict(list)
        for peptide_id, position, unimod_accession in query:
            modifications_by_peptide[peptide_id].append((position, unimod_accession))

        return modifications_by_peptide

    def _get_terminal_modifications(self, is_decoy: bool) -> Dict[int, Tuple[int, int]]:
        """Get terminal modifications for a peptide."""
        PeptidesTerminalModification = (
            msf.PeptidesTerminalModification if is_decoy else msf.PeptidesTerminalModificationDecoy
        )
        query = (
            self._session.query(
                PeptidesTerminalModification.PeptideID,
                msf.AminoAcidModification.PositionType,
                msf.AminoAcidModification.UnimodAccession,
            )
            .select_from(msf.AminoAcidModification)
            .join(
                PeptidesTerminalModification,
                PeptidesTerminalModification.TerminalModificationID
                == msf.AminoAcidModification.AminoAcidModificationID,
            )
        )
        terminal_modifications = defaultdict(list)
        for peptide_id, position_type, unimod_accession in query:
            terminal_modifications[peptide_id].append((position_type, unimod_accession))
        return terminal_modifications

    def _get_protein_entries(self, is_decoy: bool) -> Dict[int, List[str]]:
        """Get protein descriptions or a peptide."""
        PeptidesProtein = msf.PeptidesProteinDecoy if is_decoy else msf.PeptidesProtein
        query = (
            self._session.query(PeptidesProtein.PeptideID, msf.ProteinAnnotation.Description)
            .select_from(PeptidesProtein)
            .join(
                msf.ProteinAnnotation,
                PeptidesProtein.ProteinID == msf.ProteinAnnotation.ProteinID,
            )
        )
        proteins = defaultdict(list)
        for peptide_id, description in query:
            proteins[peptide_id].append(re.sub(r"^>", "", description))
        return proteins

    def _get_main_score(self, is_decoy: bool) -> Dict[int, Tuple[float, str]]:
        """Get main score and its name for a peptide."""
        PeptideScore = msf.PeptideScoreDecoy if is_decoy else msf.PeptideScore
        query = (
            self._session.query(
                PeptideScore.PeptideID, PeptideScore.ScoreValue, msf.ProcessingNodeScore.ScoreName
            )
            .select_from(PeptideScore)
            .join(
                msf.ProcessingNodeScore,
                msf.ProcessingNodeScore.ScoreID == PeptideScore.ScoreID,
            )
            .filter(msf.ProcessingNodeScore.IsMainScore == True)  # noqa: E712
        )
        scores = dict()
        for peptide_id, score_value, score_name in query:
            scores[peptide_id] = (score_value, score_name)
        return scores

    def _get_secondary_scores(self, is_decoy: bool) -> Dict[int, Dict[str, float]]:
        """Get secondary scores and their names for a peptide."""
        PeptideScore = msf.PeptideScoreDecoy if is_decoy else msf.PeptideScore
        query = (
            self._session.query(
                PeptideScore.PeptideID, PeptideScore.ScoreValue, msf.ProcessingNodeScore.ScoreName
            )
            .select_from(PeptideScore)
            .join(
                msf.ProcessingNodeScore,
                msf.ProcessingNodeScore.ScoreID == PeptideScore.ScoreID,
            )
            .filter(msf.ProcessingNodeScore.IsMainScore == False)  # noqa: E712
        )
        scores = defaultdict(dict)
        for peptide_id, score_value, score_name in query:
            scores[peptide_id][score_name] = score_value
        return scores

    def _compile_peptidoform(
        self,
        sequence: str,
        charge: int,
        modifications: List[Tuple[int, int]],
        terminal_modifications: List[Tuple[int, int]],
    ) -> Peptidoform:
        """
        Compile a peptidoform from a sequence, charge, and list of (terminal) modifications.

        Parameters
        ----------
        sequence
            The stripped sequence of the peptidoform.
        charge
            Precursor charge.
        modifications
            List of tuples of the form (position, unimod identifier).
        terminal_modifications
            List of tuples of the form (position type, unimod identifier).

        Notes
        -----
        The position type is either 1 (Any N-term), 2 (Any C-term), 3 (Protein N-term), or 4
        (Protein C-term). Position type 0 (Anywhere) should not be present in the
        terminal_modifications list.

        """
        modifications_dict = defaultdict(list)
        for position, unimod_id in modifications:
            modifications_dict[position].append(proforma.process_tag_tokens(f"U:{unimod_id}"))

        n_term = [
            proforma.process_tag_tokens(f"U:{unimod_id}")
            for position_type, unimod_id in terminal_modifications
            if position_type in [1, 3]  # Position types 'Any N-term' or 'Protein N-term'
        ]
        c_term = [
            proforma.process_tag_tokens(f"U:{unimod_id}")
            for position_type, unimod_id in terminal_modifications
            if position_type in [2, 4]  # Position types 'Any C-term' or 'Protein C-term'
        ]

        sequence = [(aa, modifications_dict[i] or None) for i, aa in enumerate(sequence)]
        properties = {
            "n_term": n_term,
            "c_term": c_term,
            "charge_state": proforma.ChargeState(charge),
            "unlocalized_modifications": [],
            "labile_modifications": [],
            "fixed_modifications": [],
            "intervals": [],
            "isotopes": [],
            "group_ids": [],
        }

        return Peptidoform(proforma.ProForma(sequence, properties))

    def _parse_entry(
        self,
        entry: Tuple[msf.Peptide, msf.SpectrumHeader, msf.MassPeak, msf.FileInfo],
        modifications: List[Tuple[int, int]],
        terminal_modifications: List[Tuple[int, int]],
        protein_entries: List[str],
        main_score: Tuple[float, str],
        secondary_scores: Dict[str, float],
        is_decoy: bool,
    ) -> PSM:
        """Parse an entry from the MSF file."""
        peptide = entry.PeptideDecoy if is_decoy else entry.Peptide
        return PSM(
            peptidoform=self._compile_peptidoform(
                peptide.Sequence,
                entry.SpectrumHeader.Charge,
                modifications,
                terminal_modifications,
            ),
            spectrum_id=entry.SpectrumHeader.LastScan,
            run=Path(entry.FileInfo.FileName).stem,
            is_decoy=is_decoy,
            score=main_score[0],
            qvalue=None,
            pep=None,
            precursor_mz=entry.MassPeak.Mass,
            retention_time=entry.SpectrumHeader.RetentionTime,
            ion_mobility=None,
            protein_list=protein_entries,
            rank=peptide.SearchEngineRank,
            source="proteome_discoverer",
            provenance_data={
                "scan_numbers": entry.SpectrumHeader.ScanNumbers,
            },
            metadata={
                "ms1_intensity": str(entry.MassPeak.Intensity),
                "ms1_percent_isolation_interference": str(
                    entry.MassPeak.PercentIsolationInterference
                ),
                "ms1_ion_inject_time": str(entry.MassPeak.IonInjectTime),
                "main_score_name": main_score[1],
                **secondary_scores,
            },
            rescoring_features={
                "missed_cleavages": peptide.MissedCleavages,
                "total_ions_count": peptide.TotalIonsCount,
                "matched_ions_count": peptide.MatchedIonsCount,
            },
        )
