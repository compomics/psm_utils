"""Interface with TPP pepXML PSM files."""

from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path
from typing import List, Optional, Union

from pyteomics import pepxml, proforma

from psm_utils.io._base_classes import ReaderBase
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PSM
from psm_utils.utils import mass_to_mz

logger = logging.getLogger(__name__)

STANDARD_SEARCHENGINE_SCORES = [
    "expect",
    "EValue",
    "Evalue",
    "SpecEValue",
    "xcorr",  # Fallback if no e-value is present
    "delta_dot",  # SpectraST
    "mzFidelity",
]


class PepXMLReader(ReaderBase):
    def __init__(self, filename: Union[str, Path], *args, score_key: str = None, **kwargs) -> None:
        """
        Reader for pepXML PSM files.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to PSM file.
        score_key: str, optional
            Name of the score metric to use as PSM score. If not provided, the score metric is
            inferred from a list of known search engine scores.

        """
        super().__init__(filename, *args, **kwargs)
        self.score_key = score_key or self._infer_score_name()

    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""
        with pepxml.read(str(self.filename)) as reader:
            for spectrum_query in reader:
                if "search_hit" not in spectrum_query:
                    continue
                for search_hit in spectrum_query["search_hit"]:
                    yield self._parse_psm(spectrum_query, search_hit)

    def _infer_score_name(self) -> str:
        """Infer the score from the list of known PSM scores."""
        # Get scores from first PSM
        with pepxml.read(str(self.filename)) as reader:
            for spectrum_query in reader:
                score_keys = spectrum_query["search_hit"][0]["search_score"].keys()
                break

        # Infer score name
        if not score_keys:
            logger.warning("No pepXML scores found.")
            return None
        else:
            for score in STANDARD_SEARCHENGINE_SCORES:  # Check for known scores
                if score in score_keys:
                    logger.debug(f"Using known pepXML score `{score}`.")
                    return score
            else:
                logger.warning(f"No known pepXML scores found. Defaulting to `{score_keys[0]}`.")
                return score_keys[0]  # Default to the first one if nothing found

    @staticmethod
    def _parse_peptidoform(peptide: str, modifications: List[dict], charge: Optional[int] = None):
        """Parse pepXML peptide to :py:class:`~psm_utils.peptidoform.Peptidoform`."""
        modifications_dict = defaultdict(list)
        n_term = []
        c_term = []
        for mod in modifications:
            mod_tag = proforma.process_tag_tokens(f"{mod['mass']:+}")
            if mod["position"] == 0:
                n_term.append(mod_tag)
            elif mod["position"] == len(peptide) + 1:
                c_term.append(mod_tag)
            else:
                modifications_dict[mod["position"]].append(mod_tag)

        sequence = [(aa, modifications_dict[i] or None) for i, aa in enumerate(peptide)]
        properties = {
            "n_term": n_term,
            "c_term": c_term,
            "charge_state": proforma.ChargeState(charge) if charge else None,
            "unlocalized_modifications": [],
            "labile_modifications": [],
            "fixed_modifications": [],
            "intervals": [],
            "isotopes": [],
            "group_ids": [],
        }
        return Peptidoform(proforma.ProForma(sequence, properties))

    def _parse_psm(self, spectrum_query: dict, search_hit: dict) -> PSM:
        """Parse pepXML PSM to :py:class:`~psm_utils.psm.PSM`."""
        return PSM(
            peptidoform=self._parse_peptidoform(
                search_hit["peptide"],
                search_hit["modifications"],
                spectrum_query["assumed_charge"],
            ),
            spectrum_id=spectrum_query["spectrum"],
            run=None,
            collection=None,
            spectrum=None,
            is_decoy=None,
            score=search_hit["search_score"][self.score_key],
            qvalue=None,
            pep=None,
            precursor_mz=mass_to_mz(
                spectrum_query["precursor_neutral_mass"], spectrum_query["assumed_charge"]
            ),
            retention_time=spectrum_query["retention_time_sec"]
            if "retention_time_sec" in spectrum_query
            else None,
            ion_mobility=spectrum_query["ion_mobility"]
            if "ion_mobility" in spectrum_query
            else None,
            protein_list=[p["protein"] for p in search_hit["proteins"]],
            rank=search_hit["hit_rank"],
            source=None,
            provenance_data={
                "pepxml_index": str(spectrum_query["index"]),
                "start_scan": str(spectrum_query["start_scan"]),
                "end_scan": str(spectrum_query["end_scan"]),
            },
            metadata={
                "num_matched_ions": str(search_hit["num_matched_ions"]),
                "tot_num_ions": str(search_hit["tot_num_ions"]),
                "num_missed_cleavages": str(search_hit["num_missed_cleavages"]),
            }.update(
                {
                    f"search_score_{key.lower()}": str(search_hit["search_score"][key])
                    for key in search_hit["search_score"]
                }
            ),
            rescoring_features=None,
        )
