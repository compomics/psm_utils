"""
Reader for PSM files from the AlphaDIA search engine.

Reads the AlphaDIA ``precursor.tsv`` file as defined on the
`TODO: NOT YET A LINK`_.

"""

from __future__ import annotations

import csv
from abc import ABC
from typing import Iterable, Optional

from psm_utils.io._base_classes import ReaderBase
from psm_utils.io._utils import set_csv_field_size_limit
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

set_csv_field_size_limit()


class AlphaDIAReader(ReaderBase, ABC):
    def __init__(self, filename, score_column: str = "score", *args, **kwargs):
        """
        Reader for AlphaDIA ``precursor.tsv`` file.

        Parameters
        ----------
        filename : str or Path
            Path to PSM file.
        score_column: str, optional
            Name of the column that holds the primary PSM score. Default is
            ``score``.

        """
        super().__init__(filename, *args, **kwargs)
        self.filename = filename
        self.score_column = score_column

    def __iter__(self) -> Iterable[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        with open(self.filename) as msms_in:
            reader = csv.DictReader(msms_in, delimiter="\t")
            for row in reader:
                yield self._get_peptide_spectrum_match(row)

    def _get_peptide_spectrum_match(self, psm_dict) -> PSM:
        """Parse a single PSM from a AlphaDIA PSM file."""
        rescoring_features = {}
        for ft in RESCORING_FEATURES:
            try:
                rescoring_features[ft] = psm_dict[ft]
            except KeyError:
                continue

        return PSM(
            peptidoform=self._parse_peptidoform(
                psm_dict["sequence"], psm_dict["mods"], psm_dict["mod_sites"], psm_dict["charge"]
            ),
            spectrum_id=psm_dict["frame_start"],  # TODO: needs to be checked
            run=psm_dict["run"],
            spectrum=psm_dict["frame_start"],  # TODO: needs to be checked
            is_decoy=bool(int(psm_dict["decoy"])),
            score=psm_dict[self.score_column],
            qvalue=psm_dict["qval"],
            pep=psm_dict[
                "proba"
            ],  # TODO: needs to be checked, assumption because if it is 1-proba than it's really bad
            precursor_mz=psm_dict["mz_observed"],
            retention_time=psm_dict["rt_observed"],
            ion_mobility=psm_dict["mobility_observed"],
            protein_list=psm_dict["proteins"].split(";"),
            rank=int(psm_dict["rank"]) + 1,  # AlphaDIA ranks are 0-based
            source="alphadia",
            provenance_data=({"alphadia_filename": str(self.filename)}),
            metadata={},
            rescoring_features=rescoring_features,
        )

    @staticmethod
    def _parse_peptidoform(sequence: str, mods: str, mod_sites, charge: Optional[str]) -> str:
        if mods:
            mods = mods.split(";")
            mod_sites = mod_sites.split(";")
            for mod, site in reversed(sorted(zip(mods, mod_sites), key=lambda x: int(x[1]))):
                if int(site) == 0:
                    sequence = (
                        sequence[: int(site)] + f"[{mod.split('@')[0]}]-" + sequence[int(site) :]
                    )
                else:
                    sequence = (
                        sequence[: int(site)] + f"[{mod.split('@')[0]}]" + sequence[int(site) :]
                    )
        if charge:
            sequence += f"/{int(float(charge))}"
        return sequence

    @classmethod
    def from_dataframe(cls, dataframe) -> PSMList:
        """Create a PSMList from a AlphaDIA Pandas DataFrame."""
        return PSMList(
            psm_list=[
                cls._get_peptide_spectrum_match(cls(""), entry)
                for entry in dataframe.to_dict(orient="records")
            ]
        )


# TODO: check
RESCORING_FEATURES = [
    "rt_observed",
    "mobility_observed",
    "mz_observed",
    "score",
    "charge",
    "delta_rt",
]
