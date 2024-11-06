"""Reader for PSM files from the AlphaDIA search engine."""

from __future__ import annotations

import csv
from abc import ABC
from typing import Iterable, Optional

from psm_utils.io._base_classes import ReaderBase
from psm_utils.io._utils import set_csv_field_size_limit
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

set_csv_field_size_limit()

# TODO: check
RESCORING_FEATURES = [
    "rt_observed",
    "mobility_observed",
    "mz_observed",
    "charge",
    "delta_rt",
]


class AlphaDIAReader(ReaderBase, ABC):
    def __init__(self, filename, *args, **kwargs):
        """
        Reader for AlphaDIA ``precursor.tsv`` file.

        Parameters
        ----------
        filename : str or Path
            Path to PSM file.

        """
        super().__init__(filename, *args, **kwargs)
        self.filename = filename

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
            score=psm_dict["score"],
            qvalue=psm_dict["qval"],
            pep=psm_dict["proba"],
            precursor_mz=psm_dict["mz_observed"],
            retention_time=psm_dict["rt_observed"],
            ion_mobility=psm_dict["mobility_observed"],
            protein_list=psm_dict["proteins"].split(";"),
            rank=int(psm_dict["rank"]) + 1,  # AlphaDIA ranks are 0-based
            source="AlphaDIA",
            provenance_data=({"alphadia_filename": str(self.filename)}),
            metadata={},
            rescoring_features=rescoring_features,
        )

    @staticmethod
    def _parse_peptidoform(sequence: str, mods: str, mod_sites, charge: Optional[str]) -> str:
        """Parse a peptidoform from a AlphaDIA PSM file."""
        # Parse modifications
        if mods:
            sequence_list = [""] + list(sequence) + [""]  # N-term, sequence, C-term
            for mod, site in zip(mods.split(";"), mod_sites.split(";")):
                site = int(site)
                name = mod.split("@")[0]
                # N-terminal modification
                if site == 0:
                    sequence_list[0] = f"[{name}]-"
                # C-terminal modification
                elif site == -1:
                    sequence_list[-1] = f"-[{name}]"
                # Sequence modification
                else:
                    sequence_list[site] = f"{sequence_list[site]}[{name}]"
            sequence = "".join(sequence_list)

        # Add charge
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
