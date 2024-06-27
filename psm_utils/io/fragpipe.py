"""
Reader for PSM files from the Fragpipe platform.

Reads the Philosopher ``psm.tsv`` file as defined on the
`Fragpipe documentation page <https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html>`_.

"""

from __future__ import annotations

import csv
from abc import ABC
from pathlib import Path
from typing import Iterable, Optional

from psm_utils.io._base_classes import ReaderBase
from psm_utils.io._utils import set_csv_field_size_limit
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

set_csv_field_size_limit()


class FragpipeReader(ReaderBase, ABC):
    def __init__(
        self,
        filename,
        score_column: str = "Hyperscore",
        mz_column: str = "Observed M/Z",
        *args,
        **kwargs,
    ) -> None:
        """
        Reader for MSFragger ``psm.tsv`` file.

        Parameters
        ----------
        filename : str or Path
            Path to PSM file.
        score_column: str, optional
            Name of the column that holds the primary PSM score. Default is
            ``Hyperscore``.
        mz_column: str, optional
            Name of the column that holds the precursor m/z. Default is
            ``Observed M/Z``.

        """
        super().__init__(filename, *args, **kwargs)
        self.filename = filename
        self.score_column = score_column
        self.mz_column = mz_column

    def __iter__(self) -> Iterable[PSM]:
        """Iterate over file and return PSMs one-by-one."""
        with open(self.filename) as msms_in:
            reader = csv.DictReader(msms_in, delimiter="\t")
            for row in reader:
                yield self._get_peptide_spectrum_match(row)

    def _get_peptide_spectrum_match(self, psm_dict) -> PSM:
        """Parse a single PSM from a MSFragger PSM file."""
        rescoring_features = {}
        for ft in RESCORING_FEATURES:
            try:
                rescoring_features[ft] = psm_dict[ft]
            except KeyError:
                continue

        return PSM(
            peptidoform=self._parse_peptidoform(
                psm_dict["Modified Peptide"], psm_dict["Peptide"], psm_dict["Charge"]
            ),
            spectrum_id=self._parse_spectrum_id(psm_dict["Spectrum"]),  # TODO: needs to be checked
            run=self._parse_run(psm_dict["Spectrum File"]),
            is_decoy=False,
            qvalue=None,  # Q-value is not outputted by Philosopher
            pep=1
            - float(
                psm_dict["Probability"]
            ),  # PeptideProphet Probability, not explicitely stated if this is the inverse of PEP
            # But I'm assuming it is
            score=psm_dict[self.score_column],
            precursor_mz=psm_dict[
                self.mz_column
            ],  # Allows use of both calibrated and uncalibrated Observed M/Z?
            retention_time=float(psm_dict["Retention"]),
            ion_mobility=float(psm_dict["Ion Mobility"]) if "Ion Mobility" in psm_dict else None,
            protein_list=self._parse_protein_list(
                psm_dict["Protein"], psm_dict["Mapped Proteins"]
            ),
            source="fragpipe",
            rank=1,
            provenance_data=({"fragpipe_filename": str(self.filename)}),
            rescoring_features=rescoring_features,
            metadata={},
        )

    @staticmethod
    def _parse_peptidoform(mod_peptide: str, peptide: str, charge: Optional[str]) -> str:
        if mod_peptide:
            peptide = mod_peptide
        if charge:
            peptide += f"/{int(float(charge))}"
        if peptide.startswith("n"):
            peptide = peptide[1:]
            # A hyphen needs to be added after the N-terminal modification, thus after the ]
            peptide = peptide.replace("]", "]-", 1)
        return peptide

    @staticmethod
    def _parse_spectrum_id(spectrum: str) -> str:
        return spectrum.split(".")[1]

    @staticmethod
    def _parse_protein_list(razor_protein: str, mapped_proteins) -> list[str]:
        if mapped_proteins:
            mapped_proteins_list = mapped_proteins.split(", ")
            return [razor_protein] + mapped_proteins_list
        else:
            return [razor_protein]

    # Dependent on the fragpipe workflow used the run name can be different, but in most cases
    # something like 'interact-<run_name>.pep.xml' is used
    @staticmethod
    def _parse_run(spectrum_file: str) -> str:
        if (spectrum_file.endswith(".pep.xml")) and (spectrum_file.startswith("interact-")):
            spectrum_file = spectrum_file.replace("interact-", "")
            return Path(Path(spectrum_file).stem).stem
        else:
            return Path(spectrum_file).stem

    @classmethod
    def from_dataframe(cls, dataframe) -> PSMList:
        """Create a PSMList from a pandas DataFrame."""
        return PSMList(
            ptm_list=[
                cls._get_peptide_spectrum_match(cls(""), entry)
                for entry in dataframe.to_dict(orient="records")
            ]
        )


# TODO: check
RESCORING_FEATURES = [
    "Peptide Length",
    "Retention",
    "Observed Mass",
    "Observed M/Z",
    "Calculated Peptide Mass",
    "Calculated M/Z",
    "Delta Mass",
    "Hyperscore",
    "Number of Missed Cleavages",
]
