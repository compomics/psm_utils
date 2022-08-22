"""
Interface with OpenMS idXML PSM files.


Notes
-----

* idXML supports multiple peptide hits (identifications) per spectrum. Each peptide hit
  is parsed as an individual :py:class:`psm_utils.psm.PeptideSpectrumMatch` object.

"""


from __future__ import annotations

import re

from pyteomics.openms import idxml

from psm_utils.exceptions import PSMUtilsException
from psm_utils.io._base_classes import ReaderBase
from psm_utils.psm import PeptideSpectrumMatch
from psm_utils.psm_list import PSMList

# patterns to match open and closed round brackets
MODIFICATION_PATTERN = re.compile(r"\(((?:[^)(]+|\((?:[^)(]+|\([^)(]*\))*\))*)\)")
MODIFICATION_PATTERN_NTERM = re.compile(
    r"^\.\(((?:[^)(]+|\((?:[^)(]+|\([^)(]*\))*\))*)\)"
)


class IdXMLReader(ReaderBase):
    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""
        with idxml.read(str(self.filename)) as reader:
            for entry in reader:
                for peptide_hit in entry["PeptideHit"]:
                    yield self._parse_psm(entry, peptide_hit)

    def read_file(self) -> PSMList:
        """Read full PSM file into a PSMList object."""
        return PSMList([psm for psm in self.__iter__()])

    @staticmethod
    def _parse_peptidoform(sequence: str, charge: int):
        """Parse idXML peptide to :py:class:`psm_utils.peptidoform.Peptidoform`."""
        # N-terminal modification
        if sequence[0] == ".":
            sequence = MODIFICATION_PATTERN_NTERM.sub(r"[\1]-", sequence)
        # Sequence modifications
        sequence = MODIFICATION_PATTERN.sub(r"[\1]", sequence)

        # TODO: C-terminal modifications?

        # Add charge
        sequence += f"/{charge}"
        return sequence

    @staticmethod
    def _parse_is_decoy(target_decoy: str):
        if target_decoy == "target":
            return False
        elif target_decoy == "decoy":
            return True
        else:
            return None

    def _parse_psm(self, entry: dict, peptide_hit: dict) -> PeptideSpectrumMatch:
        """Parse idXML PSM to :py:class:`psm_utils.psm.PeptideSpectrumMatch`."""
        return PeptideSpectrumMatch(
            peptide=self._parse_peptidoform(
                peptide_hit["sequence"], peptide_hit["charge"]
            ),
            spectrum_id=entry["spectrum_reference"],
            is_decoy=self._parse_is_decoy(peptide_hit["target_decoy"]),
            score=peptide_hit["score"],
            precursor_mz=entry["MZ"],
            retention_time=entry["RT"],
            protein_list=[protein["accession"] for protein in peptide_hit["protein"]],
            source="idXML",
            metadata={
                "idxml:score_type": entry["score_type"],
                "idxml:higher_score_better": entry["higher_score_better"],
                "idxml:significance_threshold": entry["significance_threshold"],
                "idxml:consensus_support": peptide_hit["consensus_support"],
            },
        )


class XTandemException(PSMUtilsException):
    pass


class XTandemModificationException(XTandemException):
    pass
