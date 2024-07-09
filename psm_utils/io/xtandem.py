"""
Interface with X!Tandem XML PSM files.


Notes
-----

* In X!Tandem XML, N/C-terminal modifications are encoded as normal modifications and
  are therefore parsed accordingly. Any information on which modifications are
  N/C-terminal is therefore lost.

  N-terminal modification in X!Tandem XML:

  .. code-block::

      <aa type="M" at="1" modified="42.01057" />


* Consecutive modifications, i.e., a modified residue that is modified further, is
  encoded in X!Tandem XML as two distinctive modifications on the same site. However,
  in :py:mod:`psm_utils`, multiple modifications on the same site are not supported.
  While parsing X!Tandem XML PSMs, the mass shift labels of these two modifications will
  therefore be summed into a single modification.

  For example, carbamidomethylation of cystein (57.02200) plus ammonia-loss (-17.02655)
  will be parsed as one modification with mass shift 39.994915, which matches the
  combined modification Pyro-carbamidomethyl:

  .. code-block:: xml

      <aa type="C" at="189" modified="57.02200" />
      <aa type="C" at="189" modified="-17.02655" />


  .. code-block::

      [+39,99545]
"""

from __future__ import annotations

import logging
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Union

import numpy as np
from pyteomics import tandem

from psm_utils.exceptions import PSMUtilsException
from psm_utils.io._base_classes import ReaderBase
from psm_utils.peptidoform import Peptidoform, format_number_as_string
from psm_utils.psm import PSM

logger = logging.getLogger(__name__)


class XTandemReader(ReaderBase):
    def __init__(
        self,
        filename: Union[str, Path],
        *args,
        decoy_prefix="DECOY_",
        score_key="expect",
        **kwargs,
    ) -> None:
        """
        Reader for X!Tandem XML PSM files.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to PSM file.
        decoy_prefix: str, optional
            Protein name prefix used to denote decoy protein entries. Default:
            ``"DECOY_"``.
        score_key: str, optional
            Key of score to use as PSM score. One of ``"expect"``, ``"hyperscore"``,
            ``"delta"``, or ``"nextscore"``. Default: ``"expect"``. The ``"expect"`` score
            (e-value) is converted to its negative natural logarithm to facilitate downstream
            analysis.

        Examples
        --------

        :py:class:`XTandemReader` supports iteration:

        >>> from psm_utils.io.xtandem import XTandemReader
        >>> for psm in XTandemReader("pyro.t.xml"):
        ...     print(psm.peptidoform.proforma)
        WFEELSK
        NDVPLVGGK
        GANLGEMTNAGIPVPPGFC[+57.022]VTAEAYK
        ...

        Or a full file can be read at once into a :py:class:`~psm_utils.psm_list.PSMList`
        object:

        >>> reader = XTandemReader("pyro.t.xml")
        >>> psm_list = reader.read_file()

        """
        super().__init__(filename)
        self.decoy_prefix = decoy_prefix
        self.score_key = score_key

    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""

        with tandem.read(str(self.filename)) as reader:
            run = self._parse_run(self.filename)
            for entry in reader:
                for psm in self._parse_entry(entry, run):
                    yield psm

    @staticmethod
    def _parse_peptidoform(peptide_entry, charge):
        if "aa" in peptide_entry:
            # Parse modifications
            seq_list = list(peptide_entry["seq"])
            unmodified_seq = seq_list.copy()

            for mod_entry in peptide_entry["aa"]:
                # Locations are encoded relative to position in protein
                mod_loc = mod_entry["at"] - peptide_entry["start"]
                mass_shift = float(mod_entry["modified"])

                # Check if site matches amino acid
                if not mod_entry["type"] == unmodified_seq[mod_loc]:
                    raise XTandemModificationException(
                        f"Found unexpected residue `{seq_list[mod_loc]}` at "
                        f"modification location for `{mod_entry}`."
                    )

                # Add to sequence in ProForma format
                seq_list[mod_loc] += f"[{format_number_as_string(mass_shift)}]"

            proforma_seq = "".join(seq_list)

        else:
            # No modifications to parse
            proforma_seq = peptide_entry["seq"]

        proforma_seq += f"/{charge}"

        return Peptidoform(proforma_seq)

    def _parse_entry(self, entry, run: str) -> list:
        """Parse X!Tandem XML entry to a list of :py:class:`~psm_utils.psm.PSM`."""
        pepform_to_psms = dict()

        for protein_entry in entry["protein"]:
            peptide_entry = protein_entry["peptide"]
            peptidoform = self._parse_peptidoform(peptide_entry, entry["z"])

            if peptidoform not in pepform_to_psms:
                psm = PSM(
                    peptidoform=self._parse_peptidoform(peptide_entry, entry["z"]),
                    spectrum_id=entry["support"]["fragment ion mass spectrum"]["note"],
                    is_decoy=protein_entry["label"].startswith(self.decoy_prefix),
                    score=(
                        -np.log(peptide_entry[self.score_key])
                        if self.score_key == "expect"
                        else peptide_entry[self.score_key]
                    ),
                    precursor_mz=entry["mh"] / entry["z"],
                    retention_time=entry["rt"],
                    run=run,
                    protein_list=[protein_entry["note"]],
                    source="X!Tandem",
                    provenance_data={
                        "xtandem_filename": str(self.filename),
                        "xtandem_id": str(entry["id"]),
                    },
                    metadata={
                        "xtandem_hyperscore": str(peptide_entry["hyperscore"]),
                        "xtandem_delta": str(peptide_entry["delta"]),
                        "xtandem_nextscore": str(peptide_entry["nextscore"]),
                    },
                )
                pepform_to_psms[peptidoform] = psm
            else:
                pepform_to_psms[peptidoform].protein_list.append(protein_entry["note"])

        return list(pepform_to_psms.values())

    def _parse_run(self, filepath):
        """Parse X!Tandem XML run to :py:class:`~psm_utils.psm.PSM`."""

        tree = ET.parse(str(filepath))
        root = tree.getroot()
        full_label = root.attrib["label"]
        run_match = re.search(r"\/(?P<run>[^\s\/\\]+)\.(?P<filetype>mgf|mzML|mzml)", full_label)
        if run_match:
            run = run_match.group("run")
        else:
            run = Path(filepath).stem
            logger.warning(
                f"Could not parse run from X!Tandem XML label entry. Setting PSM filename `{run}` "
                "as run."
            )

        return run


class XTandemException(PSMUtilsException):
    pass


class XTandemModificationException(XTandemException):
    pass
