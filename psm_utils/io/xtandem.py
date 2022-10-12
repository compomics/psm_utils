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

* Although X!Tandem XML allows multiple peptide/protein identifications per entry, only
  the first peptide/protein per entry is parsed.

"""


from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np
from pyteomics import mass, tandem

from psm_utils.exceptions import PSMUtilsException
from psm_utils.io._base_classes import ReaderBase
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList


class XTandemReader(ReaderBase):
    def __init__(self, filename: Union[str, Path], decoy_prefix="DECOY_") -> None:
        """
        Reader for X!Tandem XML PSM files.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to PSM file.
        decoy_prefix: str, optional
            Protein name prefix used to denote decoy protein entries. Default:
            ``"DECOY_"``.

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

    def __iter__(self):
        """Iterate over file and return PSMs one-by-one."""
        with tandem.read(str(self.filename)) as reader:
            for entry in reader:
                psm = self._parse_entry(entry)
                yield psm

    def read_file(self) -> PSMList:
        """Read full PSM file into a PSMList object."""
        psm_list = []
        with tandem.read(str(self.filename)) as reader:
            for entry in reader:
                psm_list.append(self._parse_entry(entry))
        return PSMList(psm_list=psm_list)

    def _parse_peptidoform(self, peptide_entry, charge: int) -> Peptidoform:
        """Parse X!Tandem XML peptide entry to :py:class:`~psm_utils.peptidoform.Peptidoform`."""
        if "aa" in peptide_entry:
            # Parse modifications
            seq_list = list(peptide_entry["seq"])
            mod_dict = {}
            for mod_entry in peptide_entry["aa"]:
                # Locations are encoded relative to position in protein
                mod_loc = mod_entry["at"] - peptide_entry["start"]

                # Check if site matches amino acid
                if not mod_entry["type"] == seq_list[mod_loc]:
                    raise XTandemModificationException(
                        f"Found unexpected residue `{seq_list[mod_loc]}` at "
                        f"modification location for `{mod_entry}`."
                    )

                # Add modifications to dict
                if mod_loc not in mod_dict:
                    mod_dict[mod_loc] = float(mod_entry["modified"])
                else:
                    # "sum" multiple modifications per site, e.g.,
                    # cmm + ammonia-loss = pyro-cmm
                    mod_dict[mod_loc] += float(mod_entry["modified"])

            # Add modification in ProForma format
            for mod_loc, mass_shift in mod_dict.items():
                seq_list[mod_loc] += f"[{mass_shift:+g}]"
            proforma_seq = "".join(seq_list)

        else:
            # No modifications to parse
            proforma_seq = peptide_entry["seq"]

        proforma_seq += f"/{charge}"

        return Peptidoform(proforma_seq)

    def _parse_entry(self, entry) -> PSM:
        """Parse X!Tandem XML entry to :py:class:`~psm_utils.psm.PSM`."""
        peptide_entry = entry["protein"][0]["peptide"]
        psm = PSM(
            peptidoform=self._parse_peptidoform(peptide_entry, entry["z"]),
            spectrum_id=entry["support"]["fragment ion mass spectrum"]["note"].split(
                " "
            )[0],
            is_decoy=entry["protein"][0]["label"].startswith(self.decoy_prefix),
            score=-np.log(peptide_entry["expect"]),
            precursor_mz=entry["mh"] - mass.nist_mass["H"][0][0],
            retention_time=entry["rt"],
            protein_list=[entry["protein"][0]["label"]],
            source="X!Tandem",
            provenance_data={
                "xtandem_filename": str(self.filename),
                "xtandem_id": entry["id"],
            },
            metadata={
                "xtandem_hyperscore": peptide_entry["hyperscore"],
                "xtandem_delta": peptide_entry["delta"],
                "xtandem_nextscore": peptide_entry["nextscore"],
            },
        )
        return psm


class XTandemException(PSMUtilsException):
    pass


class XTandemModificationException(XTandemException):
    pass
