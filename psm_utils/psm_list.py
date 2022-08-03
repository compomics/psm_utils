from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Union

from psm_utils.psm import PeptideSpectrumMatch


@dataclass
class PSMList:
    """
    Data class representing a list of PSMs.

    Parameters
    ----------
    psm_list : list[PeptideSpectrumMatch]
        List of PeptideSpectrumMatch instances.

    Examples
    --------
    Initiate a :py:class:`PSMList` from a list of PeptideSpectrumMatch objects:

    >>> psm_list = PSMList([
    ...     PeptideSpectrumMatch(peptide="ACDK", spectrum_id=1),
    ...     PeptideSpectrumMatch(peptide="CDEFR", spectrum_id=2),
    ...     PeptideSpectrumMatch(peptide="DEM[Oxidation]K", spectrum_id=3),
    ... ])

    :py:class:`PSMList` directly supports iteration:

    >>> for psm in psm_list:
    ...     print(psm.peptide.theoretical_mass)
    436.12639936491996
    512.1576994932
    454.15222018994

    :py:class:`PSMList` supports indexing and slicing:

    >>> psm_list[1].peptide.theoretical_mass
    512.1576994932


    """

    # TODO
    # __slots__ = ()

    psm_list: list[PeptideSpectrumMatch]

    def __iter__(self) -> Iterable[PeptideSpectrumMatch]:
        return self.psm_list.__iter__()

    def __getitem__(
        self, idx
    ) -> Union[PeptideSpectrumMatch, list[PeptideSpectrumMatch]]:
        # TODO: Expand usage? E.g. index by spectrum_id? Return new PSMList for slice?
        return self.psm_list[idx]

    def __len__(self) -> int:
        return self.psm_list.__len__()

    def rename_modifications(self, mapping: dict[str, str]) -> None:
        """
        Apply mapping to rename modification tags for all PSMs.

        Applies :py:meth:`psm_utils.peptidoform.Peptidoform.rename_modifications` on
        all PSM peptidoforms in the :py:class:`PSMList`.

        Parameters
        ----------
        mapping : dict[str, str]
            Mapping of ``old label`` → ``new label`` for each modification that
            requires renaming. Modification labels that are not in the mapping will not
            be renamed.

        """
        for psm in self.psm_list:
            psm.peptide.rename_modifications(mapping)

    @classmethod
    def from_csv(cls) -> "PSMList":
        """Read PSMList from comma-separated values file."""
        raise NotImplementedError

    def to_csv(self) -> None:  # Or "to_dataframe"?
        """Write PSMList to comma-separated values file."""
        raise NotImplementedError
