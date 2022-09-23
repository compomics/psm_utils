from __future__ import annotations

import re
from typing import Iterable, List, Sequence, Union

import numpy as np
import pandas as pd
import pyteomics
from pydantic import BaseModel
from itertools import compress

from psm_utils.psm import PeptideSpectrumMatch


class PSMList(BaseModel):
    """
    Data class representing a list of PSMs.

    Parameters
    ----------
    psm_list : list[PeptideSpectrumMatch]
        List of PeptideSpectrumMatch instances.

    Examples
    --------
    Initiate a :py:class:`PSMList` from a list of PeptideSpectrumMatch objects:

    >>> psm_list = PSMList(psm_list=[
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

    psm_list: List[PeptideSpectrumMatch]

    def __iter__(self) -> Iterable[PeptideSpectrumMatch]:
        return self.psm_list.__iter__()

    def __len__(self) -> int:
        return self.psm_list.__len__()

    def __getitem__(
        self, item
    ) -> Union[PeptideSpectrumMatch, list[PeptideSpectrumMatch]]:
        # TODO: Expand usage? E.g. index by spectrum_id? Return new PSMList for slice?
        if isinstance(item, int):
            return self.psm_list[item]
        elif isinstance(item, slice):
            return PSMList(psm_list=self.psm_list[item])
        elif isinstance(item, str):
            return np.array([psm[item] for psm in self.psm_list])
        elif _is_iterable_of_ints(item):
            # Return new PSMList with selection of PSMs by list indices
            return PSMList(psm_list=[self.psm_list[i] for i in item])
        elif _is_iterable_of_strings(item):
            # Return 2D numpy array of multiple PSM properties along full PSMList
            return np.array([self[inner_item] for inner_item in item]).transpose()
        else:
            raise TypeError(f"Unsupported indexing type: {type(item)}")

    def __setitem__(self, item, values: Sequence) -> None:
        if not len(values) == len(self):
            raise ValueError(f"Expected value with same length as PSMList: {len(self)}")
        for value, psm in zip(values, self):
            psm[item] = value

    @property
    def collections(self) -> list:
        """List of collections in :py:class:`PSMList`."""
        if (self["collection"] != None).any():
            return list(np.unique(self["collection"]))
        else:
            return [None]

    @property
    def runs(self) -> list:
        """List of runs in :py:class:`PSMList`."""
        if (self["run"] != None).any():
            return list(np.unique(self["run"]))
        else:
            return [None]

    def get_psm_dict(self):
        """Get nested dictionary of PSMs by collection, run, and spectrum_id."""
        psm_dict = {}
        for psm in self:
            try:
                psm_dict[psm.collection][psm.run][psm.spectrum_id].append(psm)
            except KeyError:
                try:
                    psm_dict[psm.collection][psm.run][psm.spectrum_id] = [psm]
                except KeyError:
                    try:
                        psm_dict[psm.collection][psm.run] = {psm.spectrum_id: [psm]}
                    except KeyError:
                        psm_dict[psm.collection] = {psm.run: {psm.spectrum_id: [psm]}}
        return psm_dict

    def get_rank1_psms(self, lower_score_better=False) -> PSMList:
        """Return new PSMList with only first-rank PSMs"""
        columns = ["collection", "run", "spectrum_id", "score"]
        df_rank1 = (
            pd.DataFrame(self[columns], columns=columns)
            .sort_values("score", ascending=lower_score_better)
            .drop_duplicates(["collection", "run", "spectrum_id"])
        )
        first_rank_indices = df_rank1.index
        return self[first_rank_indices]

    def find_decoys(self, decoy_pattern: str) -> None:
        """
        Use regular expression pattern to find decoy PSMs by protein name(s).

        This method allows a regular expression pattern to be applied on
        :py:obj:`~psm_utils.psm.PeptideSpectrumMatch`
        :py:attr:`~psm_utils.psm.PeptideSpectrumMatch.protein_list` items to set the
        :py:attr:`~psm_utils.psm.PeptideSpectrumMatch.is_decoy` attribute.
        Decoy protein entries are commonly marked with a prefix or suffix, e.g.
        ``DECOY_``, or ``_REVERSED``. If ``decoy_pattern`` matches to a substring of all
        entries in :py:attr:`protein_list`, the PSM is interpreted as a decoy. Existing
        :py:attr:`is_decoy` entries are overwritten.

        Parameters
        ----------
        decoy_pattern: str
            Regular expression pattern to match decoy protein entries.

        Examples
        --------
        >>> psm_list.find_decoys(r"^DECOY_")

        """
        decoy_pattern = re.compile(decoy_pattern)
        for psm in self:
            psm.is_decoy = all(
                [decoy_pattern.search(p) is not None for p in psm.protein_list]
            )

    def calculate_qvalues(self, reverse: bool = True, **kwargs) -> None:
        """
        Calculate q-values using the target-decoy approach.

        Q-values are calculated for all PSMs from the target and decoy scores. This
        requires that all PSMs have a :py:attr:`score` and a target/decoy state
        (:py:attr:`is_decoy`) assigned. Any existing q-values will be overwritten.

        Parameters
        ----------
        reverse: boolean, optional
            If True (default), a higher score value indicates a better PSM.
        **kwargs : dict, optional
            Additional arguments to be passed to
            `pyteomics.auxiliary.target_decoy.qvalues <https://pyteomics.readthedocs.io/en/latest/api/auxiliary.html#pyteomics.auxiliary.target_decoy.qvalues>`_.

        """
        for key in ["score", "is_decoy"]:
            if (self[key] == None).any():
                raise ValueError(
                    f"Cannot calculate q-values if not all PSMs have `{key}` assigned."
                )

        qvalues = pyteomics.auxiliary.qvalues(
            self,
            key="score",
            is_decoy="is_decoy",
            remove_decoy=False,
            reverse=reverse,
            full_output=True,
            **kwargs,
        )
        for score, is_decoy, qvalue, psm in qvalues:
            psm.qvalue = qvalue

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

    def set_ranks(self):
        "Set the ranks for all the psm in the psm list"

        def rank_simple(vector):
            return sorted(range(len(vector)), key=vector.__getitem__)

        spectrum_list = self.__getitem__("spectrum_id")
        scores = self.__getitem__("score")
        spectrum_ranks = [None] * len(self.psm_list)

        for spectrum_id in set(spectrum_list):
            # Get indices of spectrum ids
            indices = spectrum_list == spectrum_id
            indices = list(compress(range(len(indices)), indices))
            # Get corresponding scores
            spectrum_id_scores = scores[indices]
            # Get ranks for each score
            spectrum_id_ranks = rank_simple(spectrum_id_scores)
            for (index, rank) in zip(indices, spectrum_id_ranks):
                spectrum_ranks[index] = rank

        self.__setitem__("rank", spectrum_ranks)

    @classmethod
    def from_csv(cls) -> "PSMList":
        """Read PSMList from comma-separated values file."""
        raise NotImplementedError

    def to_csv(self) -> None:  # Or "to_dataframe"?
        """Write PSMList to comma-separated values file."""
        raise NotImplementedError

    def to_dataframe(self) -> pd.DataFrame:
        """Convert :py:class:`PSMList` to :py:class:`pandas.DataFrame`."""
        return pd.DataFrame.from_records([psm.__dict__ for psm in self])


def _is_iterable_of_ints(obj):
    try:
        if not all(isinstance(x, int) for x in obj):
            return False
        else:
            return True
    except (TypeError, ValueError):
        return False


def _is_iterable_of_strings(obj):
    try:
        if not all(isinstance(x, str) for x in obj):
            return False
        else:
            return True
    except (TypeError, ValueError):
        return False
