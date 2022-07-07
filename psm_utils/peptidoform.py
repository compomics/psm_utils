from dataclasses import dataclass, field
from typing import Any, Optional, Union, List, Dict

import pandas as pd
from pyteomics import proforma

from psm_utils._exceptions import PeptidoformError

@dataclass
class Peptidoform:
    """
    Data Class representing a peptidoform.

    Parameters
    ----------
    sequence : str
        Peptidoform sequence in ProForma v2 notation.

    Attributes
    ----------
    parsed_sequence : list
        List of tuples with residue and modifications for each location.
    modifiers : dict[str, Any]
        Dict with sequence-wide properties.


    """
    sequence: str
    parsed_sequence: list = field(init=False, repr=False, compare=False)
    modifiers: Dict[str, Any] = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        # TODO: Add parsing of custom modifications
        self.parsed_sequence, self.modifiers = proforma.parse(self.sequence)

    @property
    def stripped_sequence(self) -> str:
        """Unmodified peptide sequence."""
        return "".join(position[0] for position in self.parsed_sequence)

    @property
    def ms2pip_modifications(self) -> str:
        """Peptide modifications in MS²PIP PEPREC notation."""

        def _mod_to_ms2pip(mod_site: list, location: int):
            """Proforma modification site (list) to MS²PIP modification."""
            if len(mod_site) > 1:
                raise PeptidoformError(
                    "Multiple modifications per site not supported."
                )
            return "|".join([str(location), mod_site[0].name])

        ms2pip_mods = []
        if self.modifiers["n_term"]:
            ms2pip_mods.append(_mod_to_ms2pip(self.modifiers["n_term"], 0))
        for i, (aa, mod) in enumerate(self.parsed_sequence):
            if mod:
                ms2pip_mods.append(_mod_to_ms2pip(mod, i + 1))
        if self.modifiers["c_term"]:
            ms2pip_mods.append(_mod_to_ms2pip(self.modifiers["n_term"], -1))

        return "|".join(ms2pip_mods) if ms2pip_mods else "-"





# # Copied from spectrum_utils
# def parse(proforma: str) -> List[ProteoForm]:
#     """
#     Parse a ProForma-encoded string.
#     The string is parsed by building an abstract syntax tree to interpret the
#     ProForma language. Next, modifications are localized to residue positions
#     and translated into Python objects.
#     Parameters
#     ----------
#     proforma : str
#         The ProForma string.
#     Returns
#     -------
#     List[Proteoform]
#         A list of proteoforms with their modifications (and additional)
#         information specified by the ProForma string.
#     Raises
#     ------
#     URLError
#         If a controlled vocabulary could not be retrieved from its online
#         resource.
#     KeyError
#         If a term was not found in its controlled vocabulary.
#     NotImplementedError
#         Mass calculation using a molecular formula that contains isotopes (not
#         supported by Pyteomics mass calculation).
#     ValueError
#         If no mass was specified for a GNO term or its parent terms.
#     """
#     dir_name = os.path.dirname(os.path.realpath(__file__))
#     with open(os.path.join(dir_name, "proforma.ebnf")) as f_in:
#         parser = lark.Lark(
#             f_in.read(),
#             start="proforma",
#             parser="earley",
#             lexer="dynamic_complete",
#             import_paths=[dir_name],
#         )
#     # noinspection PyUnresolvedReferences
#     try:
#         return ProFormaTransformer().transform(parser.parse(proforma))
#     except lark.visitors.VisitError as e:
#         raise e.orig_exc
# """