from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any

from pyteomics import proforma, mass

from psm_utils._exceptions import PSMUtilsException


class PeptidoformException(PSMUtilsException):
    """Error while handling :py:class:`Peptidoform`."""

    pass


class AmbiguousResidueException(PeptidoformException):
    """Error while handling ambiguous residue."""

    pass


class ModificationException(PeptidoformException):
    """Error while handling amino acid modification."""

    pass


@dataclass
class Peptidoform:
    """
    Data Class representing a peptidoform.

    Parameters
    ----------
    proforma : str
        Peptidoform sequence in ProForma v2 notation.

    Attributes
    ----------
    parsed_sequence : list
        List of tuples with residue and modifications for each location.
    properties : dict[str, Any]
        Dict with sequence-wide properties.


    """

    proforma: str
    parsed_sequence: list = field(init=False, repr=False, compare=False)
    properties: dict[str, Any] = field(init=False, repr=False, compare=False)
    sequence: str = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        self.parsed_sequence, self.properties = proforma.parse(self.proforma)
        if self.properties["isotopes"]:
            raise NotImplementedError(
                "Peptidoforms with isotopes are currently not supported."
            )
        self.sequence = "".join(pos[0] for pos in self.parsed_sequence)

    @property
    def ms2pip_modifications(self) -> str:
        """Peptide modifications in MS²PIP PEPREC notation."""

        def _mod_to_ms2pip(mod_site: list, location: int):
            """Proforma modification site (list) to MS²PIP modification."""
            if len(mod_site) > 1:
                raise PeptidoformException(
                    "Multiple modifications per site not supported."
                )
            return "|".join([str(location), mod_site[0].name])

        ms2pip_mods = []
        if self.properties["n_term"]:
            ms2pip_mods.append(_mod_to_ms2pip(self.properties["n_term"], 0))
        for i, (aa, mod) in enumerate(self.parsed_sequence):
            if mod:
                ms2pip_mods.append(_mod_to_ms2pip(mod, i + 1))
        if self.properties["c_term"]:
            ms2pip_mods.append(_mod_to_ms2pip(self.properties["n_term"], -1))

        return "|".join(ms2pip_mods) if ms2pip_mods else "-"

    @property
    def sequential_composition(self) -> list[mass.Composition]:
        """Atomic compositions of both termini and each (modified) residue."""
        # Get compositions for fixed modifications by amino acid
        fixed_rules = {}
        for rule in self.properties["fixed_modifications"]:
            for aa in rule.targets:
                fixed_rules[aa] = rule.modification_tag.composition

        comp_list = []

        # N-terminus
        n_term = mass.Composition({"H": 1})
        if self.properties["n_term"]:
            for tag in self.properties["n_term"]:
                try:
                    n_term += tag.composition
                except (AttributeError, KeyError):
                    raise ModificationException(
                        f"Cannot resolve composition for modification {tag.name}."
                    )
        comp_list.append(n_term)

        # Sequence
        for aa, tags in self.parsed_sequence:
            # Amino acid
            try:
                position_comp = mass.std_aa_comp[aa].copy()
            except (AttributeError, KeyError):
                raise AmbiguousResidueException(
                    f"Cannot resolve composition for amino acid {aa}."
                )
            # Fixed modifications
            if aa in fixed_rules:
                position_comp += fixed_rules[aa]
            # Localized modifications
            if tags:
                for tag in tags:
                    try:
                        position_comp += tag.composition
                    except (AttributeError, KeyError):
                        raise ModificationException(
                            "Cannot resolve composition for modification "
                            f"{tag.name}."
                        )
            comp_list.append(position_comp)

        # C-terminus
        c_term = mass.Composition({"H": 1, "O": 1})
        if self.properties["n_term"]:
            for tag in self.properties["c_term"]:
                try:
                    c_term += tag.composition
                except (AttributeError, KeyError):
                    raise ModificationException(
                        f"Cannot resolve composition for modification {tag.name}."
                    )
        comp_list.append(c_term)

        return comp_list

    @property
    def composition(self) -> mass.Composition:
        """Atomic composition of the full (modified) peptide."""
        comp = mass.Composition()
        for position_comp in self.sequential_composition:
            comp += position_comp
        for tag in self.properties["labile_modifications"]:
            try:
                comp += tag.composition
            except (AttributeError, KeyError):
                raise ModificationException(
                    f"Cannot resolve composition for modification {tag.name}."
                )
        for tag in self.properties["unlocalized_modifications"]:
            try:
                comp += tag.composition
            except (AttributeError, KeyError):
                raise ModificationException(
                    f"Cannot resolve composition for modification {tag.name}."
                )
        return comp

    @property
    def sequential_theoretical_mass(self) -> float:
        """Monoisotopic mass of both termini and each (modified) residue."""
        return [comp.mass() for comp in self.sequential_composition]

    @property
    def theoretical_mass(self) -> float:
        """Monoisotopic mass of the full uncharged (modified) peptide."""
        return sum(self.sequential_theoretical_mass)
