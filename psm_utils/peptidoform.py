from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from pyteomics import mass, proforma

from psm_utils.exceptions import PSMUtilsException


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
                        f"Cannot resolve composition for modification {tag.value}."
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
                            f"{tag.value}."
                        )
            comp_list.append(position_comp)

        # C-terminus
        c_term = mass.Composition({"H": 1, "O": 1})
        if self.properties["c_term"]:
            for tag in self.properties["c_term"]:
                try:
                    c_term += tag.composition
                except (AttributeError, KeyError):
                    raise ModificationException(
                        f"Cannot resolve composition for modification {tag.value}."
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
                    f"Cannot resolve composition for modification {tag.value}."
                )
        for tag in self.properties["unlocalized_modifications"]:
            try:
                comp += tag.composition
            except (AttributeError, KeyError):
                raise ModificationException(
                    f"Cannot resolve composition for modification {tag.value}."
                )
        return comp

    @property
    def sequential_theoretical_mass(self) -> float:
        """Monoisotopic mass of both termini and each (modified) residue."""
        fixed_rules = {}
        for rule in self.properties["fixed_modifications"]:
            for aa in rule.targets:
                fixed_rules[aa] = rule.modification_tag.mass

        mass_list = []

        # N-terminus
        n_term = mass.Composition({"H": 1}).mass()
        if self.properties["n_term"]:
            for tag in self.properties["n_term"]:
                try:
                    n_term += tag.mass
                except (AttributeError, KeyError):
                    raise ModificationException(
                        f"Cannot resolve mass for modification {tag.value}."
                    )
        mass_list.append(n_term)

        # Sequence
        for aa, tags in self.parsed_sequence:
            # Amino acid
            try:
                position_mass = mass.std_aa_mass[aa]
            except (AttributeError, KeyError):
                raise AmbiguousResidueException(
                    f"Cannot resolve mass for amino acid {aa}."
                )
            # Fixed modifications
            if aa in fixed_rules:
                position_mass += fixed_rules[aa]
            # Localized modifications
            if tags:
                for tag in tags:
                    try:
                        position_mass += tag.mass
                    except (AttributeError, KeyError):
                        raise ModificationException(
                            "Cannot resolve mass for modification " f"{tag.value}."
                        )
            mass_list.append(position_mass)

        # C-terminus
        c_term = mass.Composition({"H": 1, "O": 1}).mass()
        if self.properties["c_term"]:
            for tag in self.properties["c_term"]:
                try:
                    c_term += tag.mass
                except (AttributeError, KeyError):
                    raise ModificationException(
                        f"Cannot resolve mass for modification {tag.value}."
                    )
        mass_list.append(c_term)

        return mass_list

    @property
    def theoretical_mass(self) -> float:
        """Monoisotopic mass of the full uncharged (modified) peptide."""
        mass = sum(self.sequential_theoretical_mass)
        for tag in self.properties["labile_modifications"]:
            try:
                mass += tag.mass
            except (AttributeError, KeyError):
                raise ModificationException(
                    f"Cannot resolve mass for modification {tag.value}."
                )
        for tag in self.properties["unlocalized_modifications"]:
            try:
                mass += tag.mass
            except (AttributeError, KeyError):
                raise ModificationException(
                    f"Cannot resolve mass for modification {tag.value}."
                )
        return mass

    def rename_modifications(self, mapping: dict[str, str]) -> None:
        """
        Apply mapping to rename modification tags.

        Parameters
        ----------
        mapping : dict[str, str]
            Mapping of ``old label`` → ``new label`` for each modification that
            requires renaming. Modification labels that are not in the mapping will not
            be renamed.

        Examples
        --------
        >>> peptide = Peptidoform('[Acedyl]-AC[Carbamidomedyl]DEFGHIK')
        >>> peptide.rename_modifications({
        ...     "Acedyl": "Acetyl",
        ...     "Carbamidomedyl": "Carbamidomethyl"
        ... })
        >>> peptide.proforma
        '[Acetyl]-AC[Carbamidomethyl]DEFGHIK'

        """

        def _rename_modification_list(mods):
            new_mods = []
            for mod in mods:
                if mod.value in mapping:
                    new_mods.append(proforma.process_tag_tokens(mapping[mod.value]))
                else:
                    new_mods.append(mod)
            return new_mods

        # Sequential modifications
        for i, (aa, mods) in enumerate(self.parsed_sequence):
            if mods:
                new_mods = _rename_modification_list(mods)
                self.parsed_sequence[i] = (aa, new_mods)

        # Non-sequence modifications
        for mod_type in [
            "n_term",
            "c_term",
            "unlocalized_modifications",
            "labile_modifications",
            "fixed_modifications",
        ]:
            if self.properties[mod_type]:
                self.properties[mod_type] = _rename_modification_list(
                    self.properties[mod_type]
                )

        # Rebuild proforma string with new labels
        self.proforma = proforma.to_proforma(self.parsed_sequence, **self.properties)


class PeptidoformException(PSMUtilsException):
    """Error while handling :py:class:`Peptidoform`."""

    pass


class AmbiguousResidueException(PeptidoformException):
    """Error while handling ambiguous residue."""

    pass


class ModificationException(PeptidoformException):
    """Error while handling amino acid modification."""

    pass
