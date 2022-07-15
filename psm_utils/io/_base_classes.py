"""Abstract base classes for psm_utils.io."""

from __future__ import annotations

import warnings
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Union

from pyteomics import proforma

from psm_utils.io.exceptions import (
    InvalidModificationDefinitionError,
    InvalidModificationError,
    UnresolvableModificationError,
)
from psm_utils.psm import PeptideSpectrumMatch
from psm_utils.psm_list import PSMList


class ReaderBase(ABC):
    """Abstract base class for PSM file readers."""

    def __init__(self,
        filename: Union[str, Path],
        modification_definitions: list[dict[str, str]] = None,
    ) -> None:
        """
        Reader for PSM file.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to PSM file.
        modification_definitions: list[dict[str, str]]
            List of modification definition dictionaries. See
            `Parsing of modification labels`_ for more information.

        """
        super().__init__()

        self.filename = Path(filename)

        if modification_definitions:
            self.modification_definitions = modification_definitions
        else:
            self.modification_definitions = {}
        self.validate_modification_definitions(self.modification_definitions)

    def __iter__(self):
        return self

    @abstractmethod
    def __next__(self) -> PeptideSpectrumMatch:
        raise NotImplementedError()

    @staticmethod
    def validate_modification_definitions(
        modification_definitions: list[dict[str, str]]
    ) -> None:
        """
        Validate resolvability of ProForma labels in modification definitions.

        Validate that ``proforma_label`` values in ``modification_definitions`` list can be
        resolved and that it has working `.composition` and `.mass` attributes.

        Parameters
        ----------
        modification_definitions: list[dict[str, str]]
            List of modification definition dictionaries. See
            `Parsing of modification labels`_ for more information.

        Raises
        ------
        :py:class:`InvalidModificationException`
            If the ProForma label could not be parsed
        :py:class:`UnresolvableModificationException`
            If the ProForma label could not be resolved against the supported
            controlled vocabularies.

        Warns
        -----
        `Atomic composition for modification could not be resolved.`
            If the mass of the modification could be resolved, but not its atomic
            composition.
        """
        for mod_def in modification_definitions:
            # Validate required keys
            for key in ["site", "search_engine_label", "proforma_label"]:
                if key not in mod_def:
                    raise InvalidModificationDefinitionError(
                        f"`{key}` missing from modification definition."
                    )

            # Validate that proforma_label in modification definitions resolvable
            proforma_label = mod_def["proforma_label"]

            # Parse modification
            try:
                proforma_modification = proforma.process_tag_tokens(proforma_label)
            except ValueError:
                raise InvalidModificationError(
                    f"Could not process modification ProForma label `{proforma_label}`"
                )

            # Validate that Tag is in fact a Modification, not another supported ProForma tag
            if not isinstance(
                proforma_modification,
                (proforma.MassModification, proforma.ModificationBase),
            ):
                raise InvalidModificationError(
                    f"ProForma label `{proforma_label}` was not recognized as a valid ProForma modification."
                )

            # "Duck type" whether composition and mass attributes are available; avoid
            # hardcoding Modification types to support other (future) ModificationBase objects
            try:
                _ = proforma_modification.composition
            except (KeyError, AttributeError):
                try:
                    _ = proforma_modification.mass
                except (KeyError, AttributeError):
                    raise UnresolvableModificationError(
                        f"Atomic composition nor mass shift could be resolved for modification "
                        f"with ProForma label `{proforma_label}`."
                    )
                else:
                    warnings.warn(
                        f"Atomic composition for modification `{proforma_label}` could "
                        "not be resolved.",
                        UserWarning,
                    )

    @abstractmethod
    def read_file() -> PSMList:
        """Read full PSM file into a PSMList object."""
        raise NotImplementedError()


class WriterBase(ABC):
    """Abstract base class for PSM file writers."""

    def __init__(self, filename):
        self.filename = filename

    @abstractmethod
    def write_psm(self, psm: PeptideSpectrumMatch):
        """Write a single PSM to the PSM file."""

    @abstractmethod
    def write_file(self, psm_list: PSMList):
        """Write an entire PSMList to the PSM file."""
        raise NotImplementedError()
