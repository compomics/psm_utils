"""Abstract base classes for psm_utils.io."""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Union

from psm_utils.psm import PeptideSpectrumMatch
from psm_utils.psm_list import PSMList


class ReaderBase(ABC):
    """Abstract base class for PSM file readers."""

    def __init__(
        self,
        filename: Union[str, Path],
    ) -> None:
        """
        Reader for PSM file.

        Parameters
        ----------
        filename: str, pathlib.Path
            Path to PSM file.

        """
        super().__init__()

        self.filename = Path(filename)

        # # Create easy-to-use mapping of modification SE label to ProForma label
        # self._label_mapping = {
        #     mod_def["search_engine_label"]: mod_def["proforma_label"]
        #     for mod_def
        #     in self.modification_definitions
        # }

    @abstractmethod
    def __iter__(self):
        raise NotImplementedError()

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
