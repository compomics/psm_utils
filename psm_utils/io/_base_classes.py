"""Abstract base classes for psm_utils.io."""

from abc import ABC, abstractmethod
from pathlib import Path

from psm_utils.psm import PeptideSpectrumMatch
from psm_utils.psm_list import PSMList


class ReaderBase(ABC):
    """Abstract base class for PSM file readers."""

    def __init__(self, filename) -> None:
        super().__init__()

        self.filename = Path(filename)
        self._validate_filepath()

    def __iter__(self):
        return self

    @abstractmethod
    def __next__(self) -> PeptideSpectrumMatch:
        raise NotImplementedError()

    @abstractmethod
    def read_file() -> PSMList:
        """Read full PSM file into a PSMList object."""
        raise NotImplementedError()

    def _validate_filepath(self):
        """Check if given filepath exists"""

        if not self.filename.exists():
            raise FileNotFoundError("Please give a valid filepath")


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
