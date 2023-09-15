"""Abstract base classes for psm_utils.io."""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path

from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList


class ReaderBase(ABC):
    """Abstract base class for PSM file readers."""

    def __init__(
        self,
        filename: str | Path,
        *args,
        **kwargs,
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

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        pass

    @abstractmethod
    def __iter__(self):
        raise NotImplementedError()

    def read_file(self) -> PSMList:
        """Read full PSM file into a PSMList object."""
        return PSMList(psm_list=[psm for psm in self.__iter__()])


class WriterBase(ABC):
    """Abstract base class for PSM file writers."""

    def __init__(self, filename, *args, **kwargs):
        super().__init__()
        self.filename = Path(filename)

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        pass

    @abstractmethod
    def write_psm(self, psm: PSM):
        """Write a single PSM to the PSM file."""
        raise NotImplementedError()

    @abstractmethod
    def write_file(self, psm_list: PSMList):
        """Write an entire PSMList to the PSM file."""
        raise NotImplementedError()
