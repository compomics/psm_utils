"""Abstract base classes for psm_utils.io."""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Union

from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList


class ReaderBase(ABC):
    """Abstract base class for PSM file readers."""

    def __init__(
        self,
        filename: Union[str, Path],
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

    @abstractmethod
    def __iter__(self):
        raise NotImplementedError()

    @abstractmethod
    def read_file() -> PSMList:
        """Read full PSM file into a PSMList object."""
        raise NotImplementedError()


class WriterBase(ABC):
    """Abstract base class for PSM file writers."""

    def __init__(self, filename, *args, **kwargs):
        super().__init__()
        self.filename = Path(filename)

    @abstractmethod
    def write_psm(self, psm: PSM):
        """Write a single PSM to the PSM file."""

    @abstractmethod
    def write_file(self, psm_list: PSMList):
        """Write an entire PSMList to the PSM file."""
        raise NotImplementedError()
