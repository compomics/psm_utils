import csv
import logging
import os
import re
from functools import cmp_to_key
from itertools import compress
from pathlib import Path
from typing import Dict, List, Tuple, Union

import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

from psm_utils.exceptions import PSMUtilsException
from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.psm import PeptideSpectrumMatch
from psm_utils.psm_list import PSMList

logger = logging.getLogger(__name__)


class MzidReader(ReaderBase):
    """Reader for MaxQuant msms.txt PSM files."""

    def __init__(self, filename: Union[str, Path]) -> None:
        super().__init__(filename)

        self.filename = filename
        self.source = self._infer_source()

    def __iter__(self):
        raise NotImplementedError()

    def __next__(self) -> PeptideSpectrumMatch:
        raise NotImplementedError()

    def read_file() -> PSMList:
        raise NotImplementedError()

    def _get_peptide_spectrum_match() -> PeptideSpectrumMatch:
        """Return a PeptideSpectrumMatch object from maxquat msms PSM"""
        raise NotImplementedError()

    def _get_peptidoform():
        """Return a peptido form"""
        raise NotImplementedError()

    def _infer_source(self):
        """Get the source of the mzid file"""

        mzid_xml = ET.parse(self.filename)
        root = mzid_xml.getroot()
        name_space = self._get_xml_namespace(root.tag)

        return root.find(f".//{name_space}AnalysisSoftware").attrib["name"]

    @staticmethod
    def _get_xml_namespace(root_tag):
        """Get the namespace of the xml root"""

        m = re.match(r"\{.*\}", root_tag)
        return m.group(0) if m else ""


class MzidWriter(WriterBase):
    """Writer for MaxQuant msms.txt files (not implemented)."""

    def __init__(self, filename):
        raise NotImplementedError()

    def write_psm(self, psm: PeptideSpectrumMatch):
        """Write a single PSM to the MaxQuant msms.txt PSM file."""
        raise NotImplementedError()

    def write_file(self, psm_list: PSMList):
        """Write entire PSMList to the MaxQuant msms.txt PSM file."""
        raise NotImplementedError()
