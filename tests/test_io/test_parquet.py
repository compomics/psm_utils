"""Tests for psm_utils.io.tsv."""

import os

from psm_utils.io.parquet import ParquetReader, ParquetWriter
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

test_cases = [
    {"peptidoform": "ACDE", "spectrum_id": "1"},
    {
        "peptidoform": "ACDE",
        "spectrum_id": "2",
        "run": None,
        "collection": None,
        "ion_mobility": None,
        "is_decoy": None,
        "pep": None,
        "precursor_mz": None,
        "protein_list": None,
        "qvalue": None,
        "rank": None,
        "retention_time": None,
        "score": None,
        "source": None,
        "spectrum": None,
        "provenance_data": {"source": "test"},
        "metadata": {},
        "rescoring_features": {"feature": 2.0},
    },
]


class TestParquetWriter:
    def test_write_psm(self):
        with ParquetWriter("test.pq") as writer:
            for test_case in test_cases:
                writer.write_psm(PSM(**test_case))

        with ParquetReader("test.pq") as reader:
            for i, psm in enumerate(reader):
                assert psm == PSM(**test_cases[i])

        os.remove("test.pq")

    def test_write_file(self):
        with ParquetWriter("test.pq") as writer:
            writer.write_file(PSMList(psm_list=[PSM(**t) for t in test_cases]))

        with ParquetReader("test.pq") as reader:
            for i, psm in enumerate(reader):
                assert psm == PSM(**test_cases[i])

        os.remove("test.pq")


class TestParquetReader:
    def test_iter(self):
        # Read test cases from file
        with ParquetReader("tests/test_data/test.pq") as reader:
            for i, psm in enumerate(reader):
                assert psm == PSM(**test_cases[i])
