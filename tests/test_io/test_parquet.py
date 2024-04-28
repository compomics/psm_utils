"""Tests for psm_utils.io.tsv."""

import os
import hashlib

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



def compute_checksum(filename):
    hash_func = hashlib.sha256()
    with open(filename, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b''):
            hash_func.update(chunk)
    return hash_func.hexdigest()

class TestParquetWriter:
    expected_checksum = "c0782793f8c6fd52e39d5ec1cf5567fb0a7e7e245d795f4f1f720337f756b44c"
    def test_write_psm(self):
        with ParquetWriter("test.pq") as writer:
            for test_case in test_cases:
                writer.write_psm(PSM(**test_case))
        actual_checksum = compute_checksum("test.pq")
        assert actual_checksum == self.expected_checksum, "Checksums do not match"
        os.remove("test.pq")

    def test_write_file(self):
        with ParquetWriter("test.pq") as writer:
            writer.write_file(PSMList(psm_list=[PSM(**t) for t in test_cases]))
        actual_checksum = compute_checksum("test.pq")
        assert actual_checksum == self.expected_checksum, "Checksums do not match"
        # os.remove("test.pq")

class TestParquetReader:
    def test_iter(self):
        # Write test cases to file
        ParquetWriter("test.pq").write_file(PSMList(psm_list=[PSM(**t) for t in test_cases]))

        # Read test cases from file
        for i, psm in enumerate(ParquetReader("test.pq")):
            assert psm == PSM(**test_cases[i])

        os.remove("test.pq")
