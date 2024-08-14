"""Tests for psm_utils.io.tsv."""

import pytest

from psm_utils.io.exceptions import PSMUtilsIOException  # noqa: F401
from psm_utils.io.tsv import TSVReader, TSVWriter  # noqa: F401

test_cases = [
    (
        {"peptidoform": "ACDE", "spectrum_id": "1"},
        {
            "peptidoform": "ACDE",
            "spectrum_id": "1",
            "provenance_data": {},
            "metadata": {},
            "rescoring_features": {},
        },
    ),
    (
        {"peptidoform": "ACDE", "spectrum_id": "1", "provenance:test": "value"},
        {
            "peptidoform": "ACDE",
            "spectrum_id": "1",
            "provenance_data": {"test": "value"},
            "metadata": {},
            "rescoring_features": {},
        },
    ),
]


class TestTSVReader:
    def test__parse_entry(self):
        for test_in, expected_out in test_cases:
            assert TSVReader._parse_entry(test_in) == expected_out

    def test_iter(self):
        reader = TSVReader("tests/test_data/test.tsv")
        for psm in reader:
            assert psm.peptidoform == "ACDEK/2"
            assert psm.spectrum_id == "peptide1"
            assert psm.provenance_data == {}
            assert psm.metadata == {}
            assert psm.rescoring_features == {}
            break

    def test_iter_raises(self):
        with TSVReader("tests/test_data/peprec.tsv") as reader:
            with pytest.raises(PSMUtilsIOException):
                for psm in reader:
                    pass
