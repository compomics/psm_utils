"""Tests for psm_utils.io.tsv."""

from psm_utils.io.tsv import TSVReader, TSVWriter

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
