"""Tests for psm_utils.io.xtandem."""

from psm_utils.io.xtandem import XTandemReader


class TestXTandemReader:

    reader = XTandemReader("path")

    def test__parse_peptidoform(self):
        test_cases = [
            {
                "test_in": {
                    "start": 556,
                    "seq": "KMDYPPKR",
                },
                "expected_out": "KMDYPPKR",
            },
            {
                "test_in": {
                    "start": 556,
                    "seq": "KMDYPPKR",
                    "aa": [{"type": "M", "at": 557, "modified": 15.994}],
                },
                "expected_out": "KM[+15.994]DYPPKR",
            },
            {
                "test_in": {
                    "start": 189,
                    "seq": "CWASLWTAR",
                    "aa": [
                        {"type": "C", "at": 189, "modified": 57.022},
                        {"type": "C", "at": 189, "modified": -17.02655},
                    ],
                },
                "expected_out": "C[+39.9954]WASLWTAR",
            },
        ]

        for case in test_cases:
            test_out = self.reader._parse_peptidoform(case["test_in"]).proforma
            assert test_out == case["expected_out"]
