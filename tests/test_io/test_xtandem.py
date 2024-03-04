"""Tests for psm_utils.io.xtandem."""

from psm_utils.io.xtandem import XTandemReader


class TestXTandemReader:
    reader = XTandemReader("path")

    def test__parse_peptidoform(self):
        test_cases = [
            {
                "test_in": ({"start": 556, "seq": "KMDYPPKR"}, 2),
                "expected_out": "KMDYPPKR/2",
            },
            {
                "test_in": (
                    {
                        "start": 556,
                        "seq": "KMDYPPKR",
                        "aa": [{"type": "M", "at": 557, "modified": 15.994}],
                    },
                    3,
                ),
                "expected_out": "KM[+15.994]DYPPKR/3",
            },
            {
                "test_in": (
                    {
                        "start": 189,
                        "seq": "CWASLWTAR",
                        "aa": [
                            {"type": "C", "at": 189, "modified": 57.022},
                            {"type": "C", "at": 189, "modified": -17.02655},
                        ],
                    },
                    2,
                ),
                "expected_out": "C[+57.022][-17.02655]WASLWTAR/2",
            },
        ]

        for case in test_cases:
            test_out = self.reader._parse_peptidoform(*case["test_in"]).proforma
            assert test_out == case["expected_out"]
