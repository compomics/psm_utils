"""Tests for psm_utils.io.percolator."""

from psm_utils.io.percolator import PercolatorTabReader


class TestPercolatorTabReader:
    def test__infer_charge_columns(self):
        test_cases = [
            (["Charge"], ("Charge", {})),
            (["charge"], ("charge", {})),
            (["Charge", "Charge2"], ("Charge", {2: "Charge2"})),
            (
                ["Charge2", "Charge3", "Charge4"],
                (None, {2: "Charge2", 3: "Charge3", 4: "Charge4"}),
            ),
            (
                ["charge2", "charge3", "charge4"],
                (None, {2: "charge2", 3: "charge3", 4: "charge4"}),
            ),
        ]

        for test_in, expected_out in test_cases:
            assert expected_out == PercolatorTabReader._infer_charge_columns(test_in)

    def test_parse_peptidoform(self):
        test_cases = [
            (("ACDEFGHR", None), "ACDEFGHR"),
            (("K.ACDEFGHR.I", 1), "ACDEFGHR/1"),
            (("K.ACDEFGHR.-", 2), "ACDEFGHR/2"),
            (("-.ACDEFGHR.I", 3), "ACDEFGHR/3"),
            (("-.ACDEFGHR.-", None), "ACDEFGHR"),
        ]
        for test_in, expected_out in test_cases:
            assert expected_out == PercolatorTabReader._parse_peptidoform(*test_in).proforma
