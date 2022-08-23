"""Tests for psm_utils.io.idxml."""

from psm_utils.io.idxml import IdXMLReader


class TestIdXMLReader:
    def test__parse_peptidoform(self):
        test_cases = [
            (("DFPIAMGER", 2), "DFPIAMGER/2"),
            (("DFPIAM(Oxidation)GER", 2), "DFPIAM[Oxidation]GER/2"),
            ((".DFPIAMGER.", 2), "DFPIAMGER/2"),
            ((".DFPIAM(Oxidation)GER.", 2), "DFPIAM[Oxidation]GER/2"),
            ((".DFPIAM(UniMod:35)GER.", 2), "DFPIAM[UniMod:35]GER/2"),
            ((".DFPIAM[+16]GER.", 2), "DFPIAM[+16]GER/2"),
            ((".DFPIAM[147]GER.", 2), "DFPIAM[147]GER/2"),
            ((".(Dimethyl)DFPIAMGER.", 2), "[Dimethyl]-DFPIAMGER/2"),
            ((".DFPIAMGER.(Label:18O(2))", 2), "DFPIAMGER-[Label:18O(2)]/2"),
            ((".DFPIAMGER(Phospho).", 2), "DFPIAMGER[Phospho]/2"),
        ]

        for test_in, expected_out in test_cases:
            assert IdXMLReader._parse_peptidoform(*test_in) == expected_out
