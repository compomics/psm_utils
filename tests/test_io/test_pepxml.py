from psm_utils.peptidoform import Peptidoform
from psm_utils.io.pepxml import PepXMLReader

class TestPepXMLReader:
    def test_parse_peptidoform(self):
        test_cases = [
            {
                "in": {
                    "peptide": "ACDEK",
                    "modifications": [],
                    "charge": 2,
                },
                "out": Peptidoform("ACDEK/2"),
            },
            {
                "in": {
                    "peptide": "STEEQNGGGQK",
                    "modifications": [
                        {"position": 0, "mass": 43.017841151532004},
                        {"position": 2, "mass": 181.014009},
                    ],
                    "charge": 2,
                },
                "out": Peptidoform("[+43.017841151532004]-STE[+181.014009]EQNGGGQK/2"),
            },
            {
                "in": {
                    "peptide": "ACDEK",
                    "modifications": [{"position": 6, "mass": 181.014009}],
                    "charge": 3,
                },
                "out": Peptidoform("ACDEK-[+181.014009]/3"),
            },
        ]

        for test_case in test_cases:
            assert test_case["out"] == PepXMLReader._parse_peptidoform(**test_case["in"])
