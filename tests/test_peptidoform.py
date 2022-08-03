from cProfile import label

from psm_utils.peptidoform import Peptidoform


class TestPeptidoform:
    def test_rename_modifications(self):
        label_mapping = {
            "Acedyl": "Acetyl",
            "Carbamidomedyl": "Carbamidomethyl"
        }

        test_cases = [
            ("ACDEFGHIK", "ACDEFGHIK"),
            ("[Acedyl]-AC[Carbamidomedyl]DEFGHIK", "[Acetyl]-AC[Carbamidomethyl]DEFGHIK"),
            ("[Acedyl]-AC[Carbamidomethyl]DEFGHIK", "[Acetyl]-AC[Carbamidomethyl]DEFGHIK"),
            ("[Acetyl]-AC[Carbamidomedyl]DEFGHIK", "[Acetyl]-AC[Carbamidomethyl]DEFGHIK"),
        ]

        for test_case_in, expected_out in test_cases:
            peptide = Peptidoform(test_case_in)
            peptide.rename_modifications(label_mapping)
            assert peptide.proforma == expected_out
