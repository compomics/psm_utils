from psm_utils.peptidoform import Peptidoform


class TestPeptidoform:
    def test_rename_modifications(self):
        label_mapping = {"ac": "Acetyl", "cm": "Carbamidomethyl"}

        test_cases = [
            ("ACDEFGHIK", "ACDEFGHIK"),
            ("[ac]-AC[cm]DEFGHIK", "[Acetyl]-AC[Carbamidomethyl]DEFGHIK"),
            ("[ac]-AC[Carbamidomethyl]DEFGHIK", "[Acetyl]-AC[Carbamidomethyl]DEFGHIK"),
            ("[Acetyl]-AC[cm]DEFGHIK", "[Acetyl]-AC[Carbamidomethyl]DEFGHIK"),
            ("<[cm]@C>[Acetyl]-ACDEFGHIK", "<[Carbamidomethyl]@C>[Acetyl]-ACDEFGHIK"),
            ("<[Carbamidomethyl]@C>[ac]-ACDEFGHIK", "<[Carbamidomethyl]@C>[Acetyl]-ACDEFGHIK"),
        ]

        for test_case_in, expected_out in test_cases:
            peptidoform = Peptidoform(test_case_in)
            peptidoform.rename_modifications(label_mapping)
            assert peptidoform.proforma == expected_out
