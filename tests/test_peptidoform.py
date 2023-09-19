from pyteomics import proforma

from psm_utils.peptidoform import Peptidoform


class TestPeptidoform:

    def test__len__(self):
        test_cases = [
            ("ACDEFGHIK", 9),
            ("[ac]-AC[cm]DEFGHIK", 9),
            ("[ac]-AC[Carbamidomethyl]DEFGHIK", 9),
            ("[Acetyl]-AC[cm]DEFGK", 7),
            ("<[cm]@C>[Acetyl]-ACDK", 4),
            ("<[Carbamidomethyl]@C>[ac]-ACDEFGHIK", 9),
        ]

        for test_case_in, expected_out in test_cases:
            peptidoform = Peptidoform(test_case_in)
            assert len(peptidoform) == expected_out

    def test__iter__(self):
        for aa, mods in Peptidoform("ACDEM[U:35]K"):
            assert isinstance(aa, str)
            if mods is not None:
                assert isinstance(mods, list)
                for mod in mods:
                    assert isinstance(mod, proforma.TagBase)


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
