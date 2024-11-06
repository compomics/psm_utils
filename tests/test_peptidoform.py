import pytest
from pyteomics import proforma

from psm_utils.peptidoform import Peptidoform, format_number_as_string


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

    def test__eq__(self):
        test_cases = [
            ("ACDEFGHIK", "ACDEFGHIK", True),
            ("ACDEFGHIK", "ACDEFGHI", False),
            ("ACDEFGHIK/2", "ACDEFGHIK/2", True),
            ("ACDEFGHIK/2", "ACDEFGHIK/3", False),
            ("[ac]-AC[cm]DEFGHIK", "[ac]-AC[cm]DEFGHIK", True),
            ("[ac]-AC[cm]DEFGHIK", "[ac]-AC[cm]DEFGH", False),
            ("[ac]-AC[cm]DEFGHIK", "[ac]-AC[cm]DEFGH", False),
            ("[ac]-AC[cm]DEFGHIK", "[ac]-AC[cm]DEFGH", False),
        ]

        for test_case_in_1, test_case_in_2, expected_out in test_cases:
            assert (Peptidoform(test_case_in_1) == test_case_in_2) == expected_out
            assert (Peptidoform(test_case_in_1) == Peptidoform(test_case_in_2)) == expected_out

    with pytest.raises(TypeError):
        Peptidoform("ACDEFGHIK") == 1

    def test__getitem__(self):
        test_cases = [
            ("ACDEFGHIK", 0, ("A", None)),
            ("ACDEFGHIK", 8, ("K", None)),
            ("[ac]-AC[cm]DEFGHIK", 0, ("A", None)),
            ("[ac]-AC[cm]DEFGHIK", 1, ("C", [proforma.GenericModification("cm")])),
            ("[ac]-AC[cm]DEFGHIK", 8, ("K", None)),
        ]

        for test_case_in, index, expected_out in test_cases:
            peptidoform = Peptidoform(test_case_in)
            assert peptidoform[index] == expected_out

    def test__iter__(self):
        for aa, mods in Peptidoform("ACDEM[U:35]K"):
            assert isinstance(aa, str)
            if mods is not None:
                assert isinstance(mods, list)
                for mod in mods:
                    assert isinstance(mod, proforma.TagBase)

    def test_sequence(self):
        test_cases = [
            ("ACDEFGHIK", "ACDEFGHIK"),
            ("[ac]-AC[cm]DEFGHIK", "ACDEFGHIK"),
            ("[ac]-AC[Carbamidomethyl]DEFGHIK", "ACDEFGHIK"),
            ("[Acetyl]-AC[cm]DEFGK", "ACDEFGK"),
            ("<[cm]@C>[Acetyl]-ACDK", "ACDK"),
            ("<[Carbamidomethyl]@C>[ac]-ACDEFGHIK", "ACDEFGHIK"),
        ]

        for test_case_in, expected_out in test_cases:
            peptidoform = Peptidoform(test_case_in)
            assert peptidoform.sequence == expected_out

    def test_modified_sequence(self):
        test_cases = [
            ("ACDEFGHIK", "ACDEFGHIK"),
            ("ACDEFGHIK/3", "ACDEFGHIK"),
            ("[ac]-AC[cm]DEFGHIK", "[ac]-AC[cm]DEFGHIK"),
            ("[ac]-AC[cm]DEFGHIK/3", "[ac]-AC[cm]DEFGHIK"),
            ("<[cm]@C>[Acetyl]-ACDK/3", "<[cm]@C>[Acetyl]-ACDK"),
        ]

        for test_case_in, expected_out in test_cases:
            peptidoform = Peptidoform(test_case_in)
            assert peptidoform.modified_sequence == expected_out

    def test_precursor_charge(self):
        test_cases = [
            ("ACDEFGHIK", None),
            ("ACDEFGHIK/2", 2),
            ("ACDEFGHIK/3", 3),
        ]

        for test_case_in, expected_out in test_cases:
            peptidoform = Peptidoform(test_case_in)
            assert peptidoform.precursor_charge == expected_out

    def test_rename_modifications(self):
        label_mapping = {
            "ac": "Acetyl",
            "cm": "Carbamidomethyl",
            "+57.021": "Carbamidomethyl",
            "-18.010565": "Glu->pyro-Glu",
        }

        test_cases = [
            ("ACDEFGHIK", "ACDEFGHIK"),
            ("[ac]-AC[cm]DEFGHIK", "[Acetyl]-AC[Carbamidomethyl]DEFGHIK"),
            ("[ac]-AC[Carbamidomethyl]DEFGHIK", "[Acetyl]-AC[Carbamidomethyl]DEFGHIK"),
            ("[Acetyl]-AC[cm]DEFGHIK", "[Acetyl]-AC[Carbamidomethyl]DEFGHIK"),
            ("<[cm]@C>[Acetyl]-ACDEFGHIK", "<[Carbamidomethyl]@C>[Acetyl]-ACDEFGHIK"),
            ("<[Carbamidomethyl]@C>[ac]-ACDEFGHIK", "<[Carbamidomethyl]@C>[Acetyl]-ACDEFGHIK"),
            ("[ac]-AC[cm]DEFGHIK", "[Acetyl]-AC[Carbamidomethyl]DEFGHIK"),
            ("AC[+57.021]DEFGHIK", "AC[Carbamidomethyl]DEFGHIK"),
            ("E[-18.010565]DEK", "E[Glu->pyro-Glu]DEK"),
        ]

        for test_case_in, expected_out in test_cases:
            peptidoform = Peptidoform(test_case_in)
            peptidoform.rename_modifications(label_mapping)
            assert peptidoform.proforma == expected_out

    def test_add_apply_fixed_modifications(self):
        test_cases = [
            ("ACDEK", [("Cmm", ["C"])], "AC[Cmm]DEK"),
            ("AC[Cmm]DEK", [("SecondMod", ["C"])], "AC[Cmm][SecondMod]DEK"),
            ("ACDEK", [("TMT6plex", ["K", "N-term"])], "[TMT6plex]-ACDEK[TMT6plex]"),
            ("ACDEK-[CT]", [("TMT6plex", ["K", "N-term"])], "[TMT6plex]-ACDEK[TMT6plex]-[CT]"),
        ]

        for test_case_in, fixed_modifications, expected_out in test_cases:
            peptidoform = Peptidoform(test_case_in)
            peptidoform.add_fixed_modifications(fixed_modifications)
            peptidoform.apply_fixed_modifications()
            assert peptidoform.proforma == expected_out


def test_format_number_as_string():
    test_cases = [
        (1212.12, "+1212.12"),
        (-1212.12, "-1212.12"),
        (0.1, "+0.1"),
        (-0.1, "-0.1"),
        (1212.000, "+1212"),
        (1212.1200, "+1212.12"),
    ]

    for test_case_in, expected_out in test_cases:
        assert format_number_as_string(test_case_in) == expected_out
