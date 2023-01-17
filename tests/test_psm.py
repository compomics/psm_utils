from psm_utils.peptidoform import Peptidoform


class TestPeptidoform:
    def test_parse_simple(self):
        pep = Peptidoform("ACDE")
        assert pep.sequence == "ACDE"

    def test_is_modified(self):
        for pep, expected in [
            ("ACDE", False),
            ("AC[U:4]DE", True),
            ("[iTRAQ4plex]-AC[U:4]DE", True),
            ("[iTRAQ4plex]-ACDE", True),
            ("ACDE-[Amidation]", True),
            ("{Glycan:Hex}ACDE", True),
            ("<[U:4]@C>ACDE", True),
        ]:
            assert Peptidoform(pep).is_modified == expected
