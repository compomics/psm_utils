from psm_utils.peptidoform import Peptidoform


class TestPeptidoform:
    def test_parse_simple(self):
        pep = Peptidoform("ACDE")
        assert pep.sequence == "ACDE"
