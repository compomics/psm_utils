from psm_utils.peptidoform import Peptidoform


class TestPeptideSpectrumMatch:
    def test_parse_simple(self):
        pep = Peptidoform("ACDE")
        assert pep.sequence == "ACDE"
