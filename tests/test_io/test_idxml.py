"""Tests for psm_utils.io.idxml."""

from psm_utils.io.idxml import IdXMLReader
from psm_utils.psm import PSM

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


    def test__parse_psm(self):
        test_psm = PSM(
            peptidoform=Peptidoform('LVPIWKK/2'),
            spectrum_id='controllerType=0 controllerNumber=1 scan=294',
            run='HepG2_rep3_small',
            collection=None,
            spectrum=None,
            is_decoy=True,
            score=1.0099999904632568,
            qvalue=None,
            pep=None,
            precursor_mz=442.29595853189994,
            retention_time=849.0,
            ion_mobility=None,
            protein_list=['DECOY_sp|Q5TH74|STPG1_HUMAN'],
            rank=1,
            source='idXML',
            provenance_data={'LVPIWKK/2': 'LVPIWKK'},
            metadata={
                'idxml:score_type': 'expect',
                'idxml:higher_score_better': 'False',
                'idxml:significance_threshold': '0.0'
            },
            rescoring_features={'MS:1002252': 0.693,
                                'COMET:deltaCn': 1.0,
                                'MS:1002255': 35.9, 'COMET:deltaLCn': 0.0, 'COMET:lnExpect': 0.009950330853168092, 'COMET:lnNumSP': 3.555348061489414, 'COMET:lnRankSP': 0.0, 'COMET:IonFrac': 0.25}
        )

        reader = IdXMLReader("./tests/test_data/test.idXML")
        psm = IdXMLReader._parse_psm(reader.protein_ids, reader.peptide_ids[0], reader.peptide_ids[0])
        assert psm.peptidoform == test_psm.peptidoform
        assert psm.spectrum_id == test_psm.spectrum_id
        assert psm.run == test_psm.run
        assert psm.collection == test_psm.collection
        assert psm.spectrum == test_psm.spectrum
        assert psm.is_decoy == test_psm.is_decoy
        assert psm.score == test_psm.score
        assert psm.qvalue == test_psm.qvalue
        assert psm.pep == test_psm.pep
        assert psm.precursor_mz == test_psm.precursor_mz
        assert psm.retention_time == test_psm.retention_time
        assert psm.ion_mobility == test_psm.ion_mobility
        assert psm.protein_list == test_psm.protein_list
        assert psm.rank == test_psm.rank
        assert psm.source == test_psm.source
        assert psm.provenance_data == test_psm.provenance_data
        assert psm.metadata == test_psm.metadata
        assert psm.rescoring_features == test_psm.rescoring_features


    def test__get_run():
        reader = IdXMLReader("./tests/test_data/test.idXML")
        expected_output = "HepG2_rep3_small"
        assert IdXMLReader._get_run(reader.protein_ids, reader.peptide_ids[0]) == expected_output
    

    def _get_rescoring_features():
        reader = IdXMLReader("./tests/test_data/test.idXML")
        expected_output = ['MS:1002252', 'COMET:deltaCn', 'MS:1002255', 'COMET:deltaLCn', 'COMET:lnExpect', 'COMET:lnNumSP', 'COMET:lnRankSP', 'COMET:IonFrac']
        assert IdXMLReader._get_rescoring_features(reader.peptide_ids[0].getHits()[0]) == expected_output
