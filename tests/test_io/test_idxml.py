"""Tests for psm_utils.io.idxml."""

import hashlib

from psm_utils.io.idxml import IdXMLReader, IdXMLWriter
from psm_utils.io.sage import SageReader
from psm_utils.psm import PSM
from psm_utils.peptidoform import Peptidoform


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
            peptidoform=Peptidoform("LVPIWKK/2"),
            spectrum_id="controllerType=0 controllerNumber=1 scan=294",
            run="HepG2_rep3_small",
            collection=None,
            spectrum=None,
            is_decoy=True,
            score=1.0099999904632568,
            qvalue=None,
            pep=None,
            precursor_mz=442.29595853189994,
            retention_time=849.0,
            ion_mobility=None,
            protein_list=["DECOY_sp|Q5TH74|STPG1_HUMAN"],
            rank=1,
            source="idXML",
            provenance_data={"LVPIWKK/2": "LVPIWKK"},
            metadata={
                "idxml:score_type": "expect",
                "idxml:higher_score_better": "False",
                "idxml:significance_threshold": "0.0",
                "target_decoy": "decoy",
                "protein_references": "unique",
            },
            rescoring_features={
                "MS:1002252": 0.693,
                "COMET:deltaCn": 1.0,
                "MS:1002255": 35.9,
                "COMET:deltaLCn": 0.0,
                "COMET:lnExpect": 0.009950330853168092,
                "COMET:lnNumSP": 3.555348061489414,
                "COMET:lnRankSP": 0.0,
                "COMET:IonFrac": 0.25,
            },
        )
        reader = IdXMLReader("./tests/test_data/test_in.idXML")
        psm = reader._parse_psm(
            reader.protein_ids, reader.peptide_ids[0], reader.peptide_ids[0].getHits()[0]
        )
        assert psm == test_psm

    def test__get_run(self):
        expected_output = "HepG2_rep3_small"
        reader = IdXMLReader("./tests/test_data/test_in.idXML")
        assert reader._get_run(reader.protein_ids, reader.peptide_ids[0]) == expected_output

    def test__get_rescoring_features(self):
        expected_output = [
            "MS:1002252",
            "COMET:deltaCn",
            "MS:1002255",
            "COMET:deltaLCn",
            "COMET:lnExpect",
            "COMET:lnNumSP",
            "COMET:lnRankSP",
            "COMET:IonFrac",
        ]
        reader = IdXMLReader("./tests/test_data/test_in.idXML")
        assert (
            reader._get_rescoring_features(reader.peptide_ids[0].getHits()[0]) == expected_output
        )


class TestIdXMLWriter:
    def test_write_file_with_pyopenms_objects(self):
        expected_sha = "8d8cb6d8194c5c296f0f5ee8be83d2072be125547b2d51b88100859b001f47fa"
        reader = IdXMLReader("./tests/test_data/test_in.idXML")
        psm_list = reader.read_file()
        writer = IdXMLWriter(
            "./tests/test_data/test_out.idXML", reader.protein_ids, reader.peptide_ids
        )
        writer.write_file(psm_list)
        sha = hashlib.sha256(open("./tests/test_data/test_out.idXML", "rb").read()).hexdigest()
        assert sha == expected_sha

    def test_write_file_without_pyopenms_objects(self):
        expected_sha = "b81addaf8ef1f5cb5007f14a914bee508c54d59f34f8857a5770d3db9aa2c15b"
        reader = SageReader("./tests/test_data/results.sage.tsv")
        psm_list = reader.read_file()
        writer = IdXMLWriter("./tests/test_data/test_out_sage.idXML")
        writer.write_file(psm_list)
        sha = hashlib.sha256(
            open("./tests/test_data/test_out_sage.idXML", "rb").read()
        ).hexdigest()
        assert sha == expected_sha
