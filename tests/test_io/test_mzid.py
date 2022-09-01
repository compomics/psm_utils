import xml.etree.ElementTree as ET

import psm_utils.io.mzid as mzid
from psm_utils import peptidoform, psm, psm_list

PEAKS_TEST_FILE = "./tests/test_data/test_peaks.mzid"
MSGF_TEST_FILE = "./tests/test_data/test_msgf.mzid"


class TestMzIdentMlReader:
    def test_infer_source(self):

        msgf = mzid.MzidReader(MSGF_TEST_FILE)
        assert msgf.source == "MS-GF+"

        peaks = mzid.MzidReader(PEAKS_TEST_FILE)
        assert peaks.source == "PEAKS Studio"

    def test_get_namespace(self):

        root = ET.parse(MSGF_TEST_FILE).getroot()
        root_tag = root.tag
        assert (
            mzid.MzidReader._get_xml_namespace(root_tag)
            == "{http://psidev.info/psi/pi/mzIdentML/1.1}"
        )

    def test_parse_peptidoform(self):

        test_cases = {
            "TMKQNAVS": [
                {
                    "monoisotopicMassDelta": 15.994915,
                    "location": 2,
                    "residues": ["M"],
                    "name": "Oxidation",
                }
            ],
            "KEMVGRIRYGKHRSPMRKQEKT": [
                {
                    "monoisotopicMassDelta": 15.994915,
                    "location": 3,
                    "residues": ["M"],
                    "name": "Oxidation",
                },
                {
                    "monoisotopicMassDelta": 15.994915,
                    "location": 16,
                    "residues": ["M"],
                    "name": "Oxidation",
                },
            ],
            "NGFLMRKV": [
                {
                    "monoisotopicMassDelta": 42.010565,
                    "location": 0,
                    "name": "Acetyl",
                }
            ],
            "GETKMVETAL": [
                {
                    "monoisotopicMassDelta": 15.994915,
                    "location": 5,
                    "residues": ["M"],
                    "name": "Oxidation",
                },
                {
                    "monoisotopicMassDelta": 42.010565,
                    "location": 0,
                    "name": "Acetyl",
                },
            ],
            "RESRSGQSSGY": [
                {
                    "monoisotopicMassDelta": -0.984016,
                    "location": 12,
                    "name": "Amidated",
                }
            ],
        }
        expected_output = {
            "TMKQNAVS": "TM[Oxidation]KQNAVS",
            "KEMVGRIRYGKHRSPMRKQEKT": "KEM[Oxidation]VGRIRYGKHRSPM[Oxidation]RKQEKT",
            "NGFLMRKV": "[Acetyl]-NGFLMRKV",
            "GETKMVETAL": "[Acetyl]-GETKM[Oxidation]VETAL",
            "RESRSGQSSGY": "RESRSGQSSGY-[Amidated]",
        }

        for seq in test_cases.keys():

            assert (
                mzid.MzidReader._parse_peptidoform(seq, test_cases[seq]).proforma
                == expected_output[seq]
            )

    def test_parse_peptide_evidence_ref(self):

        test_cases = [
            [
                {
                    "isDecoy": True,
                    "end": 3242,
                    "start": 3235,
                    "pre": "Q",
                    "post": "F",
                    "PeptideSequence": "TMKQNAVS",
                    "Modification": [
                        {
                            "monoisotopicMassDelta": 15.994915,
                            "location": 2,
                            "residues": ["M"],
                            "name": "Oxidation",
                        }
                    ],
                    "accession": "#DECOY#Q86UQ4|ABCAD_HUMAN",
                    "length": 5058,
                    "Seq": "prot_seq",
                    "protein description": "protein_descr",
                    "numDatabaseSequences": 20844,
                    "location": "fast_location",
                    "name": "human_2022_02",
                    "FileFormat": "FASTA format",
                    "DatabaseName": {"human_2022_02": ""},
                }
            ],
            [
                {
                    "isDecoy": False,
                    "end": 1063,
                    "start": 1053,
                    "pre": "H",
                    "post": "G",
                    "PeptideSequence": "RESRSGQSSGY",
                    "accession": "#CONTAM#Q86YZ3",
                    "length": 2850,
                    "Seq": "prot_seq",
                    "protein description": "protein_descr",
                    "numDatabaseSequences": 20844,
                    "location": "fasta_location",
                    "name": "human_2022_02",
                    "FileFormat": "FASTA format",
                    "DatabaseName": {"human_2022_02": ""},
                },
                {
                    "isDecoy": False,
                    "end": 1063,
                    "start": 1053,
                    "pre": "H",
                    "post": "G",
                    "PeptideSequence": "RESRSGQSSGY",
                    "accession": "Q86YZ3|HORN_HUMAN",
                    "length": 2850,
                    "Seq": "prot_seq",
                    "protein description": "protein_descr",
                    "numDatabaseSequences": 20844,
                    "location": "fasta_location",
                    "name": "human_2022_02",
                    "FileFormat": "FASTA format",
                    "DatabaseName": {"human_2022_02": ""},
                },
            ],
        ]

        expected_output = [
            (True, ["#DECOY#Q86UQ4|ABCAD_HUMAN"]),
            (False, ["#CONTAM#Q86YZ3", "Q86YZ3|HORN_HUMAN"]),
        ]
        for i, test_case in enumerate(test_cases):
            assert (
                mzid.MzidReader._parse_peptide_evidence_ref(test_case)
                == expected_output[i]
            )

    def test_get_rawfile_name(self):

        test_cases = {
            "file:/F:/FromE/H088333_1p_WNE15_trap24_CMB-867_3_LFQ_PNE1904-039.mgf": "H088333_1p_WNE15_trap24_CMB-867_3_LFQ_PNE1904-039",
            "C:/Users/robbin/Desktop/mgf/ignore/YPIC_Run1_CH2-1.mgf": "YPIC_Run1_CH2-1",
            "/Users/kims336/Research/Data/QCShew/test.mgf": "test",
        }

        for test_case in test_cases.keys():
            assert mzid.MzidReader._get_rawfile_name(test_case) == test_cases[test_case]

    # TODO
    def test_get_peptide_spectrum_match(self):
        pass

    # TODO 2
    def test_read_file(self):
        pass
