import xml.etree.ElementTree as ET

import psm_utils.io.mzid as mzid
from psm_utils import peptidoform, psm, psm_list

PEAKS_TEST_FILE = "./tests/test_data/test_peaks.mzid"
MSGF_TEST_FILE = "./tests/test_data/test_msgf.mzid"


class TestMzIdentMlReader:
    def test_infer_source(self):
        msgf = mzid.MzidReader(MSGF_TEST_FILE)
        assert msgf._source == "MS-GF+"

        peaks = mzid.MzidReader(PEAKS_TEST_FILE)
        assert peaks._source == "PEAKS Studio"

    def test_get_namespace(self):
        root = ET.parse(MSGF_TEST_FILE).getroot()
        root_tag = root.tag
        assert (
            mzid.MzidReader._get_xml_namespace(root_tag)
            == "{http://psidev.info/psi/pi/mzIdentML/1.1}"
        )

    def test_parse_peptidoform(self):
        test_cases = {
            "TM[Oxidation]KQNAVS/2": (
                "TMKQNAVS",
                [
                    {
                        "monoisotopicMassDelta": 15.994915,
                        "location": 2,
                        "residues": ["M"],
                        "name": "Oxidation",
                    }
                ],
                2,
            ),
            "KEM[Oxidation]VGRIRYGKHRSPM[Oxidation]RKQEKT/3": (
                "KEMVGRIRYGKHRSPMRKQEKT",
                [
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
                3,
            ),
            "[Acetyl]-NGFLMRKV/4": (
                "NGFLMRKV",
                [
                    {
                        "monoisotopicMassDelta": 42.010565,
                        "location": 0,
                        "name": "Acetyl",
                    }
                ],
                4,
            ),
            "[Acetyl]-GETKM[Oxidation]VETAL/2": (
                "GETKMVETAL",
                [
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
                2,
            ),
            "RESRSGQSSGY-[Amidated]/2": (
                "RESRSGQSSGY",
                [
                    {
                        "monoisotopicMassDelta": -0.984016,
                        "location": 12,
                        "name": "Amidated",
                    }
                ],
                2,
            ),
        }

        for expected_out, test_in in test_cases.items():
            test_out = mzid.MzidReader._parse_peptidoform(*test_in).proforma
            assert test_out == expected_out

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
            assert mzid.MzidReader._parse_peptide_evidence_ref(test_case) == expected_output[i]

    # TODO
    def test_get_peptide_spectrum_match(self):
        pass

    # TODO 2
    def test_read_file(self):
        pass
