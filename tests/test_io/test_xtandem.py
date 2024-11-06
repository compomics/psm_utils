"""Tests for psm_utils.io.xtandem."""

from psm_utils.io.xtandem import XTandemReader
from psm_utils.peptidoform import Peptidoform
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList


class TestXTandemReader:
    def test__parse_peptidoform(self):
        reader = XTandemReader("path")
        test_cases = [
            {
                "test_in": ({"start": 556, "seq": "KMDYPPKR"}, 2),
                "expected_out": "KMDYPPKR/2",
            },
            {
                "test_in": (
                    {
                        "start": 556,
                        "seq": "KMDYPPKR",
                        "aa": [{"type": "M", "at": 557, "modified": 15.994}],
                    },
                    3,
                ),
                "expected_out": "KM[+15.994]DYPPKR/3",
            },
            {
                "test_in": (
                    {
                        "start": 189,
                        "seq": "CWASLWTAR",
                        "aa": [
                            {"type": "C", "at": 189, "modified": 57.022},
                            {"type": "C", "at": 189, "modified": -17.02655},
                        ],
                    },
                    2,
                ),
                "expected_out": "C[+57.022][-17.02655]WASLWTAR/2",
            },
        ]

        for case in test_cases:
            test_out = reader._parse_peptidoform(*case["test_in"]).proforma
            assert test_out == case["expected_out"]

    def test_reader(self):
        test_cases = PSMList(
            psm_list=[
                PSM(
                    peptidoform=Peptidoform('RGFNPDFDVKVDIPER/3'),
                    spectrum_id='635.328491210938_3990.64720000002 RTINSECONDS=3990.64720000002',
                    run='Velos005137',
                    is_decoy=False,
                    score=-4.3694478524670215,
                    precursor_mz=634.656974,
                    retention_time=3990.64720000002,
                    protein_list=[
                        'tr|Q8U2N0|Q8U2N0_PYRFU Uncharacterized protein OS=Pyrococcus furiosus (strain ATCC 43587 / DSM 3638 / JCM 8422 / Vc1) OX=186497 GN=PF0803 PE=4 SV=1'
                    ],
                    source='X!Tandem',
                    provenance_data={
                        'xtandem_filename': 'tests/test_data/test.t.xml',
                        'xtandem_id': '10487'
                    },
                    metadata={
                        'xtandem_hyperscore': '8.8',
                        'xtandem_delta': '0.0048',
                        'xtandem_nextscore': '8.0'
                    },
                ),
                PSM(
                    peptidoform=Peptidoform('GGIGWTGMGIATVKDDGIR/3'),
                    spectrum_id='635.328491210938_3990.64720000002 RTINSECONDS=3990.64720000002',
                    run='Velos005137',
                    is_decoy=True,
                    score=-4.3694478524670215,
                    precursor_mz=634.656974,
                    retention_time=3990.64720000002,
                    protein_list=[
                        'DECOY_tr|Q8TZM9_REVERSED|Q8TZM9_PYRFU Aldose reductase OS=Pyrococcus furiosus (strain ATCC 43587 / DSM 3638 / JCM 8422 / Vc1) OX=186497 GN=PF1960 PE=4 SV=1'
                    ],
                    source='X!Tandem',
                    provenance_data={
                        'xtandem_filename': 'tests/test_data/test.t.xml',
                        'xtandem_id': '10487'
                    },
                    metadata={
                        'xtandem_hyperscore': '8.8',
                        'xtandem_delta': '0.0015',
                        'xtandem_nextscore': '8.0'
                    },
                )
            ]
        )

        psms = XTandemReader("./tests/test_data/test.t.xml").read_file()
        assert psms == test_cases
