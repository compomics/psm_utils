import pytest

import psm_utils.io.msamanda as msamanda

TEST_COL = [
    "Title",
    "Sequence",
    "Modifications",
    "Protein Accessions",
    "Amanda Score",
    "m/z",
    "Charge",
    "RT",
    "Filename",
    "number of missed cleavages",
    "number of residues",
    "number of considered fragment ions",
    "delta M",
    "avg MS2 error[ppm]",
    "assigned intensity fraction",
    "binom score",
    "Id",
]


class TestMSMSReader:
    def test_evaluate_columns(self):
        columns = TEST_COL.copy()
        # Test with the right column names
        msamanda.MSAmandaReader("test.file")._evaluate_columns(columns)

        # Test when column name is missing
        columns.remove("m/z")
        with pytest.raises(msamanda.MSAmandaParsingError):
            msamanda.MSAmandaReader("test.file")._evaluate_columns(columns)

    def test_parse_peptidoform(self):
        test_cases = {
            "input": [
                ("ENRFcSWVDQK", "C5(Carbamidomethyl|57.021464|fixed)", 3),
                ("NVSmFPGNPER", "M4(Oxidation|15.994915|variable)", 3),
                (
                    "LRmDTcLQK",
                    "M3(Oxidation|15.994915|variable);C6(Carbamidomethyl|57.021464|fixed)",
                    4,
                ),
                (
                    "LRDTcLQK",
                    "N-Term(Acetyl|40|variable);C5(Carbamidomethyl|57.021464|fixed)",
                    4,
                ),
                ("TLPMFHDEEHAR", "", 3),
                ("VSAGEIAVTGAGR", "C-Term(Amidated|-0.984016|variable)", 2),
                ("VQAELDETK", "", 2),
            ],
            "expected_output": [
                "ENRFC[Carbamidomethyl]SWVDQK/3",
                "NVSM[Oxidation]FPGNPER/3",
                "LRM[Oxidation]DTC[Carbamidomethyl]LQK/4",
                "[Acetyl]-LRDTC[Carbamidomethyl]LQK/4",
                "TLPMFHDEEHAR/3",
                "VSAGEIAVTGAGR-[Amidated]/2",
                "VQAELDETK/2",
            ],
        }

        for test_in, expected_out in zip(test_cases["input"], test_cases["expected_output"]):
            output = msamanda.MSAmandaReader._parse_peptidoform(*test_in).proforma
            assert output == expected_out
