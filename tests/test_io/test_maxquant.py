import pytest

import psm_utils.io.maxquant as maxquant
from psm_utils import peptidoform, psm, psm_list

TEST_COL = [
    "Raw file",
    "Scan number",
    "Charge",
    "Modified sequence",
    "Proteins",
    "m/z",
    "Reverse",
    "Retention time",
    "PEP",
    "Score",
]


class TestMSMSReader:
    def test_evaluate_columns(self):
        columns = TEST_COL.copy()
        # Test with the right column names
        maxquant.MSMSReader._evaluate_columns(columns)

        # Test with right columns names but lowercase columnname
        columns[0] = "raw file"
        maxquant.MSMSReader._evaluate_columns(columns)

        # Test when column name is missing
        columns.remove("m/z")
        with pytest.raises(maxquant.MSMSParsingError):
            maxquant.MSMSReader._evaluate_columns(columns)

    def test_parse_peptidoform(self):
        test_cases = {
            "input": [
                ("_VGVGFGR_", 2),
                ("_MCK_", 2),
                ("_(ac)EEEIAALVIDNGSGMCK_", 2),
                ("_(gl)QYDADLEQILIQWITTQCRK_", 2),
                ("_LAM(ox)QEFMILPVGAANFR_", 2),
                ("_VGVN(de)GFGR_", 2),
                ("_(ac)EEEIAALVIDNGSGM(ox)CK_", 2),
                ("_(ac)SDKPDM(ox)AEIEK_", 2),
                ("_YYWGGHYSWDM(Ox)AK_", 2),
                ("_YYWGGHYSWDM(Oxidation (M))AK_", 2),
                ("_YYWGGHYM(ox)WDM(ox)AK_", 2),
                ("_(Acetyl (Protein N-term))ATGPM(ox)SFLK_", 2),
                ("_ACDE(Amidated (Peptide C-term))_", 2),
                ("_ACM(Ox)DE(Amidated (Peptide C-term))_", 2),
                # "_(Acetyl (Protein N-term))M(Ox)ACM(Ox)DEM(Ox)(Amidated (Peptide C-term))_",  # See levitsky/pyteomics/#77
            ],
            "expected_output": [
                "VGVGFGR/2",
                "MCK/2",
                "[ac]-EEEIAALVIDNGSGMCK/2",
                "[gl]-QYDADLEQILIQWITTQCRK/2",
                "LAM[ox]QEFMILPVGAANFR/2",
                "VGVN[de]GFGR/2",
                "[ac]-EEEIAALVIDNGSGM[ox]CK/2",
                "[ac]-SDKPDM[ox]AEIEK/2",
                "YYWGGHYSWDM[Ox]AK/2",
                "YYWGGHYSWDM[Oxidation (M)]AK/2",
                "YYWGGHYM[ox]WDM[ox]AK/2",
                "[Acetyl (Protein N-term)]-ATGPM[ox]SFLK/2",
                "ACDE-[Amidated (Peptide C-term)]/2",
                "ACM[Ox]DE-[Amidated (Peptide C-term)]/2",
                # "[Acetyl (Protein N-term)]-M[Ox]ACM[Ox]DEM[Ox]-[Amidated (Peptide C-term)]",
            ],
        }

        msms_reader = maxquant.MSMSReader("./tests/test_data/test_msms.txt")

        for test_in, expected_out in zip(test_cases["input"], test_cases["expected_output"]):
            output = msms_reader._parse_peptidoform(*test_in)
            assert output.proforma == expected_out
