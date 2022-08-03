import pytest

import psm_utils.io.maxquant as maxquant
from psm_utils import peptidoform, psm, psm_list

TEST_COL = [
    "Raw file",
    "Scan number",
    "Charge",
    "Length",
    "Sequence",
    "Modified sequence",
    "Proteins",
    "Missed cleavages",
    "Mass",
    "Mass error [Da]",
    "Mass error [ppm]",
    "Reverse",
    "Retention time",
    "PEP",
    "Score",
    "Delta score",
    "Localization prob",
    "Matches",
    "Intensities",
    "Mass deviations [Da]",
    "Mass deviations [ppm]",
    "Intensity coverage",
    "id",
]


class TestMaxQuantReader:
    def test_evaluate_columns(self):

        columns = TEST_COL.copy()
        # Test with the right column names
        assert maxquant.MaxQuantReader._evaluate_columns(columns) == True

        # Test with right columns names but lowercase columnname
        columns[0] = "raw file"
        assert maxquant.MaxQuantReader._evaluate_columns(columns) == True

        # Test when column name is missing
        columns.remove("Mass")
        with pytest.raises(maxquant.MSMSParsingError):
            maxquant.MaxQuantReader._evaluate_columns(columns)

    def test_fix_column_case(self):

        expected_rename_dict = {
            "Raw file": "Raw file",
            "Scan number": "Scan number",
            "Sequence": "Sequence",
            "Length": "Length",
            "Missed cleavages": "Missed cleavages",
            "Modified sequence": "Modified sequence",
            "Proteins": "Proteins",
            "Charge": "Charge",
            "Mass": "Mass",
            "Mass error [ppm]": "Mass error [ppm]",
            "Mass error [Da]": "Mass error [Da]",
            "Retention time": "Retention time",
            "PEP": "PEP",
            "Score": "Score",
            "Delta score": "Delta score",
            "Localization prob": "Localization prob",
            "Matches": "Matches",
            "Intensities": "Intensities",
            "Mass Deviations [Da]": "Mass deviations [Da]",
            "Mass Deviations [ppm]": "Mass deviations [ppm]",
            "Intensity coverage": "Intensity coverage",
            "Reverse": "Reverse",
            "id": "id",
        }

        columns = TEST_COL.copy()

        # Test to get rename dict with default msms
        assert maxquant.MaxQuantReader._fix_column_case(columns) == expected_rename_dict

    def test_set_mass_error_unit(self):
        msms_reader = maxquant.MaxQuantReader("./tests/test_data/test_msms.txt")
        # Test dalton mass error case
        assert msms_reader._mass_error_unit == "Da"

        # Test ppm mass error case
        columns = TEST_COL.copy()
        columns.remove("Mass error [Da]")
        msms_reader._set_mass_error_unit(columns)
        assert msms_reader._mass_error_unit == "ppm"

        # Test NotImplementedError on mass error unit
        columns.remove("Mass error [ppm]")
        with pytest.raises(NotImplementedError):
            msms_reader._set_mass_error_unit(columns)

    def test_parse_peptidoform(self):
        test_cases = {
            "input_modified_sequence": [
                "_VGVGFGR_",
                "_MCK_",
                "_(ac)EEEIAALVIDNGSGMCK_",
                "_(gl)QYDADLEQILIQWITTQCRK_",
                "_LAM(ox)QEFMILPVGAANFR_",
                "_VGVN(de)GFGR_",
                "_(ac)EEEIAALVIDNGSGM(ox)CK_",
                "_(ac)SDKPDM(ox)AEIEK_",
                "_YYWGGHYSWDM(Ox)AK_",
                "_YYWGGHYSWDM(Oxidation (M))AK_",
                "_YYWGGHYM(ox)WDM(ox)AK_",
                "_(Acetyl (Protein N-term))ATGPM(ox)SFLK_",
                "_ACDE(Amidated (Peptide C-term))_",
                "_ACM(Ox)DE(Amidated (Peptide C-term))_",
                "_(Acetyl (Protein N-term))M(Ox)ACM(Ox)DEM(Ox)(Amidated (Peptide C-term))_",
            ],
            "expected_output": [
                "VGVGFGR",
                "MCK",
                "[ac]-EEEIAALVIDNGSGMCK",
                "[gl]-QYDADLEQILIQWITTQCRK",
                "LAM[ox]QEFMILPVGAANFR",
                "VGVN[de]GFGR",
                "[ac]-EEEIAALVIDNGSGM[ox]CK",
                "[ac]-SDKPDM[ox]AEIEK",
                "YYWGGHYSWDM[Ox]AK",
                "YYWGGHYSWDM[Oxidation (M)]AK",
                "YYWGGHYM[ox]WDM[ox]AK",
                "[Acetyl (Protein N-term)]-ATGPM[ox]SFLK",
                "ACDE-[Amidated (Peptide C-term)]",
                "ACM[Ox]DE-[Amidated (Peptide C-term)]",
                "[Acetyl (Protein N-term)]-M[Ox]ACM[Ox]DEM[Ox]-[Amidated (Peptide C-term)]",
            ],
        }

        msms_reader = maxquant.MaxQuantReader("./tests/test_data/test_msms.txt")

        for test_in, expected_out in zip(
            test_cases["input_modified_sequence"], test_cases["expected_output"]
        ):
            output = msms_reader._parse_peptidoform(test_in)
            assert output.proforma == expected_out
