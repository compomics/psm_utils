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

MODIFICATION_DEFINITIONS = [
    {
        "site": "S|T|Y",
        "search_engine_label": "Phospho",
        "proforma_label": "U:21",
    },
    {
        "site": "M",
        "search_engine_label": "Oxidation (M)",
        "proforma_label": "U:35",
    },
    {
        "site": "M",
        "search_engine_label": "ox",
        "proforma_label": "U:35",
    },
    {
        "site": "M",
        "search_engine_label": "Ox",
        "proforma_label": "U:35",
    },
    {
        "site": "N-term",
        "search_engine_label": "Acetyl (Protein N-term)",
        "proforma_label": "U:1",
    },
    {
        "site": "N-term",
        "search_engine_label": "ac",
        "proforma_label": "U:1",
    },
    {
        "site": "N-term",
        "search_engine_label": "gl",
        "proforma_label": "U:28",
    },
    {
        "site": "N",
        "search_engine_label": "de",
        "proforma_label": "U:7",
    },
    {
        "site": "C-term",
        "search_engine_label": "Amidated (Peptide C-term)",
        "proforma_label": "U:2",
    },
    {
        "site": "K",
        "search_engine_label": "Delta:H(4)C(3)",
        "proforma_label": "U:256",
    },
]


class TestMaxquantReader:
    def test_evaluate_columns(self):

        columns = TEST_COL.copy()
        # Test with the right column names
        assert maxquant.MaxquantReader._evaluate_columns(columns) == True

        # Test with right columns names but lowercase columnname
        columns[0] = "raw file"
        assert maxquant.MaxquantReader._evaluate_columns(columns) == True

        # Test when column name is missing
        columns.remove("Mass")
        with pytest.raises(maxquant.MsmsParsingError):
            maxquant.MaxquantReader._evaluate_columns(columns)

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
        assert maxquant.MaxquantReader._fix_column_case(columns) == expected_rename_dict

    def test_set_mass_error_unit(self):
        msms_reader = maxquant.MaxquantReader(
            "./tests/test_data/test_msms.txt", MODIFICATION_DEFINITIONS
        )
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

    def test_parse_maxquant_modification(self):

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
                "_DFK(Delta:H(4)C(3))SK_",
                "_(Acetyl (Protein N-term))ATGPM(ox)SFLK_",
                "_ACDE(Amidated (Peptide C-term))_",
                "_ACM(Ox)DE(Amidated (Peptide C-term))_",
                "_(Acetyl (Protein N-term))M(Ox)ACM(Ox)DEM(Ox)(Amidated (Peptide C-term))_",
            ],
            "expected_output": [
                "VGVGFGR",
                "MCK",
                "[U:1]-EEEIAALVIDNGSGMCK",
                "[U:28]-QYDADLEQILIQWITTQCRK",
                "LAM[U:35]QEFMILPVGAANFR",
                "VGVN[U:7]GFGR",
                "[U:1]-EEEIAALVIDNGSGM[U:35]CK",
                "[U:1]-SDKPDM[U:35]AEIEK",
                "YYWGGHYSWDM[U:35]AK",
                "YYWGGHYSWDM[U:35]AK",
                "YYWGGHYM[U:35]WDM[U:35]AK",
                "DFK[U:256]SK",
                "[U:1]-ATGPM[U:35]SFLK",
                "ACDE-[U:2]",
                "ACM[U:35]DE-[U:2]",
                "[U:1]-M[U:35]ACM[U:35]DEM[U:35]-[U:2]",
            ],
        }

        msms_reader = maxquant.MaxquantReader(
            "./tests/test_data/test_msms.txt", MODIFICATION_DEFINITIONS
        )
        output = [
            msms_reader.parse_maxquant_modification(x)
            for x in test_cases["input_modified_sequence"]
        ]
        assert output == test_cases["expected_output"]

    def test_get_peptidoform(self):
        msms_reader = maxquant.MaxquantReader(
            "./tests/test_data/test_msms.txt", MODIFICATION_DEFINITIONS
        )

        test_cases = {
            "_VGVGFGR_": peptidoform.Peptidoform("VGVGFGR"),
            "_MCK_": peptidoform.Peptidoform("MCK"),
            "_(ac)EEEIAALVIDNGSGMCK_": peptidoform.Peptidoform(
                "[U:1]-EEEIAALVIDNGSGMCK"
            ),
            "_(gl)QYDADLEQILIQWITTQCRK_": peptidoform.Peptidoform(
                "[U:28]-QYDADLEQILIQWITTQCRK"
            ),
            "_LAM(ox)QEFMILPVGAANFR_": peptidoform.Peptidoform(
                "LAM[U:35]QEFMILPVGAANFR"
            ),
            "_VGVN(de)GFGR_": peptidoform.Peptidoform("VGVN[U:7]GFGR"),
            "_(ac)EEEIAALVIDNGSGM(ox)CK_": peptidoform.Peptidoform(
                "[U:1]-EEEIAALVIDNGSGM[U:35]CK"
            ),
            "_(ac)SDKPDM(ox)AEIEK_": peptidoform.Peptidoform("[U:1]-SDKPDM[U:35]AEIEK"),
            "_YYWGGHYSWDM(Ox)AK_": peptidoform.Peptidoform("YYWGGHYSWDM[U:35]AK"),
            "_YYWGGHYSWDM(Oxidation (M))AK_": peptidoform.Peptidoform(
                "YYWGGHYSWDM[U:35]AK"
            ),
            "_YYWGGHYM(ox)WDM(ox)AK_": peptidoform.Peptidoform(
                "YYWGGHYM[U:35]WDM[U:35]AK"
            ),
            "_DFK(Delta:H(4)C(3))SK_": peptidoform.Peptidoform("DFK[U:256]SK"),
            "_(Acetyl (Protein N-term))ATGPM(ox)SFLK_": peptidoform.Peptidoform(
                "[U:1]-ATGPM[U:35]SFLK"
            ),
        }

        for test_string in test_cases.keys():
            assert msms_reader._get_peptidoform(test_string) == test_cases[test_string]

    def test_get_peptide_spectrum_match(self):
        pass

    def test_read_file(self):
        pass
