from cgi import test

import pytest

import psm_utils.io.maxquant as maxquant

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
            "Mass deviations [Da]": "Mass Deviations [Da]",
            "Mass deviations [ppm]": "Mass Deviations [ppm]",
            "Intensity coverage": "Intensity coverage",
            "Reverse": "Reverse",
            "id": "id",
        }

        columns = TEST_COL.copy()

        # Test to get rename dict with default msms
        assert maxquant.MaxquantReader._fix_column_case(columns) == expected_rename_dict

    def test_set_mass_error_unit(self):
        msms_reader = maxquant.MaxquantReader(
            "./tests/test_data/test_msms.txt"
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
        