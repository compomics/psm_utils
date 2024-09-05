from unittest.mock import mock_open, patch

import pytest

from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.io.flashlfq import FlashLFQReader, FlashLFQWriter
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList

# Sample data for testing
sample_tsv_data = """File Name\tBase Sequence\tFull Sequence\tPeptide Monoisotope Mass\tScan Retention Time\tPrecursor Charge\tProtein Accessions
sample1.raw\tPEPTIDE\tPEPTIDE\t1000.0\t5.0\t2\tP12345;P67890
sample2.raw\tPEPTIDE\tPEPTIDE\t1000.0\t10.0\t2\tP23456|P78901
"""


@pytest.fixture
def valid_psm():
    return PSM(
        peptidoform="PEPTIDE/2",
        spectrum_id="0",
        run="sample1.raw",
        retention_time=5.0,
        protein_list=["P12345", "P67890"],
    )


@pytest.fixture
def invalid_psm_entry():
    return {
        "Full Sequence": "PEPTIDE",
        "Precursor Charge": "2",
        "File Name": "sample1.raw",
        "Scan Retention Time": "5.0",
        "Protein Accessions": None,
    }


def test_flashlfqreader_parse_entry(valid_psm):
    reader = FlashLFQReader("dummy_file.tsv")
    entry = {
        "Full Sequence": "PEPTIDE",
        "Precursor Charge": "2",
        "File Name": "sample1.raw",
        "Scan Retention Time": "5.0",
        "Protein Accessions": "P12345;P67890",
    }
    psm = reader._parse_entry(entry, spectrum_id="0")
    assert psm == valid_psm


def test_flashlfqreader_iterate_over_file():
    with patch("builtins.open", mock_open(read_data=sample_tsv_data)):
        reader = FlashLFQReader("dummy_file.tsv")
        psms = list(reader)
        assert len(psms) == 2
        assert psms[0].run == "sample1.raw"
        assert psms[1].run == "sample2.raw"


def test_flashlfqreader_invalid_entry_handling():
    invalid_data = """File Name\tBase Sequence\tPeptide Monoisotope Mass\tScan Retention Time\tPrecursor Charge\tProtein Accessions
sample1.raw\tPEPTIDE\t1000.0\t5.0
sample2.raw\tPEPTIDE\t1000.0\t10.0
sample3.raw\tPEPTIDE\t1000.0\t15.0
"""
    with patch("builtins.open", mock_open(read_data=invalid_data)):
        reader = FlashLFQReader("dummy_file.tsv")
        with pytest.raises(PSMUtilsIOException):
            psms = list(reader)  # noqa: F841


def test_flashlfqwriter_write_psm(valid_psm):
    with patch("builtins.open", mock_open()) as mocked_file:
        with FlashLFQWriter("dummy_file.tsv") as writer:
            writer.write_psm(valid_psm)
        mocked_file().write.assert_called()
        assert "PEPTIDE" in mocked_file().write.call_args[0][0]


def test_flashlfqwriter_write_file(valid_psm):
    psm_list = PSMList(psm_list=[valid_psm, valid_psm])
    with patch("builtins.open", mock_open()) as mocked_file:
        with FlashLFQWriter("dummy_file.tsv") as writer:
            writer.write_file(psm_list)
        mocked_file().write.assert_called()
        assert mocked_file().write.call_count == 4  # One for header, two for entries, final newline


def test_flashlfqwriter_existing_file():
    # Simulate a file that already exists
    with patch("builtins.open", mock_open(read_data=sample_tsv_data)):
        with FlashLFQWriter("dummy_file.tsv") as writer:
            assert writer.fieldnames == [
                "File Name",
                "Base Sequence",
                "Full Sequence",
                "Peptide Monoisotope Mass",
                "Scan Retention Time",
                "Precursor Charge",
                "Protein Accessions",
            ]


def test_flashlfqwriter_context_manager():
    with patch("builtins.open", mock_open()):
        with FlashLFQWriter("dummy_file.tsv") as writer:
            assert writer._open_file is not None
        assert writer._open_file is None  # Ensure file is closed after context exit
