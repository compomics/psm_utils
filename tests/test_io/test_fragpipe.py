"""Tests for psm_utils.io.fragpipe."""

from psm_utils.io.fragpipe import FragpipeReader
from psm_utils.psm import PSM

test_psm = PSM(
    peptidoform="TGAPNNGQYGADNGNPNGER/2",
    spectrum_id="00001",
    run="LFQ_Orbitrap_AIF_Yeast_01_Q1",
    collection=None,
    spectrum=None,
    is_decoy=False,
    score=57.2940,
    qvalue=None,
    pep=1 - 1.0000,
    precursor_mz=1001.9342,
    retention_time=2432.1640,
    ion_mobility=None,
    protein_list=["sp|P40159|YNU8_YEAST"],
    rank=1,
    source="fragpipe",
    metadata={},
    rescoring_features={
        "Peptide Length": 20,
        "Retention": 2432.1640,
        "Observed Mass": 2001.8539,
        "Observed M/Z": 1001.9342,
        "Calculated Peptide Mass": 2001.8524,
        "Calculated M/Z": 1001.9335,
        "Delta Mass": 0.0002,
        "Hyperscore": 57.2940,
        "Number of Missed Cleavages": 0,
    },
)


class TestFragpipeReader:
    def test_iter(self):
        with FragpipeReader("./tests/test_data/test_fragpipe.tsv") as reader:
            for psm in reader:
                psm.provenance_data = {}
                assert psm == test_psm
