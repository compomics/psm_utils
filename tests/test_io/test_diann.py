"""Tests for psm_utils.io.diann."""

from psm_utils.io.diann import DIANNReader
from psm_utils.psm import PSM

test_psm = PSM(
    peptidoform="AAAAEINVKDPKEDLETSVVDEGR/4",
    spectrum_id="116903",
    run="LFQ_Orbitrap_AIF_Yeast_03",
    collection=None,
    spectrum=None,
    is_decoy=False,
    score=0.995107,
    qvalue=0.000548193,
    pep=0.0104343,
    precursor_mz=None,
    retention_time=75.2574,
    ion_mobility=0,
    protein_list=["P38156"],
    rank=1,
    source="diann",
    metadata={},
    rescoring_features={
        "CScore": 0.995107,
        "RT": 75.2574,
        "Predicted.RT": 75.2713,
        "iRT": 33.9222,
        "Predicted.iRT": 33.8999,
        "Ms1.Profile.Corr": 0.347567,
        "Ms1.Area": 849940,
        "IM": 0,
        "iIM": 0,
        "Predicted.IM": 0,
        "Predicted.iIM": 0,
    },
)


class TestDIANNReader:
    def test_iter(self):
        with DIANNReader("./tests/test_data/test_diann.tsv") as reader:
            for psm in reader:
                psm.provenance_data = {}
                assert psm == test_psm
