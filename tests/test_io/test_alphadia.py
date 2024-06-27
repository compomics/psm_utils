"""Tests for psm_utils.io.alphadia."""

from psm_utils.io.alphadia import AlphaDIAReader
from psm_utils.psm import PSM

test_psm = PSM(
    peptidoform="LMIEM[Oxidation]DGTENK/2",
    spectrum_id="79426",
    run="LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_03",
    collection=None,
    spectrum="79426",
    is_decoy=False,
    score=170.287918,
    qvalue=0.000000,
    pep=0.000005,
    precursor_mz=648.794128,
    retention_time=3111.141602,
    ion_mobility=0.000001,
    protein_list=["P06733"],
    rank=1,
    source="alphadia",
    metadata={},
    rescoring_features={
        "rt_observed": 3111.141602,
        "mobility_observed": 0.000001,
        "mz_observed": 648.794128,
        "score": 170.287918,
        "charge": 2,
        "delta_rt": -39.528809,
    },
)


class TestAlphaDIAReader:
    def test_iter(self):
        with AlphaDIAReader("./tests/test_data/test_alphadia.tsv") as reader:
            for psm in reader:
                psm.provenance_data = {}
                assert psm == test_psm
