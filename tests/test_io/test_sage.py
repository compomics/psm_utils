"""Tests for psm_utils.io.sage."""

from psm_utils.io.sage import SageReader
from psm_utils.psm import PSM

test_psm = PSM(
    peptidoform="LQSRPAAPPAPGPGQLTLR/3",
    spectrum_id="controllerType=0 controllerNumber=1 scan=30069",
    run="LQSRPAAPPAPGPGQLTLR",
    collection=None,
    spectrum=None,
    is_decoy=False,
    score=1.2944585,
    qvalue=1.0,
    pep=None,
    precursor_mz=643.0349916987367,
    retention_time=108.2854,
    ion_mobility=None,
    protein_list=["sp|Q99536|VAT1_HUMAN"],
    rank=1,
    source="sage",
    metadata={},
    rescoring_features={
        "expmass": 1926.0815,
        "calcmass": 1926.08,
        "peptide_len": 19.0,
        "missed_cleavages": 0.0,
        "isotope_error": 0.0,
        "precursor_ppm": 0.8239083,
        "fragment_ppm": 0.5347518,
        "hyperscore": 71.78844460255384,
        "delta_next": 71.78844460255384,
        "delta_best": 0.0,
        "delta_rt_model": 0.0,
        "aligned_rt": 0.0,
        "predicted_rt": 0.0,
        "matched_peaks": 22.0,
        "longest_b": 9.0,
        "longest_y": 12.0,
        "longest_y_pct": 0.6315789,
        "matched_intensity_pct": 50.785,
        "scored_candidates": 1.0,
        "poisson": -1.9562811911083433,
        "ms1_intensity": 306146180.0,
        "ms2_intensity": 56930696.0,
    },
)


class TestSageReader:
    def test_iter(self):
        with SageReader("./tests/test_data/results.sage.tsv") as reader:
            for psm in reader:
                psm.provenance_data = {}
                assert psm == test_psm
