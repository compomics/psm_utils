"""Tests for psm_utils.io.sage."""

import pytest

from psm_utils.io.sage import SageParquetReader, SageTSVReader
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
    ion_mobility=0.0,
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
        "fragment_ppm": 0.503857,
        "hyperscore": 72.26591573806016,
        "delta_next": 72.26591573806016,
        "delta_best": 0.0,
        "delta_rt_model": 0.993444,
        "aligned_rt": 0.993444,
        "predicted_rt": 0.0,
        "predicted_mobility": 0.0,
        "matched_peaks": 22.0,
        "longest_b": 9.0,
        "longest_y": 12.0,
        "longest_y_pct": 0.6315789,
        "matched_intensity_pct": 64.770966,
        "scored_candidates": 1.0,
        "poisson": -1.9562811911083433,
        "ms2_intensity": 72609170.0,
    },
)


class TestSageTSVReader:
    def test_iter(self):
        with SageTSVReader("./tests/test_data/results.sage.tsv") as reader:
            for psm in reader:
                psm.provenance_data = {}
                assert psm == test_psm


class TestSageParquetReader:
    def test_iter(self):
        with SageParquetReader("./tests/test_data/results.sage.parquet") as reader:
            # Parquet results in float precision differences, so pytest.approx is used, which does
            # not support objects with nested dicts.
            for psm in reader:
                psm_dict = dict(psm)
                test_psm_dict = dict(test_psm)

                # Nested dicts
                assert psm_dict.pop("rescoring_features", {}) == pytest.approx(
                    test_psm_dict.pop("rescoring_features", {})
                )
                assert psm_dict.pop("metadata", {}) == test_psm_dict.pop("metadata", {})
                psm_dict.pop("provenance_data", {})

                # Remaining keys
                for k, v in psm_dict.items():
                    if isinstance(v, float):
                        assert v == pytest.approx(test_psm_dict[k])
                    else:
                        assert v == test_psm_dict[k]
