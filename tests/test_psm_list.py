import numpy as np
import pytest

from psm_utils import PSM, PSMList

sample_psm_list = [
    PSM(peptidoform="ACDK", spectrum_id=1, score=140.2),
    PSM(peptidoform="CDEFR", spectrum_id=2, score=132.9),
    PSM(peptidoform="DEM[Oxidation]K", spectrum_id=3, score=55.7),
]


class TestPSMList:
    def test___init__(self):
        psm_list = PSMList(psm_list=sample_psm_list)

    def test___iter__(self):
        for psm in PSMList(psm_list=sample_psm_list):
            print(psm)

    def test___len__(self):
        psm_list = PSMList(psm_list=sample_psm_list)
        assert len(psm_list) == 3

    def test___get_item__(self):
        psm_list = PSMList(psm_list=sample_psm_list)

        # Single index
        assert psm_list[0] == PSM(peptidoform="ACDK", spectrum_id=1, score=140.2)

        # Slice
        assert psm_list[0:2] == PSMList(
            psm_list=[
                PSM(peptidoform="ACDK", spectrum_id=1, score=140.2),
                PSM(peptidoform="CDEFR", spectrum_id=2, score=132.9),
            ]
        )

        # PSM property as array
        np.testing.assert_equal(psm_list["spectrum_id"], np.array([1, 2, 3]))

        # Multiple PSM properties as 2D array
        np.testing.assert_equal(
            psm_list[["spectrum_id", "score"]],
            np.array([[1, 140.2], [2, 132.9], [3, 55.7]]),
        )

        # Index by multiple indices
        psm_list[0, 2] == PSMList(
            psm_list=[
                PSM(peptidoform="ACDK", spectrum_id=1, score=140.2),
                PSM(peptidoform="DEM[Oxidation]K", spectrum_id=3, score=55.7),
            ]
        )

        # Unsupported indexing type
        with pytest.raises(TypeError):
            psm_list[0.21]

    def test___set_item__(self):
        psm_list = PSMList(psm_list=sample_psm_list)

        # Correct usage
        psm_list["score"] = [10, 20, 30]
        np.testing.assert_equal(psm_list["score"], np.array([10, 20, 30]))

        # Wrong length
        with pytest.raises(ValueError):
            psm_list["score"] = [10, 20, 30, 40]
