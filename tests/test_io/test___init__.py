"""Tests for psm_utils.io.__init__."""

import pytest

from psm_utils.io import _infer_filetype
from psm_utils.io.exceptions import PSMUtilsIOException


def test__infer_filetype():
    test_cases = [
        ("name.idxml", "idxml"),
        ("name.idXML", "idxml"),
        ("msms.txt", "msms"),
        ("name.peprec", "peprec"),
        ("name.peprec.txt", "peprec"),
        ("name.pin", "percolator"),
        ("name.pout", "percolator"),
        ("name.t.xml", "xtandem"),
        ("name.msamanda.csv", "msamanda"),
        ("name_msamanda.csv", "msamanda"),
    ]
    for test_in, expected_out in test_cases:
        assert _infer_filetype(test_in) == expected_out

    test_cases = ["name.txt", "name.idxml.txt"]
    for test_in in test_cases:
        with pytest.raises(PSMUtilsIOException):
            _infer_filetype(test_in)
