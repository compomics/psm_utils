import numpy as np

from psm_utils.io.proteoscape import _parse_peptidoform


def test_parse_peptidoform():
    test_cases = [
        ("ACDMEK", np.array([]), np.array([]), 2, "ACDMEK/2"),
        ("ACDMEK", np.array([15.99]), np.array([3]), 2, "ACDM[+15.99]EK/2"),
        ("ACDMEK", np.array([57.02, 15.99]), np.array([1, 3]), 2, "AC[+57.02]DM[+15.99]EK/2"),
        ("ACDMEK", np.array([42.01]), np.array([-1]), 2, "[+42.01]-ACDMEK/2"),
        ("ACDMEK", np.array([-0.98]), np.array([6]), 2, "ACDMEK-[-0.98]/2"),
        ("ACDMEK", np.array([42.01, -0.98]), np.array([-1, 6]), 2, "[+42.01]-ACDMEK-[-0.98]/2"),
    ]

    for peptide, ptms, ptm_locations, precursor_charge, expected in test_cases:
        assert _parse_peptidoform(peptide, ptms, ptm_locations, precursor_charge) == expected
