"""Tests for psm_utils.io.peptide_record."""

import pytest

from psm_utils.io.peptide_record import InvalidPeprecModificationError, PeptideRecord


class TestPeptideRecord:
    def test_peprec_to_proforma(self):
        # Valid cases
        valid_test_cases = [
            (("ACDMEK", ""), "ACDMEK"),
            (("ACDMEK", "-"), "ACDMEK"),
            (("ACDMEK", "4|Oxidation"), "ACDM[Oxidation]EK"),
            (("ACDMEK", "0|Acetylation"), "[Acetylation]-ACDMEK"),
            (("ACDMEK", "-1|Amidation"), "ACDMEK-[Amidation]"),
            (("ACDMEK", "4|Oxidation|-1|Amidation"), "ACDM[Oxidation]EK-[Amidation]"),
            (("ACDMEK", "0|Acetyl|4|Ox|-1|Amide"), "[Acetyl]-ACDM[Ox]EK-[Amide]"),
            (("ACDMEK", "6|Methylation|-1|Amide"), "ACDMEK[Methylation]-[Amide]"),
            (("MCDMEK", "0|Acetyl|1|Oxidation"), "[Acetyl]-M[Oxidation]CDMEK"),
        ]
        for test_in, expected_out in valid_test_cases:
            test_out = PeptideRecord.peprec_to_proforma(*test_in)
            assert test_out == expected_out

        # Invalid case
        with pytest.raises(InvalidPeprecModificationError):
            PeptideRecord.peprec_to_proforma("ACDE", "|")

class TestPeptideRecordReader:
    # TODO!!!
    def test_read_file(self):
        pass