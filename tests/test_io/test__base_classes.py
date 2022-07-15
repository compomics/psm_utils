"""Tests for psm_utils.io._base_classes."""

import pytest

from psm_utils.io._base_classes import ReaderBase
from psm_utils.io.exceptions import (
    InvalidModificationDefinitionError,
    InvalidModificationError,
    UnresolvableModificationError,
)


class TestBaseReader:
    def test_validate_modification_definitions(self):
        # Key missing
        with pytest.raises(InvalidModificationDefinitionError):
            ReaderBase.validate_modification_definitions(
                [{"search_engine_label": "TMT6", "proforma_label": "UNIMOD:TMT6plex"}]
            )

        # Invalid label
        with pytest.raises(InvalidModificationError):
            ReaderBase.validate_modification_definitions(
                [{
                    "site": "K|N-term",
                    "search_engine_label": "TMT6",
                    "proforma_label": "+23.0+"
                }]
        )

        # Unresolvable label: Neither mass nor composition can be resolved
        with pytest.raises(UnresolvableModificationError):
            ReaderBase.validate_modification_definitions(
                [{
                    "site": "K|N-term",
                    "search_engine_label": "TMT6",
                    "proforma_label": "U:TMTMTTM"
                }]
        )

        # Invalid label
        with pytest.raises(InvalidModificationError):
            ReaderBase.validate_modification_definitions(
                [{
                    "site": "K|N-term",
                    "search_engine_label": "TMT6",
                    "proforma_label": "+23.0+"
                }]
        )

        # Composition cannot be resolved
        with pytest.warns(
            UserWarning,
            match="Atomic composition for modification `.*` could not be resolved."
        ):
            ReaderBase.validate_modification_definitions(
                [{
                    "site": "K|N-term",
                    "search_engine_label": "TMT6",
                    "proforma_label": "+229.162932"
                }]
        )


