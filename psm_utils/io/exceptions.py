"""Exceptions for psm_utils.io modules."""

from psm_utils.exceptions import PSMUtilsException


class PSMUtilsIOException(PSMUtilsException):
    """Exception in psm_utils.io."""

    pass


class ModificationException(PSMUtilsIOException):
    """Exception while handling or parsing peptide modifications."""

    pass


class InvalidModificationDefinitionError(ModificationException):
    """Exception while handling or parsing modification_definition dictionary."""

    pass


class InvalidModificationError(ModificationException):
    """ProForma modification tag could not be parsed."""

    pass


class UnresolvableModificationError(ModificationException):
    """ProForma modification tag could not be resolved to a controlled vocabulary."""

    pass
