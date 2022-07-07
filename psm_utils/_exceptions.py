"""psm_utils exceptions."""


class PSMUtilsException(Exception):
    """General psm_utils exception."""

    pass


class PeptidoformError(PSMUtilsException):
    """Error while parsing peptidoform."""

    pass


class ModificationParsingException(PSMUtilsException):
    """Identification file parsing error."""
