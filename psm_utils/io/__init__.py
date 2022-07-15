"""Parsers for proteomics search results from various search engines."""

from pathlib import Path
from typing import Union

from rich.progress import track

import psm_utils.io.maxquant as maxquant
import psm_utils.io.peptide_record as peptide_record

# TODO: to be completed
READERS = {
    "peprec": peptide_record.PeptideRecordReader,
    "msms": maxquant.MaxquantReader,
}

# TODO: to be completed
WRITERS = {
    "peprec": peptide_record.PeptideRecordWriter,
}


def _infer_filetype(filename: str):
    """Infer filetype from file extension."""
    raise NotImplementedError()


def convert(
    input_filename: Union[str, Path],
    output_filename: Union[str, Path],
    input_filetype: str = "infer",
    output_filetype: str = "infer",
    show_progressbar: bool = False,
):
    """
    Convert a PSM file from one format into another.

    Parameters
    ----------
    input_filename: str
        Path to input file.
    output_filename: str
        Path to output file.
    input_filetype: str, optional
        Input file type. Any of {#TODO: add list}.
    output_filetype: str, optional
        Output file type. Any of {#TODO: add list}.
    show_progressbar: bool, optional
        Show progress bar for conversion process. (default: False)


    Examples
    --------

    Convert a MaxQuant msms.txt file to a MS²PIP peprec file, while inferring
    the applicable file types from the file extensions:

    >>> from psm_utils.io import convert
    >>> convert("msms.txt", "filename_out.peprec")

    Convert a MaxQuant msms.txt file to a MS²PIP peprec file, while explicitly
    specifying both file types:

    >>> convert(
    ...     "filename_in.msms",
    ...     "filename_out.peprec",
    ...     input_filetype="msms",
    ...     output_filetype="peprec"
    ... )

    Note that filetypes can only be inferred for select specific file names and/or
    extensions, such as ``msms.txt`` or `*.peprec`.

    """

    # If needed, infer input and output filetypes
    if input_filetype == "infer":
        input_filetype = _infer_filetype(input_filetype)
    if output_filetype == "infer":
        output_filetype = _infer_filetype(output_filetype)

    reader_cls = READERS[input_filetype]
    writer_cls = WRITERS[output_filetype]

    reader = reader_cls(input_filename)
    writer = writer_cls(output_filename)

    iterator = track(reader) if show_progressbar else reader

    for psm in iterator:
        writer.write(psm)
