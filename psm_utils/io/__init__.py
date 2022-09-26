"""Parsers for proteomics search results from various search engines."""

import re
from pathlib import Path
from typing import Union

from rich.progress import track

import psm_utils.io.idxml as idxml
import psm_utils.io.maxquant as maxquant
import psm_utils.io.msamanda as msamanda
import psm_utils.io.mzid as mzid
import psm_utils.io.peptide_record as peptide_record
import psm_utils.io.percolator as percolator
import psm_utils.io.tsv as tsv
import psm_utils.io.xtandem as xtandem
from psm_utils.io.exceptions import PSMUtilsIOException

# TODO: to be completed
FILETYPES = {
    "idxml": {
        "reader": idxml.IdXMLReader,
        "writer": None,
        "extension": ".idXML",
        "filename_pattern": r"^.*\.idxml$",
    },
    "msms": {
        "reader": maxquant.MSMSReader,
        "writer": None,
        "extension": "_msms.txt",
        "filename_pattern": r"^msms\.txt$",
    },
    "mzid": {
        "reader": mzid.MzidReader,
        "writer": mzid.MzidWriter,
        "extension": ".mzid",
        "filename_pattern": r"^.*\.(?:(?:mzidentml)|(?:mzid))$",
    },
    "peprec": {
        "reader": peptide_record.PeptideRecordReader,
        "writer": peptide_record.PeptideRecordWriter,
        "extension": ".peprec.txt",
        "filename_pattern": r"(^.*\.peprec(?:\.txt)?$)|(?:^peprec\.txt$)",
    },
    "percolator": {
        "reader": percolator.PercolatorTabReader,
        "writer": percolator.PercolatorTabWriter,
        "extension": ".percolator.txt",
        "filename_pattern": r"^.*\.(?:(?:pin)|(?:pout))$",
    },
    "tsv": {
        "reader": tsv.TSVReader,
        "writer": tsv.TSVWriter,
        "extension": ".tsv",
        "filename_pattern": r"^.*\.tsv$",
    },
    "xtandem": {
        "reader": xtandem.XTandemReader,
        "writer": None,
        "extension": ".t.xml",
        "filename_pattern": r"^.*\.t\.xml$",
    },
    "msamanda": {
        "reader": msamanda.MSAmandaReader,
        "writer": None,
        "extension": ".csv",
        "filename_pattern": r"^.*\.csv$",
    },
}
READERS = {k: v["reader"] for k, v in FILETYPES.items() if v["reader"]}
WRITERS = {k: v["writer"] for k, v in FILETYPES.items() if v["writer"]}


def _infer_filetype(filename: str):
    """Infer filetype from filename."""
    for filetype, properties in FILETYPES.items():
        if re.fullmatch(properties["filename_pattern"], filename, flags=re.IGNORECASE):
            return filetype
    else:
        raise PSMUtilsIOException("Could not infer filetype.")


def read_file(filename: Union[str, Path], *args, filetype: str = "infer", **kwargs):
    """
    Read PSM file into :py:class:`~psm_utils.psmlist.PSMList`.

    Parameters
    ----------
    filename: str
        Path to file.
    filetype: str, optional
        File type. Any PSM file type with read support. See
        :ref:`Supported file formats`.
    *args : tuple
        Additional arguments are passed to the :py:class:`psm_utils.io` reader.
    **kwargs : dict, optional
        Additional keyword arguments are passed to the :py:class:`psm_utils.io` reader.
    """
    if filetype == "infer":
        filetype = _infer_filetype(filetype)

    reader_cls = READERS[filetype]
    reader = reader_cls(filename, *args, **kwargs)
    psm_list = reader.read_file()
    return psm_list


def write_file():
    # TODO
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
        File type. Any PSM file type with read support. See
        :ref:`Supported file formats`.
    output_filetype: str, optional
        File type. Any PSM file type with write support. See
        :ref:`Supported file formats`.
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

    # Remove file if already exists to avoid appending:
    if Path(output_filename).is_file():
        Path(output_filename).unlink()

    reader = reader_cls(input_filename)
    if output_filetype == "mzid":
        writer = writer_cls(output_filename)
        writer.write_file(
            reader.read_file(), disable_progressbar=(not show_progressbar)
        )

    else:
        iterator = (
            track(
                reader,
                show_speed=False,
                description="[green]Converting file",
            )
            if show_progressbar
            else reader
        )

        for psm in reader:
            example_psm = psm
            break
        writer = writer_cls(output_filename, example_psm=example_psm, mode="write")

        with writer:
            for psm in iterator:
                writer.write_psm(psm)
