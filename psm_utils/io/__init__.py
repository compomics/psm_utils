"""Parsers for proteomics search results from various search engines."""

import re
from pathlib import Path
from tempfile import NamedTemporaryFile
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
from psm_utils.io._base_classes import WriterBase
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.psm_list import PSMList

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
        "filename_pattern": r"^.*msms\.txt$",
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
        "filename_pattern": r"^.*(?:_|\.)msamanda.csv$",
    },
}
READERS = {k: v["reader"] for k, v in FILETYPES.items() if v["reader"]}
WRITERS = {k: v["writer"] for k, v in FILETYPES.items() if v["writer"]}


def _infer_filetype(filename: str):
    """Infer filetype from filename."""
    for filetype, properties in FILETYPES.items():
        if re.fullmatch(
            properties["filename_pattern"], str(filename), flags=re.IGNORECASE
        ):
            return filetype
    else:
        raise PSMUtilsIOException("Could not infer filetype.")


def _supports_write_psm(writer: WriterBase):
    """Check if writer supports write_psm method."""
    try:
        with NamedTemporaryFile(delete=False) as temp_file:
            with writer(temp_file.name) as writer_instance:
                temp_file.close()
                writer_instance.write_psm(None)
                temp_file.delete()
    except NotImplementedError:
        supports_write_psm = False
    except AttributeError:  # `None` is not valid PSM
        supports_write_psm = True
    else:
        supports_write_psm = True
    return supports_write_psm


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
        filetype = _infer_filetype(filename)

    reader_cls = READERS[filetype]
    reader = reader_cls(filename, *args, **kwargs)
    psm_list = reader.read_file()
    return psm_list


def write_file(
    psm_list: PSMList,
    filename: Union[str, Path],
    *args,
    filetype: str = "infer",
    show_progressbar: bool = False,
    **kwargs
):
    """
    Write :py:class:`~psm_utils.psmlist.PSMList` to PSM file.

    Parameters
    ----------
    psm_list: PSMList
        PSM list to be written.
    filename: str
        Path to file.
    filetype: str, optional
        File type. Any PSM file type with read support. See
        :ref:`Supported file formats`.
    show_progressbar: bool, optional
        Show progress bar for conversion process. (default: False)
    *args : tuple
        Additional arguments are passed to the :py:class:`psm_utils.io` writer.
    **kwargs : dict, optional
        Additional keyword arguments are passed to the :py:class:`psm_utils.io` writer.
    """
    if filetype == "infer":
        filetype = _infer_filetype(filename)
    writer_cls = WRITERS[filetype]

    # Remove file if already exists to avoid appending:
    if Path(filename).is_file():
        Path(filename).unlink()

    # Get example PSM, instantiate writer, write
    example_psm = psm_list[0]
    with writer_cls(
        filename,
        *args,
        example_psm=example_psm,
        mode="write",
        show_progressbar=show_progressbar,
        **kwargs,
    ) as writer:
        writer.write_file(psm_list)


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
    extensions, such as ``msms.txt`` or ``*.peprec``.

    """

    # If needed, infer input and output filetypes
    if input_filetype == "infer":
        input_filetype = _infer_filetype(input_filename)
    if output_filetype == "infer":
        output_filetype = _infer_filetype(output_filename)

    reader_cls = READERS[input_filetype]
    writer_cls = WRITERS[output_filetype]

    # Remove file if already exists to avoid appending:
    if Path(output_filename).is_file():
        Path(output_filename).unlink()

    reader = reader_cls(input_filename)

    if _supports_write_psm(writer_cls):
        # Setup iterator, potentially with progress bar
        iterator = (
            track(reader, description="[green]Converting file")
            if show_progressbar
            else reader
        )

        # Get example PSM and instantiate writer
        for psm in reader:
            example_psm = psm
            break
        writer = writer_cls(output_filename, example_psm=example_psm, mode="write")

        # Convert
        with writer:
            for psm in iterator:
                writer.write_psm(psm)

    # First read full PSM list, then write file at once
    elif writer_cls == mzid.MzidWriter:
        writer = writer_cls(output_filename)
        writer.write_file(reader.read_file(), show_progressbar=show_progressbar)
    else:
        writer = writer_cls(output_filename)
        writer.write_file(reader.read_file())
