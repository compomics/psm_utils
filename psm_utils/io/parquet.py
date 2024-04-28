"""
Reader and writer for a simple, lossless psm_utils Parquet format.

Similar to the :py:mod:`psm_utils.io.tsv` module, this module provides a reader and writer
for :py:class:`~psm_utils.psm_list.PSMList` objects in a lossless manner. However, Parquet provides
better performance and storage efficiency compared to TSV, and is recommended for large datasets.

"""

from __future__ import annotations

from pathlib import Path
from typing import Union

import pyarrow as pa
import pyarrow.parquet as pq
from pydantic import ValidationError

from psm_utils.io._base_classes import ReaderBase, WriterBase
from psm_utils.io.exceptions import PSMUtilsIOException
from psm_utils.psm import PSM
from psm_utils.psm_list import PSMList


class ParquetReader(ReaderBase):
    def __init__(self, path: Union[str, Path], *args, **kwargs):
        """
        Reader for Parquet files.

        Parameters
        ----------
        path : Union[str, Path]
            Path to the Parquet file.

        """
        self.path = path

    def __iter__(self):
        with pq.ParquetFile(self.path) as reader:
            for batch in reader.iter_batches():
                for row in batch.to_pylist():
                    # Convert map columns (rendered as lists of tuples) to dictionaries
                    row["metadata"] = dict(row["metadata"] or {})
                    row["provenance_data"] = dict(row["provenance_data"] or {})
                    row["rescoring_features"] = dict(row["rescoring_features"] or {})

                    # Convert to PSM object and yield
                    try:
                        yield PSM(**row)
                    except ValidationError as e:
                        raise PSMUtilsIOException(f"Error while parsing row {row}:\n{e}")


class ParquetWriter(WriterBase):
    def __init__(self, path: Union[str, Path], chunk_size: int = 1e6, *args, **kwargs):
        """
        Writer for Parquet files.

        Parameters
        ----------
        path : Union[str, Path]
            Path to the Parquet file.
        chunk_size : int
            Number of PSMs to write in a single batch. Default is 1e6.

        """
        self.path = path
        self.chunk_size = chunk_size

        self._writer = None
        self._psm_cache = []

    def __enter__(self):
        self._writer = pq.ParquetWriter(self.path, schema=SCHEMA)
        return self

    def __exit__(self, *args, **kwargs):
        self._flush()
        self._writer.close()

    def write_psm(self, psm: PSM):
        """Write a single PSM to the Parquet file."""
        self._psm_cache.append(self._psm_to_entry(psm))
        if len(self._psm_cache) > self.chunk_size:
            self._flush()

    def write_file(self, psm_list: PSMList):
        """Write a list of PSMs to the Parquet file."""
        with self:
            for psm in psm_list:
                self.write_psm(psm)

    @staticmethod
    def _psm_to_entry(psm: PSM) -> dict:
        """Convert a PSM object to a dictionary suitable for writing to Parquet."""
        psm_dict = dict(psm)
        psm_dict["peptidoform"] = str(psm.peptidoform)
        return psm_dict

    def _flush(self):
        """Write the cached PSMs to the Parquet file."""
        if not self._psm_cache:
            return
        table = pa.Table.from_pylist(self._psm_cache, schema=SCHEMA)
        self._writer.write_table(table)
        self._psm_cache = []


SCHEMA = pa.schema(
    [
        ("peptidoform", pa.string()),
        ("spectrum_id", pa.string()),
        ("run", pa.string()),
        ("collection", pa.string()),
        ("spectrum", pa.string()),
        ("is_decoy", pa.bool_()),
        ("score", pa.float32()),
        ("qvalue", pa.float32()),
        ("pep", pa.float32()),
        ("precursor_mz", pa.float32()),
        ("retention_time", pa.float32()),
        ("ion_mobility", pa.float32()),
        ("protein_list", pa.list_(pa.string())),
        ("rank", pa.int32()),
        ("source", pa.string()),
        ("provenance_data", pa.map_(pa.string(), pa.string())),
        ("metadata", pa.map_(pa.string(), pa.string())),
        ("rescoring_features", pa.map_(pa.string(), pa.float32())),
    ]
)
