"""
Splitting a peptide set into chunk files.

netMHCpan's -p option takes a file of one peptide per line, and each process is
single threaded, so the only way to use more than one core is to hand several
processes a piece of the work each. Chunks are written into a temporary folder
and never persisted: the peptide parquet is the input of record.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from netmhcpan_multiprocessor.errors import ConfigError


@dataclass(frozen=True)
class PeptideChunk:
    """
    One chunk file handed to a single netMHCpan process.

    Attributes:
        path (Path): the chunk file, one peptide per line.
        peptide_count (int): how many peptides it holds, so a short result can
            be caught.
    """
    path: Path
    peptide_count: int


def write_chunk_files(
    peptides: Sequence[str],
    chunk_size: int,
    chunk_folder: Path,
    name: str,
) -> list[PeptideChunk]:
    """
    Split a list of peptides into chunk files small enough to spread across workers.

    Args:
        peptides (Sequence[str]): the peptides, in the order they are to be predicted.
        chunk_size (int): the maximum number of peptides per chunk.
        chunk_folder (Path): a directory to write the chunk files into.
        name (str): a name for the chunk files, unique within the folder.

    Returns:
        list[PeptideChunk]: each chunk in input order. Reassembling the chunks in
            this order reproduces the input peptide order exactly, which is what
            lets the output parquet be joined back positionally.

    Raises:
        ConfigError: if chunk_size is less than 1.
    """
    if chunk_size < 1:
        raise ConfigError(f"chunk_size must be at least 1, got {chunk_size}")

    chunk_folder = Path(chunk_folder)
    chunks = []
    for chunk_index, start in enumerate(range(0, len(peptides), chunk_size)):
        chunk_peptides = peptides[start:start + chunk_size]
        chunk_path = chunk_folder / f"{name}__chunk{chunk_index:04d}.txt"
        chunk_path.write_text('\n'.join(chunk_peptides))
        chunks.append(PeptideChunk(path=chunk_path, peptide_count=len(chunk_peptides)))
    return chunks
