"""
The proteome peptide source: a k-merised proteome, one folder per length.

    <input_folder>/8mers/<dataset>.parquet
    <input_folder>/9mers/<dataset>.parquet
    ...

The length comes from the folder name and every peptide is checked against it,
because netMHCpan only reports a length mismatch as a short result, which is a
long way to travel to find out the input was wrong.
"""

import re
from pathlib import Path

import polars as pl

from netmhcpan_multiprocessor.errors import PeptideInputError
from netmhcpan_multiprocessor.modes.base import PeptideFile


LENGTH_FOLDER = re.compile(r'^([0-9]+)mers$')


class ProteomeSource:
    """Reads the `<N>mers/<dataset>.parquet` layout the proteome preparation writes."""

    name = 'proteome'

    def discover(self, input_folder: Path, lengths: tuple[int, ...] | None) -> list[PeptideFile]:
        """
        Find the peptide parquet files to predict.

        Args:
            input_folder (Path): a folder of `<N>mers` subfolders.
            lengths (tuple[int, ...] | None): the lengths wanted, or None for
                every length present.

        Returns:
            list[PeptideFile]: the files to predict, ordered by length then
                dataset. More than one dataset per length folder is allowed;
                each keeps its own name through to the output.

        Raises:
            PeptideInputError: if the folder is missing, holds no length
                folders, or lacks a requested length.
        """
        input_folder = Path(input_folder)
        if not input_folder.is_dir():
            raise PeptideInputError(f"No peptide input folder at {input_folder}")

        found: dict[int, list[Path]] = {}
        for candidate in sorted(input_folder.iterdir()):
            length_match = LENGTH_FOLDER.match(candidate.name)
            if not candidate.is_dir() or not length_match:
                continue
            parquet_files = sorted(candidate.glob('*.parquet'))
            if parquet_files:
                found[int(length_match.group(1))] = parquet_files

        if not found:
            raise PeptideInputError(
                f"No peptide parquet files under {input_folder}. "
                "The proteome source expects <length>mers folders, e.g. 9mers/dataset.parquet"
            )

        # An explicitly requested length that is absent is an error, not a
        # silent skip: a run that quietly predicts three of the four lengths
        # asked for is worse than one that stops.
        if lengths is not None:
            missing = [length for length in lengths if length not in found]
            if missing:
                available = ', '.join(f"{length}mers" for length in sorted(found))
                raise PeptideInputError(
                    f"No peptide files for {', '.join(f'{length}mers' for length in missing)} "
                    f"under {input_folder}. Found: {available}"
                )
            wanted = sorted(set(lengths))
        else:
            wanted = sorted(found)

        return [
            PeptideFile(length=length, path=path, dataset=path.stem)
            for length in wanted
            for path in found[length]
        ]

    def read_peptides(self, peptide_file: PeptideFile, peptide_column: str) -> list[str]:
        """
        Read one peptide parquet, checking it holds what its folder claims.

        Args:
            peptide_file (PeptideFile): the file to read.
            peptide_column (str): the column holding the peptide sequences.

        Returns:
            list[str]: the peptides, in file order. That order is preserved all
                the way to the output parquet, so the two join positionally.

        Raises:
            PeptideInputError: if the column is missing, holds nulls, or holds
                peptides of a length other than the folder's.
        """
        try:
            frame = pl.read_parquet(peptide_file.path, columns=[peptide_column])
        except pl.exceptions.ColumnNotFoundError as missing_column:
            available = ', '.join(pl.scan_parquet(peptide_file.path).collect_schema().names())
            raise PeptideInputError(
                f"{peptide_file.path} has no {peptide_column!r} column. Found: {available}"
            ) from missing_column

        peptides = frame[peptide_column]
        if peptides.len() == 0:
            raise PeptideInputError(f"{peptide_file.path} holds no peptides")
        if peptides.null_count():
            raise PeptideInputError(
                f"{peptide_file.path} has {peptides.null_count()} null peptides in {peptide_column!r}"
            )

        # Checked as a frame-wide min/max rather than per row, so the cost does
        # not scale with how many peptides there are.
        peptide_lengths = peptides.str.len_chars()
        shortest, longest = peptide_lengths.min(), peptide_lengths.max()
        if shortest != peptide_file.length or longest != peptide_file.length:
            observed = f"{shortest}" if shortest == longest else f"{shortest}-{longest}"
            raise PeptideInputError(
                f"{peptide_file.path} is in a {peptide_file.length}mers folder but holds "
                f"peptides of length {observed}"
            )
        return peptides.to_list()
