"""
Building a proteome peptide set: every k-mer of every sequence in a FASTA.

This is the one part of the package that knows what a proteome is, and the one
part that needs the `proteome` install extra. Prediction imports nothing from
here.
"""

import os
import time

import polars as pl

from netmhcpan_multiprocessor.errors import ConfigError, PeptideInputError
from netmhcpan_multiprocessor.modes.base import PreparedFile, PrepareResult


# The twenty standard amino acids. netMHCpan cannot score a peptide containing
# anything else - X, B, Z, U and the rest - so by default they are dropped here
# rather than left to fail a prediction run much later.
STANDARD_AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'


def prepare_proteome(config) -> PrepareResult:
    """
    Build peptide parquet files from a proteome FASTA, one per requested length.

    Args:
        config (PrepareConfig): the FASTA to read, the folder to write, the
            lengths to generate and the name to give the peptide set.

    Returns:
        PrepareResult: the files written, in ascending length order. The layout
            written is exactly the layout the proteome source discovers.

    Raises:
        ConfigError: if the optional dependencies for reading FASTA are absent.
        PeptideInputError: if the FASTA yields no peptides at some length.
    """
    load_fasta_sequences, polars_kmerize = _load_proteome_dependencies()

    started = time.perf_counter()
    sequences = load_fasta_sequences(str(config.input_path))

    written = []
    for length in config.lengths:
        peptides = polars_kmerize(sequences, n=length).select('peptide_sequence', 'id', 'start', 'end')
        if not config.keep_nonstandard:
            peptides = peptides.filter(
                ~pl.col('peptide_sequence').str.contains(f"[^{STANDARD_AMINO_ACIDS}]")
            )

        # <length>mers/<dataset>.parquet: the layout the proteome source reads,
        # and the reason preparation and prediction need nothing else in common.
        length_folder = config.output_folder / f"{length}mers"
        os.makedirs(length_folder, exist_ok=True)
        output_path = length_folder / f"{config.dataset}.parquet"

        # sink_parquet streams the query straight to disk, so a whole proteome's
        # k-mers never have to be held in memory at once.
        peptides.sink_parquet(output_path)

        peptide_count = pl.scan_parquet(output_path).select(pl.len()).collect().item()
        if peptide_count == 0:
            raise PeptideInputError(
                f"{config.input_path} yielded no {length}mers. Are the sequences shorter than {length}?"
            )
        written.append(PreparedFile(length=length, path=output_path, peptide_count=peptide_count))

    return PrepareResult(
        dataset=config.dataset,
        files=written,
        elapsed_seconds=time.perf_counter() - started,
    )


def _load_proteome_dependencies():
    # Imported on use, not on import, so that an install without the 'proteome'
    # extra can still predict. The message names the extra rather than the
    # packages, because that is what has to be typed to fix it.
    try:
        from shared_helpers.kmers import polars_kmerize
        from shared_helpers.proteomes import load_fasta_sequences
    except ImportError as missing:
        raise ConfigError(
            "Preparing a proteome peptide set needs the optional FASTA dependencies. "
            "Install them with: uv sync --extra proteome "
            "(or pip install 'netmhcpan-multiprocessor[proteome]')"
        ) from missing
    return load_fasta_sequences, polars_kmerize
