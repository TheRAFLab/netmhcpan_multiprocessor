"""
The boundary between how a peptide set arrives and how it is predicted.

Nothing downstream of the peptide list cares where the peptides came from:
netMHCpan sees a list of equal-length strings either way. What differs between
kinds of peptide collection is how the set is laid out on disk, how its lengths
are established, and what has to happen before it exists at all. That is the
whole job of this layer, and it is the only part a new kind of input touches.

Two axes:

    mode    the kind of peptide collection - a proteome, an immunopeptidome,
            a display screen. It carries the sources valid for that kind, and
            its preparation step if it has one.
    source  a concrete input handler registered for a mode. A mode may hold
            several: a future pre-processed-parquet reader would be one source
            object registered under every mode, not a copy per mode.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, runtime_checkable


@dataclass(frozen=True)
class PeptideFile:
    """
    One file of peptides of a single length, found by a source.

    Attributes:
        length (int): the peptide length this file holds.
        path (Path): the file itself.
        dataset (str): the name of the peptide set, taken from the file name.
            It names the output file, so it travels with the peptides rather
            than being supplied separately.
    """
    length: int
    path: Path
    dataset: str


@runtime_checkable
class PeptideSource(Protocol):
    """A way of finding and reading a peptide set from an input folder."""

    name: str

    def discover(self, input_folder: Path, lengths: tuple[int, ...] | None) -> list[PeptideFile]:
        """
        Find the peptide files to predict.

        Args:
            input_folder (Path): the folder holding the peptide set.
            lengths (tuple[int, ...] | None): the lengths wanted, or None for
                every length present.

        Returns:
            list[PeptideFile]: the files to predict, one per length per dataset.

        Raises:
            PeptideInputError: if the folder holds nothing usable, or a
                requested length is absent.
        """
        ...

    def read_peptides(self, peptide_file: PeptideFile, peptide_column: str) -> list[str]:
        """
        Read one file's peptides, in the order they are to be predicted.

        Args:
            peptide_file (PeptideFile): the file to read.
            peptide_column (str): the column holding the peptide sequences.

        Returns:
            list[str]: the peptides, in file order.

        Raises:
            PeptideInputError: if the column is missing, or the peptides are not
                what the file claims.
        """
        ...


@runtime_checkable
class Preparer(Protocol):
    """A way of producing a peptide set that a source can then read."""

    name: str

    def prepare(self, config) -> 'PrepareResult':
        """
        Build the peptide files for one dataset.

        Args:
            config (PrepareConfig): what to prepare and where to put it.

        Returns:
            PrepareResult: the files written.
        """
        ...


@dataclass(frozen=True)
class PreparedFile:
    """
    One peptide file written by a preparation step.

    Attributes:
        length (int): the peptide length the file holds.
        path (Path): the parquet written.
        peptide_count (int): how many peptides it holds.
    """
    length: int
    path: Path
    peptide_count: int


@dataclass(frozen=True)
class PrepareResult:
    """
    What a preparation run produced.

    Attributes:
        dataset (str): the name given to the peptide set.
        files (list[PreparedFile]): the peptide files written, one per length.
        elapsed_seconds (float): how long the run took.
    """
    dataset: str
    files: list[PreparedFile]
    elapsed_seconds: float

    @property
    def peptide_count(self) -> int:
        """int: the total number of peptides written across all lengths."""
        return sum(prepared.peptide_count for prepared in self.files)


@dataclass(frozen=True)
class Mode:
    """
    A kind of peptide collection, and the input handlers registered for it.

    Attributes:
        name (str): the mode's name, as given to --mode.
        sources (dict[str, PeptideSource]): the sources valid for this mode.
        default_source (str): the source used when none is named.
        preparer (Preparer | None): the mode's preparation step, or None if it
            has none. A mode is expected to have a source before it has a
            preparer, and to keep working if its preparer is later removed.
    """
    name: str
    sources: dict[str, PeptideSource]
    default_source: str
    preparer: 'Preparer | None' = None
