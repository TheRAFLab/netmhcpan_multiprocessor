"""
The proteome mode: peptides that are every k-mer of a set of protein sequences.

Everything proteome-specific lives in this folder, and nothing outside it
imports from here except the registry. Preparation reads FASTA, which needs
dependencies prediction does not, so it is behind a lazy import as well as the
`proteome` install extra. If the preparation step is later moved into its own
tool, this folder loses prepare.py and `preparer` becomes None; nothing else
changes.
"""

from netmhcpan_multiprocessor.modes.base import Mode, PrepareResult
from netmhcpan_multiprocessor.modes.proteome.source import ProteomeSource


class ProteomePreparer:
    """Builds a peptide set by k-merising the sequences in a FASTA file."""

    name = 'proteome'

    def prepare(self, config) -> PrepareResult:
        """
        Build peptide parquet files from a FASTA.

        Args:
            config (PrepareConfig): what to prepare and where to put it.

        Returns:
            PrepareResult: the peptide files written, one per length.

        Raises:
            ConfigError: if the optional dependencies for reading FASTA are absent.
        """
        # Imported here rather than at module scope so that installing without
        # the 'proteome' extra still leaves prediction, and this registry,
        # perfectly usable.
        from netmhcpan_multiprocessor.modes.proteome.prepare import prepare_proteome

        return prepare_proteome(config)


PROTEOME_MODE = Mode(
    name='proteome',
    sources={ProteomeSource.name: ProteomeSource()},
    default_source=ProteomeSource.name,
    preparer=ProteomePreparer(),
)
