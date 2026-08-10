"""
The registry of peptide set modes and their sources.

Keyed by (mode, source). Today it holds one pair:

    proteome  ->  {proteome}

An immunopeptidome or display screen mode is a new entry here plus a source
object; an input type that is not tied to one kind of collection, such as a
pre-processed peptide parquet set, is one source object registered under every
mode. Neither case touches the predictor.

Deliberately a dict and two lookups rather than plugin machinery: there is one
pair to register, and the cost of the abstraction should stay proportional
to that.
"""

from netmhcpan_multiprocessor.errors import (
    PreparationNotSupportedError,
    UnknownModeError,
    UnknownSourceError,
)
from netmhcpan_multiprocessor.modes.base import (
    Mode,
    PeptideFile,
    PeptideSource,
    Preparer,
    PreparedFile,
    PrepareResult,
)
from netmhcpan_multiprocessor.modes.proteome import PROTEOME_MODE


MODES: dict[str, Mode] = {
    PROTEOME_MODE.name: PROTEOME_MODE,
}

DEFAULT_MODE = PROTEOME_MODE.name


def registered_modes() -> tuple[str, ...]:
    """
    List the registered modes.

    Returns:
        tuple[str, ...]: every mode name, in registration order.
    """
    return tuple(MODES)


def get_mode(mode: str | Mode) -> Mode:
    """
    Look a mode up by name.

    Args:
        mode (str | Mode): the mode name, or a Mode, which is returned as it is.

    Returns:
        Mode: the registered mode.

    Raises:
        UnknownModeError: if no mode of that name is registered.
    """
    if isinstance(mode, Mode):
        return mode
    try:
        return MODES[mode]
    except KeyError:
        raise UnknownModeError(
            f"Unknown mode {mode!r}. Registered modes: {', '.join(registered_modes())}"
        ) from None


def registered_sources(mode: str | Mode) -> tuple[str, ...]:
    """
    List the sources registered for a mode.

    Args:
        mode (str | Mode): the mode, by name or object.

    Returns:
        tuple[str, ...]: every source name valid for that mode.

    Raises:
        UnknownModeError: if no mode of that name is registered.
    """
    return tuple(get_mode(mode).sources)


def resolve_source(mode: str | Mode, source: str | None = None) -> PeptideSource:
    """
    Resolve a mode and source onto the input handler to read the peptides with.

    Args:
        mode (str | Mode): the mode, by name or object.
        source (str | None): the source name, or None for the mode's default.

    Returns:
        PeptideSource: the source to discover and read peptides with.

    Raises:
        UnknownModeError: if no mode of that name is registered.
        UnknownSourceError: if that mode has no such source. The message lists
            the sources for the mode asked for, not every source everywhere,
            so naming a source that belongs to another mode says so plainly.
    """
    resolved_mode = get_mode(mode)
    if source is None:
        return resolved_mode.sources[resolved_mode.default_source]
    try:
        return resolved_mode.sources[source]
    except KeyError:
        raise UnknownSourceError(
            f"Mode {resolved_mode.name!r} has no source {source!r}. "
            f"Sources for {resolved_mode.name!r}: {', '.join(resolved_mode.sources)}"
        ) from None


def get_preparer(mode: str | Mode) -> Preparer:
    """
    Get a mode's preparation step.

    Args:
        mode (str | Mode): the mode, by name or object.

    Returns:
        Preparer: the mode's preparer.

    Raises:
        UnknownModeError: if no mode of that name is registered.
        PreparationNotSupportedError: if the mode has no preparation step. A
            mode is expected to have a source before it has a preparer, and to
            keep working if its preparer is later removed.
    """
    resolved_mode = get_mode(mode)
    if resolved_mode.preparer is None:
        preparable = ', '.join(name for name, each in MODES.items() if each.preparer is not None)
        raise PreparationNotSupportedError(
            f"Mode {resolved_mode.name!r} has no preparation step. "
            f"Modes that can prepare: {preparable or 'none'}"
        )
    return resolved_mode.preparer


__all__ = [
    'DEFAULT_MODE',
    'MODES',
    'Mode',
    'PeptideFile',
    'PeptideSource',
    'PrepareResult',
    'PreparedFile',
    'Preparer',
    'get_mode',
    'get_preparer',
    'registered_modes',
    'registered_sources',
    'resolve_source',
]
