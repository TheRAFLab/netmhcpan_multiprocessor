"""
Configuration for the two things this package does.

One frozen dataclass per command, each the single source of truth for it. The
CLI builds one of these and calls exactly the function a pipeline would, so
there is no path through the CLI that a caller cannot also take.

Both accept a plain dictionary, because a pipeline usually already holds its
settings as one, and both validate at construction: a folder that is not there,
a mode that is not registered or an allele that netMHCpan has never heard of
should stop the run in its first second, not after the first long prediction.
"""

import json
import os
import shutil
import tomllib
from dataclasses import dataclass, fields
from pathlib import Path

from netmhcpan_multiprocessor.errors import ConfigError, NetMHCpanNotFoundError
from netmhcpan_multiprocessor.modes import DEFAULT_MODE, get_preparer, resolve_source
from netmhcpan_multiprocessor.parsing import DEFAULT_SCORES, get_score_set


# The environment variable checked when no path is configured. Named for the
# program rather than the package, since it is the program being located.
NETMHCPAN_ENVIRONMENT_VARIABLE = 'NETMHCPAN'

# Singular spellings accepted in a config dictionary, because 'allele' and
# 'length' read better than their plurals when only one is being given.
FIELD_ALIASES = {
    'allele': 'alleles',
    'length': 'lengths',
    'netmhcpan': 'netmhcpan_path',
    'input': 'input_path',
}


def resolve_netmhcpan_path(netmhcpan_path: Path | str | None = None) -> Path:
    """
    Find the netMHCpan executable to run.

    Args:
        netmhcpan_path (Path | str | None): an explicitly configured path, or
            None to search the environment and then PATH.

    Returns:
        Path: an executable netMHCpan.

    Raises:
        NetMHCpanNotFoundError: if no executable could be found, or the one
            configured is not executable.
    """
    if netmhcpan_path:
        candidate = Path(netmhcpan_path).expanduser()
        if not candidate.is_file():
            raise NetMHCpanNotFoundError(f"No netMHCpan executable at {candidate}")
        if not os.access(candidate, os.X_OK):
            raise NetMHCpanNotFoundError(f"{candidate} is not executable")
        return candidate

    from_environment = os.environ.get(NETMHCPAN_ENVIRONMENT_VARIABLE)
    if from_environment:
        # Resolved through the same checks, so a stale environment variable is
        # reported as plainly as a mistyped option.
        return resolve_netmhcpan_path(from_environment)

    on_path = shutil.which('netMHCpan')
    if on_path:
        return Path(on_path)

    raise NetMHCpanNotFoundError(
        "Could not find netMHCpan. Set it with the --netmhcpan option, the "
        f"netmhcpan_path config key, the {NETMHCPAN_ENVIRONMENT_VARIABLE} environment "
        "variable, or put netMHCpan on PATH."
    )


def _coerce_lengths(lengths) -> tuple[int, ...] | None:
    if lengths is None:
        return None
    # A string is iterable, so without this '11' would be read character by
    # character as the single length 1 - and '9' would work, which is exactly
    # how a caller ends up trusting the broken form.
    if isinstance(lengths, (int, str, bytes)):
        lengths = [lengths]
    try:
        coerced = sorted({int(length) for length in lengths})
    except (TypeError, ValueError) as bad_length:
        raise ConfigError(f"lengths must be whole numbers, got {lengths!r}") from bad_length
    if not coerced:
        raise ConfigError("lengths was given but empty; leave it unset to predict every length found")
    if any(length < 1 for length in coerced):
        raise ConfigError(f"lengths must all be at least 1, got {sorted(coerced)}")
    return tuple(coerced)


def _coerce_strings(values) -> tuple[str, ...]:
    if values is None:
        return ()
    if isinstance(values, str):
        values = [values]
    return tuple(str(value) for value in values)


def normalise_config_keys(mapping: dict, config_class) -> dict:
    """
    Resolve a config mapping's aliases onto field names, rejecting what does not fit.

    Args:
        mapping (dict): config keys as supplied, aliases and all.
        config_class (type): the config dataclass the keys are for.

    Returns:
        dict: the same values under their field names.

    Raises:
        ConfigError: if a key is not a field of that config, or an alias and its
            own field name are both present.

    Note:
        Called at each mapping's own boundary, before anything is merged. Two
        already-normalised mappings can then be merged freely - which is what
        lets a keyword argument override a config object, and a command line
        option override a config file, without either looking like a conflict.
    """
    known = {field.name for field in fields(config_class)}
    values: dict = {}
    supplied_as: dict = {}
    for key, value in dict(mapping).items():
        name = FIELD_ALIASES.get(key, key)
        if name not in known:
            raise ConfigError(
                f"Unknown config key {key!r} for {config_class.__name__}. "
                f"Valid keys: {', '.join(sorted(known))}"
            )
        if name in values:
            # An alias and its own name in one mapping, e.g. 'allele' beside
            # 'alleles'. Whichever won would be down to dictionary order, so
            # neither does.
            raise ConfigError(
                f"Config keys {supplied_as[name]!r} and {key!r} both set {name!r}; "
                "give one or the other"
            )
        values[name] = value
        supplied_as[name] = key
    return values


class _ConfigFromMapping:
    """Shared dictionary and file loading for the config dataclasses."""

    @classmethod
    def from_dict(cls, mapping: dict):
        """
        Build a config from a dictionary, as a pipeline would hold it.

        Args:
            mapping (dict): config keys and values. Paths may be strings, and
                the singular aliases 'allele', 'length', 'netmhcpan' and 'input'
                are accepted.

        Returns:
            The config, validated.

        Raises:
            ConfigError: if a key is not a config field, or a value is unusable.
        """
        try:
            return cls(**normalise_config_keys(mapping, cls))
        except TypeError as missing:
            raise ConfigError(f"{cls.__name__}: {missing}") from missing

    @classmethod
    def from_file(cls, config_file: Path | str, section: str | None = None):
        """
        Build a config from a TOML or JSON file.

        Args:
            config_file (Path | str): a .toml or .json file of config keys.
            section (str | None): a table within the file to prefer, e.g.
                'predict'. Keys in that table override same-named keys at the
                top level, so one file can configure a whole pipeline.

        Returns:
            The config, validated.

        Raises:
            ConfigError: if the file is missing, not readable as TOML or JSON,
                or holds a key that is not a config field.
        """
        return cls.from_dict(read_config_file(config_file, section))


def read_config_file(config_file: Path | str, section: str | None = None) -> dict:
    """
    Read a TOML or JSON config file into a dictionary.

    Args:
        config_file (Path | str): a .toml or .json file.
        section (str | None): a table within the file to prefer, e.g. 'predict'.
            Keys in that table override same-named keys at the top level.

    Returns:
        dict: the config keys, with any command tables other than the one asked
            for left out.

    Raises:
        ConfigError: if the file is missing or cannot be parsed.
    """
    config_file = Path(config_file).expanduser()
    if not config_file.is_file():
        raise ConfigError(f"No config file at {config_file}")

    suffix = config_file.suffix.lower()
    try:
        if suffix == '.json':
            document = json.loads(config_file.read_text())
        elif suffix in {'.toml', ''}:
            document = tomllib.loads(config_file.read_text())
        else:
            raise ConfigError(f"Config file {config_file} must be .toml or .json, not {suffix!r}")
    except (json.JSONDecodeError, ValueError) as bad_config:
        raise ConfigError(f"Could not read {config_file}: {bad_config}") from bad_config

    if not isinstance(document, dict):
        raise ConfigError(f"Config file {config_file} must hold a table of keys")

    # Tables named after a command are the sections; everything else at the top
    # level is shared. Dropping the other commands' tables is what lets one file
    # hold both a [prepare] and a [predict] block.
    section_names = {'predict', 'prepare'}
    shared = {key: value for key, value in document.items() if key not in section_names}
    chosen = document.get(section, {}) if section else {}
    if not isinstance(chosen, dict):
        raise ConfigError(f"[{section}] in {config_file} must be a table of keys")
    return {**shared, **chosen}


@dataclass(frozen=True)
class PredictionConfig(_ConfigFromMapping):
    """
    Everything one prediction run needs.

    Attributes:
        input_folder (Path): the peptide set, laid out as the source expects.
        output_folder (Path): written as <output>/<allele_slug>/<N>mers/<dataset>.parquet.
        alleles (tuple[str, ...]): allele numbers and/or slugs, mixed freely.
        lengths (tuple[int, ...] | None): the peptide lengths to predict, or
            None for every length found in the input folder.
        mode (str): the kind of peptide collection.
        source (str | None): the input handler for that mode, or None for the
            mode's default.
        scores (str): which of netMHCpan's predictions to make and keep. 'el' is
            eluted ligand likelihood, which is netMHCpan's default; 'ba' adds
            binding affinity to it rather than replacing it.
        netmhcpan_path (Path): the netMHCpan executable, resolved at construction.
        chunk_size (int): peptides per netMHCpan process. Smaller chunks balance
            the workers better.
        workers (int): concurrent netMHCpan processes.
        peptide_column (str): the column holding peptide sequences.
        raw_output_folder (Path | None): where to write netMHCpan's own text
            output, or None to write none.
        validate_alleles (bool): check the alleles against netMHCpan's allele list.
        overwrite (bool): False skips any run whose output parquet already exists.
        strict (bool): True raises on the first failed run instead of carrying on.
    """
    input_folder: Path
    output_folder: Path
    alleles: tuple[str, ...]
    lengths: tuple[int, ...] | None = None
    mode: str = DEFAULT_MODE
    source: str | None = None
    scores: str = DEFAULT_SCORES
    netmhcpan_path: Path | None = None
    chunk_size: int = 5000
    workers: int = os.cpu_count() or 1
    peptide_column: str = 'peptide_sequence'
    raw_output_folder: Path | None = None
    validate_alleles: bool = True
    overwrite: bool = True
    strict: bool = False

    def __post_init__(self):
        _set(self, 'input_folder', Path(self.input_folder).expanduser())
        _set(self, 'output_folder', Path(self.output_folder).expanduser())
        _set(self, 'alleles', _coerce_strings(self.alleles))
        _set(self, 'lengths', _coerce_lengths(self.lengths))
        _set(self, 'netmhcpan_path', resolve_netmhcpan_path(self.netmhcpan_path))
        if self.raw_output_folder is not None:
            _set(self, 'raw_output_folder', Path(self.raw_output_folder).expanduser())

        if not self.input_folder.is_dir():
            raise ConfigError(f"No input folder at {self.input_folder}")
        if not self.alleles:
            raise ConfigError("At least one allele is needed")
        if self.chunk_size < 1:
            raise ConfigError(f"chunk_size must be at least 1, got {self.chunk_size}")
        if self.workers < 1:
            raise ConfigError(f"workers must be at least 1, got {self.workers}")
        # Resolved rather than merely checked, so a mode, source or scores typo
        # is reported here rather than after the peptides have been read.
        resolve_source(self.mode, self.source)
        get_score_set(self.scores)

    @property
    def peptide_source(self):
        """PeptideSource: the input handler this run reads its peptides with."""
        return resolve_source(self.mode, self.source)

    @property
    def score_set(self):
        """ScoreSet: which predictions to ask netMHCpan for, and how to read them."""
        return get_score_set(self.scores)

    @property
    def source_name(self) -> str:
        """str: the name of the source in use, with the mode's default filled in."""
        return self.peptide_source.name


@dataclass(frozen=True)
class PrepareConfig(_ConfigFromMapping):
    """
    Everything one preparation run needs.

    Attributes:
        input_path (Path): the file to build a peptide set from; a FASTA, for
            the proteome mode.
        output_folder (Path): written as <output>/<N>mers/<dataset>.parquet, which
            is what a prediction run then takes as its input folder.
        lengths (tuple[int, ...]): the peptide lengths to generate.
        mode (str): the kind of peptide collection to prepare.
        dataset (str | None): the name for the peptide set, or None to take the
            input file's stem. It names the output files, and travels through
            prediction onto the prediction files.
        keep_nonstandard (bool): True keeps peptides containing residues outside
            the twenty standard amino acids, which netMHCpan cannot score.
    """
    input_path: Path
    output_folder: Path
    lengths: tuple[int, ...] = (8, 9, 10, 11)
    mode: str = DEFAULT_MODE
    dataset: str | None = None
    keep_nonstandard: bool = False

    def __post_init__(self):
        _set(self, 'input_path', Path(self.input_path).expanduser())
        _set(self, 'output_folder', Path(self.output_folder).expanduser())
        _set(self, 'lengths', _coerce_lengths(self.lengths))
        _set(self, 'dataset', self.dataset or self.input_path.stem)

        if not self.input_path.is_file():
            raise ConfigError(f"No input file at {self.input_path}")
        if not self.lengths:
            raise ConfigError("At least one length is needed")
        # Raises for an unregistered mode, and for a registered one that has no
        # preparation step. Both are worth knowing before the FASTA is opened.
        get_preparer(self.mode)


def _set(config, field_name: str, value) -> None:
    # Frozen dataclasses forbid assignment, including in __post_init__, so
    # coercion has to go around it. Kept to one helper so the intent is obvious.
    object.__setattr__(config, field_name, value)
