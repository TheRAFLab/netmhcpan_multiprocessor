"""
netmhcpan_multiprocessor: run netMHCpan over peptide parquet sets, in parallel.

Peptide parquet in, prediction parquet out, one file per length per allele, with
the input's row order preserved so the two join directly.

    from netmhcpan_multiprocessor import predict

    result = predict({
        'input_folder': 'peptides',
        'output_folder': 'predictions',
        'alleles': ['HLA-A*02:01', 'hla_a_01_01'],
    })

The same settings drive the `netmhcpan-mp` command, which is a thin layer over
this function. Preparing a peptide set is a separate step, `prepare`, so that
prediction's input is always a set of peptide parquet files however they
were made.
"""

from netmhcpan_multiprocessor.alleles import ResolvedAllele, resolve_alleles, slugify
from netmhcpan_multiprocessor.config import PredictionConfig, PrepareConfig, normalise_config_keys
from netmhcpan_multiprocessor.errors import (
    ConfigError,
    NetMHCpanMultiprocessorError,
    NetMHCpanNotFoundError,
    NetMHCpanRunError,
    PeptideInputError,
    PredictionFailedError,
    PreparationNotSupportedError,
    UnknownAlleleError,
    UnknownModeError,
    UnknownSourceError,
)
from netmhcpan_multiprocessor.modes import (
    DEFAULT_MODE,
    PeptideFile,
    PrepareResult,
    get_mode,
    get_preparer,
    registered_modes,
    registered_sources,
    resolve_source,
)
from netmhcpan_multiprocessor.parsing import PREDICTION_COLUMNS
from netmhcpan_multiprocessor.predict import PredictionResult, RunResult, predict
from netmhcpan_multiprocessor.progress import NullProgressReporter, ProgressReporter

__version__ = '0.2.0'


def prepare(config: PrepareConfig | dict | None = None, **settings) -> PrepareResult:
    """
    Build a peptide set that a prediction run can then read.

    Args:
        config (PrepareConfig | dict | None): the settings, as a config object
            or a dictionary.
        **settings: config keys given as keyword arguments, merged over config.

    Returns:
        PrepareResult: the peptide files written, one per length.

    Raises:
        ConfigError: if the settings are unusable, or the mode has no
            preparation step.
        PeptideInputError: if the input yields no peptides.
    """
    if isinstance(config, PrepareConfig) and not settings:
        prepare_config = config
    else:
        # Aliases resolved per mapping, before the merge, so a keyword argument
        # overrides the same field on the config rather than colliding with it.
        settings = normalise_config_keys(settings, PrepareConfig)
        if isinstance(config, PrepareConfig):
            base = {name: getattr(config, name) for name in config.__dataclass_fields__}
            prepare_config = PrepareConfig.from_dict({**base, **settings})
        else:
            prepare_config = PrepareConfig.from_dict(
                {**normalise_config_keys(config or {}, PrepareConfig), **settings})
    return get_preparer(prepare_config.mode).prepare(prepare_config)


__all__ = [
    'ConfigError',
    'DEFAULT_MODE',
    'NetMHCpanMultiprocessorError',
    'NetMHCpanNotFoundError',
    'NetMHCpanRunError',
    'NullProgressReporter',
    'PREDICTION_COLUMNS',
    'PeptideFile',
    'PeptideInputError',
    'PredictionConfig',
    'PredictionFailedError',
    'PredictionResult',
    'PrepareConfig',
    'PrepareResult',
    'PreparationNotSupportedError',
    'ProgressReporter',
    'ResolvedAllele',
    'RunResult',
    'UnknownAlleleError',
    'UnknownModeError',
    'UnknownSourceError',
    '__version__',
    'get_mode',
    'get_preparer',
    'predict',
    'prepare',
    'registered_modes',
    'registered_sources',
    'resolve_alleles',
    'resolve_source',
    'slugify',
]
