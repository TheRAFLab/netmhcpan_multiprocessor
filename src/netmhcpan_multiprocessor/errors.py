"""
The package's error hierarchy.

Everything raised deliberately derives from NetMHCpanMultiprocessorError, so a
pipeline can catch the whole package with one except clause, while the CLI can
still turn each kind into its own message.
"""


class NetMHCpanMultiprocessorError(Exception):
    """Base class for every error this package raises deliberately."""


class ConfigError(NetMHCpanMultiprocessorError):
    """A configuration value is missing, unknown or unusable."""


class NetMHCpanNotFoundError(ConfigError):
    """No netMHCpan executable could be resolved."""


class UnknownModeError(ConfigError):
    """The requested peptide set mode is not registered."""


class UnknownSourceError(ConfigError):
    """The requested source is not registered for the requested mode."""


class PreparationNotSupportedError(ConfigError):
    """The requested mode has no preparation step."""


class UnknownAlleleError(ConfigError):
    """An allele could not be resolved to a name netMHCpan knows."""


class PeptideInputError(NetMHCpanMultiprocessorError):
    """The peptide input is missing, malformed, or not what it claims to be."""


class NetMHCpanRunError(NetMHCpanMultiprocessorError):
    """A netMHCpan process failed, or returned less than it was asked for."""


class PredictionFailedError(NetMHCpanMultiprocessorError):
    """One or more prediction runs failed under strict=True."""
