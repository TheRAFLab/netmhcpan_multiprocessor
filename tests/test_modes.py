"""The mode and source registry."""

import dataclasses
import subprocess
import sys

import pytest

from netmhcpan_multiprocessor import modes
from netmhcpan_multiprocessor.errors import (
    PreparationNotSupportedError,
    UnknownModeError,
    UnknownSourceError,
)


def test_the_proteome_mode_is_registered():
    assert modes.registered_modes() == ('proteome',)
    assert modes.DEFAULT_MODE == 'proteome'


def test_a_mode_can_be_looked_up_by_name_or_object():
    mode = modes.get_mode('proteome')
    assert modes.get_mode(mode) is mode


def test_an_unknown_mode_lists_the_registered_ones():
    with pytest.raises(UnknownModeError, match='proteome'):
        modes.get_mode('immunopeptidome')


def test_the_source_defaults_to_the_modes_own():
    assert modes.resolve_source('proteome').name == 'proteome'
    assert modes.resolve_source('proteome', None).name == 'proteome'


def test_a_named_source_resolves():
    assert modes.resolve_source('proteome', 'proteome').name == 'proteome'


def test_an_unknown_source_lists_that_modes_sources():
    with pytest.raises(UnknownSourceError) as raised:
        modes.resolve_source('proteome', 'immunopeptidome')
    message = str(raised.value)
    # The message has to be about this mode's sources, so naming a source that
    # belongs to another mode says so rather than just 'not found'.
    assert "Mode 'proteome' has no source 'immunopeptidome'" in message
    assert "Sources for 'proteome': proteome" in message


def test_registered_sources():
    assert modes.registered_sources('proteome') == ('proteome',)


def test_the_proteome_mode_can_prepare():
    assert modes.get_preparer('proteome').name == 'proteome'


def test_a_mode_without_a_preparer_says_so():
    # Preparation is expected to leave this package one day, and the future
    # modes will have a source before they have a preparer, so 'no preparer' has
    # to be a supported state rather than an assertion failure.
    sourceless = dataclasses.replace(modes.get_mode('proteome'), name='display_screen', preparer=None)
    with pytest.raises(PreparationNotSupportedError, match='has no preparation step'):
        modes.get_preparer(sourceless)


def test_predicting_never_imports_the_preparation_module():
    # The fence around proteome preparation: predict must not drag in polars-bio
    # or shared-helpers, which are an optional extra. Run in a subprocess because
    # the only honest way to ask what an import pulls in is a fresh interpreter -
    # and because emptying sys.modules in this one would leave the rest of the
    # suite holding two copies of every class in the package.
    check = subprocess.run(
        [sys.executable, '-c', (
            'import netmhcpan_multiprocessor, sys;'
            'print("netmhcpan_multiprocessor.modes.proteome.prepare" in sys.modules)'
        )],
        capture_output=True,
        text=True,
        check=True,
    )
    assert check.stdout.strip() == 'False'
