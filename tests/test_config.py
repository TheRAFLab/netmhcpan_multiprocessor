"""Configuration: dictionaries, files, coercion, validation, and finding netMHCpan."""

import json

import pytest

from netmhcpan_multiprocessor.predict import _as_config
from netmhcpan_multiprocessor.config import (
    PredictionConfig,
    PrepareConfig,
    read_config_file,
    resolve_netmhcpan_path,
)
from netmhcpan_multiprocessor.errors import (
    ConfigError,
    NetMHCpanNotFoundError,
    UnknownModeError,
    UnknownSourceError,
)


def test_a_dictionary_is_enough(prediction_settings):
    config = PredictionConfig.from_dict(prediction_settings)
    assert config.alleles == ('HLA-A*02:01',)
    assert config.mode == 'proteome'
    assert config.source_name == 'proteome'


def test_paths_may_be_strings(prediction_settings):
    settings = {key: str(value) if 'folder' in key or 'path' in key else value
                for key, value in prediction_settings.items()}
    config = PredictionConfig.from_dict(settings)
    assert config.input_folder.is_dir()


def test_singular_aliases_are_accepted(prediction_settings):
    settings = {key: value for key, value in prediction_settings.items() if key != 'alleles'}
    config = PredictionConfig.from_dict({**settings, 'allele': 'HLA-A*02:01', 'length': 9})
    assert config.alleles == ('HLA-A*02:01',)
    assert config.lengths == (9,)


def test_lengths_are_sorted_and_deduplicated(prediction_settings):
    config = PredictionConfig.from_dict({**prediction_settings, 'lengths': [9, 8, 9]})
    assert config.lengths == (8, 9)


def test_lengths_default_to_everything_found(prediction_settings):
    assert PredictionConfig.from_dict(prediction_settings).lengths is None


@pytest.mark.parametrize('lengths', [[0], [-1], ['nine'], []])
def test_unusable_lengths_are_rejected(prediction_settings, lengths):
    with pytest.raises(ConfigError):
        PredictionConfig.from_dict({**prediction_settings, 'lengths': lengths})


def test_unknown_keys_are_rejected(prediction_settings):
    with pytest.raises(ConfigError, match='Unknown config key'):
        PredictionConfig.from_dict({**prediction_settings, 'protein_folder': 'x'})


def test_a_missing_input_folder_is_rejected(prediction_settings, tmp_path):
    with pytest.raises(ConfigError, match='No input folder'):
        PredictionConfig.from_dict({**prediction_settings, 'input_folder': tmp_path / 'absent'})


def test_alleles_are_required(prediction_settings):
    with pytest.raises(ConfigError, match='At least one allele'):
        PredictionConfig.from_dict({**prediction_settings, 'alleles': []})


@pytest.mark.parametrize('key', ['chunk_size', 'workers'])
def test_counts_must_be_positive(prediction_settings, key):
    with pytest.raises(ConfigError, match='at least 1'):
        PredictionConfig.from_dict({**prediction_settings, key: 0})


def test_an_unknown_mode_is_rejected(prediction_settings):
    with pytest.raises(UnknownModeError, match='Registered modes'):
        PredictionConfig.from_dict({**prediction_settings, 'mode': 'immunopeptidome'})


def test_an_unknown_source_is_rejected(prediction_settings):
    with pytest.raises(UnknownSourceError, match="Sources for 'proteome'"):
        PredictionConfig.from_dict({**prediction_settings, 'source': 'peptide_parquet'})


def test_keyword_arguments_are_not_smuggled_in_as_positional(prediction_settings):
    settings = {key: value for key, value in prediction_settings.items() if key != 'output_folder'}
    with pytest.raises(ConfigError, match='output_folder'):
        PredictionConfig.from_dict(settings)


class TestFindingNetMHCpan:
    def test_an_explicit_path_wins(self, fake_netmhcpan, monkeypatch):
        monkeypatch.setenv('NETMHCPAN', '/nowhere/netMHCpan')
        assert resolve_netmhcpan_path(fake_netmhcpan) == fake_netmhcpan

    def test_the_environment_is_next(self, fake_netmhcpan, monkeypatch):
        monkeypatch.setenv('NETMHCPAN', str(fake_netmhcpan))
        assert resolve_netmhcpan_path(None) == fake_netmhcpan

    def test_then_the_path(self, fake_netmhcpan, monkeypatch):
        monkeypatch.delenv('NETMHCPAN', raising=False)
        monkeypatch.setattr('shutil.which', lambda program: str(fake_netmhcpan))
        assert resolve_netmhcpan_path(None) == fake_netmhcpan

    def test_nothing_found_names_every_route(self, monkeypatch):
        monkeypatch.delenv('NETMHCPAN', raising=False)
        monkeypatch.setattr('shutil.which', lambda program: None)
        with pytest.raises(NetMHCpanNotFoundError) as raised:
            resolve_netmhcpan_path(None)
        message = str(raised.value)
        assert '--netmhcpan' in message and 'NETMHCPAN' in message and 'PATH' in message

    def test_a_configured_path_that_is_not_there(self, tmp_path):
        with pytest.raises(NetMHCpanNotFoundError, match='No netMHCpan executable'):
            resolve_netmhcpan_path(tmp_path / 'netMHCpan')

    def test_a_configured_path_that_is_not_executable(self, tmp_path):
        not_executable = tmp_path / 'netMHCpan'
        not_executable.write_text('#!/bin/sh\n')
        with pytest.raises(NetMHCpanNotFoundError, match='not executable'):
            resolve_netmhcpan_path(not_executable)

    def test_a_stale_environment_variable_is_reported(self, monkeypatch, tmp_path):
        monkeypatch.setenv('NETMHCPAN', str(tmp_path / 'gone'))
        with pytest.raises(NetMHCpanNotFoundError, match='No netMHCpan executable'):
            resolve_netmhcpan_path(None)


class TestConfigFiles:
    def test_toml(self, tmp_path, prediction_settings):
        config_file = tmp_path / 'run.toml'
        config_file.write_text(
            f"input_folder = '{prediction_settings['input_folder']}'\n"
            f"output_folder = '{prediction_settings['output_folder']}'\n"
            f"netmhcpan_path = '{prediction_settings['netmhcpan_path']}'\n"
            "alleles = ['HLA-A*02:01', 'HLA-B*07:02']\n"
            "lengths = [9]\n"
        )
        config = PredictionConfig.from_file(config_file)
        assert config.alleles == ('HLA-A*02:01', 'HLA-B*07:02')
        assert config.lengths == (9,)

    def test_json(self, tmp_path, prediction_settings):
        config_file = tmp_path / 'run.json'
        config_file.write_text(json.dumps({
            'input_folder': str(prediction_settings['input_folder']),
            'output_folder': str(prediction_settings['output_folder']),
            'netmhcpan_path': str(prediction_settings['netmhcpan_path']),
            'allele': 'HLA-A*02:01',
        }))
        assert PredictionConfig.from_file(config_file).alleles == ('HLA-A*02:01',)

    def test_a_command_table_overrides_the_shared_keys(self, tmp_path):
        config_file = tmp_path / 'pipeline.toml'
        config_file.write_text(
            "output_folder = 'shared'\n"
            "lengths = [8, 9]\n"
            "[prepare]\n"
            "output_folder = 'peptides'\n"
            "[predict]\n"
            "output_folder = 'predictions'\n"
        )
        assert read_config_file(config_file, 'predict') == {
            'output_folder': 'predictions', 'lengths': [8, 9],
        }
        assert read_config_file(config_file, 'prepare') == {
            'output_folder': 'peptides', 'lengths': [8, 9],
        }

    def test_other_commands_tables_are_left_out(self, tmp_path):
        config_file = tmp_path / 'pipeline.toml'
        config_file.write_text("[prepare]\ndataset = 'x'\n")
        assert read_config_file(config_file, 'predict') == {}

    def test_a_missing_file_is_reported(self, tmp_path):
        with pytest.raises(ConfigError, match='No config file'):
            read_config_file(tmp_path / 'absent.toml')

    def test_an_unreadable_file_is_reported(self, tmp_path):
        config_file = tmp_path / 'run.toml'
        config_file.write_text('this is not toml = = =\n')
        with pytest.raises(ConfigError, match='Could not read'):
            read_config_file(config_file)

    def test_an_unsupported_suffix_is_reported(self, tmp_path):
        config_file = tmp_path / 'run.yaml'
        config_file.write_text('input_folder: x\n')
        with pytest.raises(ConfigError, match='must be .toml or .json'):
            read_config_file(config_file)


class TestPrepareConfig:
    def test_the_dataset_defaults_to_the_input_stem(self, tmp_path):
        fasta = tmp_path / 'uniprot_UP000008992.fasta'
        fasta.write_text('>x\nMKV\n')
        config = PrepareConfig.from_dict({'input_path': fasta, 'output_folder': tmp_path / 'out'})
        assert config.dataset == 'uniprot_UP000008992'
        assert config.lengths == (8, 9, 10, 11)

    def test_the_dataset_can_be_named(self, tmp_path):
        fasta = tmp_path / 'uniprot_UP000008992.fasta'
        fasta.write_text('>x\nMKV\n')
        config = PrepareConfig.from_dict({
            'input_path': fasta, 'output_folder': tmp_path / 'out', 'dataset': 'UP000008992',
        })
        assert config.dataset == 'UP000008992'

    def test_a_missing_input_file_is_rejected(self, tmp_path):
        with pytest.raises(ConfigError, match='No input file'):
            PrepareConfig.from_dict({'input_path': tmp_path / 'gone.fasta', 'output_folder': tmp_path})

    def test_a_mode_that_cannot_prepare_is_rejected(self, tmp_path):
        fasta = tmp_path / 'x.fasta'
        fasta.write_text('>x\nMKV\n')
        with pytest.raises(UnknownModeError):
            PrepareConfig.from_dict({
                'input_path': fasta, 'output_folder': tmp_path, 'mode': 'display_screen',
            })

    def test_it_needs_no_netmhcpan(self, tmp_path, monkeypatch):
        # Preparation has nothing to do with netMHCpan, and must not require one.
        monkeypatch.delenv('NETMHCPAN', raising=False)
        monkeypatch.setattr('shutil.which', lambda program: None)
        fasta = tmp_path / 'x.fasta'
        fasta.write_text('>x\nMKV\n')
        assert PrepareConfig.from_dict({'input_path': fasta, 'output_folder': tmp_path})


@pytest.mark.parametrize('given, expected', [
    (9, (9,)),
    ('9', (9,)),
    ('11', (11,)),
    (['10', 11], (10, 11)),
])
def test_lengths_survive_being_written_as_text(prediction_settings, given, expected):
    # '11' iterated character by character would silently become the length 1,
    # and '9' would happen to work, which is how a caller learns to trust it.
    assert PredictionConfig.from_dict({**prediction_settings, 'lengths': given}).lengths == expected


def test_an_alias_and_its_own_name_together_is_an_error(prediction_settings):
    # Which one won would be down to dictionary order, so neither does.
    with pytest.raises(ConfigError, match="both set 'alleles'"):
        PredictionConfig.from_dict({**prediction_settings, 'allele': 'HLA-B*07:02'})


def test_an_alias_can_still_override_a_config_object(prediction_settings):
    # The merge sites resolve each mapping's aliases before merging, so this is
    # an override rather than the conflict above.
    config = PredictionConfig.from_dict(prediction_settings)
    assert _as_config(config, {'allele': 'HLA-B*07:02'}).alleles == ('HLA-B*07:02',)


def test_an_alias_in_a_config_file_can_be_overridden_by_an_option(prediction_settings):
    from netmhcpan_multiprocessor.config import normalise_config_keys
    from_file = normalise_config_keys({'allele': 'HLA-B*07:02'}, PredictionConfig)
    merged = {**prediction_settings, **from_file, 'alleles': ['HLA-A*02:01']}
    assert PredictionConfig.from_dict(merged).alleles == ('HLA-A*02:01',)
