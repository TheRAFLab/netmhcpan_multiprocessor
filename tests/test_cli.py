"""The command line interface."""

import json

import polars as pl
import pytest
from click.testing import CliRunner

from netmhcpan_multiprocessor.cli import cli


@pytest.fixture
def runner():
    return CliRunner()


def predict_arguments(peptide_folder, output_folder, fake_netmhcpan, *extra):
    return [
        'predict',
        '--input-folder', str(peptide_folder),
        '--output-folder', str(output_folder),
        '--netmhcpan', str(fake_netmhcpan),
        '--allele', 'HLA-A*02:01',
        *extra,
    ]


class TestPredict:
    def test_a_run_writes_its_parquet(self, runner, peptide_folder, output_folder, fake_netmhcpan):
        result = runner.invoke(cli, predict_arguments(peptide_folder, output_folder, fake_netmhcpan))
        assert result.exit_code == 0, result.output
        assert (output_folder / 'hla_a_02_01' / '9mers' / 'testset.parquet').exists()

    def test_the_summary_names_every_run(self, runner, peptide_folder, output_folder,
                                         fake_netmhcpan):
        result = runner.invoke(cli, predict_arguments(peptide_folder, output_folder, fake_netmhcpan))
        assert '8mers' in result.output and '9mers' in result.output

    def test_quiet_prints_nothing(self, runner, peptide_folder, output_folder, fake_netmhcpan):
        result = runner.invoke(
            cli, predict_arguments(peptide_folder, output_folder, fake_netmhcpan, '--quiet'))
        assert result.exit_code == 0
        assert result.output == ''

    def test_lengths_can_be_repeated(self, runner, peptide_folder, output_folder, fake_netmhcpan):
        runner.invoke(cli, predict_arguments(
            peptide_folder, output_folder, fake_netmhcpan, '--length', '9'))
        assert not (output_folder / 'hla_a_02_01' / '8mers').exists()

    def test_alleles_can_be_repeated_and_mixed(self, runner, peptide_folder, output_folder,
                                               fake_netmhcpan):
        result = runner.invoke(cli, predict_arguments(
            peptide_folder, output_folder, fake_netmhcpan, '--allele', 'hla_b_07_02'))
        assert result.exit_code == 0, result.output
        assert (output_folder / 'hla_b_07_02' / '9mers' / 'testset.parquet').exists()

    def test_affinity_can_be_asked_for(self, runner, peptide_folder, output_folder, fake_netmhcpan):
        runner.invoke(cli, predict_arguments(
            peptide_folder, output_folder, fake_netmhcpan, '--scores', 'ba'))
        predictions = pl.read_parquet(output_folder / 'hla_a_02_01' / '9mers' / 'testset.parquet')
        assert 'affinity_nm' in predictions.columns

    def test_eluted_ligand_is_the_default(self, runner, peptide_folder, output_folder,
                                          fake_netmhcpan):
        runner.invoke(cli, predict_arguments(peptide_folder, output_folder, fake_netmhcpan))
        predictions = pl.read_parquet(output_folder / 'hla_a_02_01' / '9mers' / 'testset.parquet')
        assert 'affinity_nm' not in predictions.columns

    def test_an_unknown_scores_value_is_refused_by_click(self, runner, peptide_folder,
                                                         output_folder, fake_netmhcpan):
        result = runner.invoke(cli, predict_arguments(
            peptide_folder, output_folder, fake_netmhcpan, '--scores', 'affinity'))
        assert result.exit_code != 0
        assert "'affinity'" in result.output

    def test_the_netmhcpan_path_can_come_from_the_environment(self, runner, peptide_folder,
                                                              output_folder, fake_netmhcpan,
                                                              monkeypatch):
        monkeypatch.setenv('NETMHCPAN', str(fake_netmhcpan))
        result = runner.invoke(cli, [
            'predict',
            '--input-folder', str(peptide_folder),
            '--output-folder', str(output_folder),
            '--allele', 'HLA-A*02:01',
        ])
        assert result.exit_code == 0, result.output

    def test_a_missing_netmhcpan_is_a_clean_message(self, runner, peptide_folder, output_folder,
                                                    monkeypatch):
        monkeypatch.delenv('NETMHCPAN', raising=False)
        monkeypatch.setattr('shutil.which', lambda program: None)
        result = runner.invoke(cli, [
            'predict',
            '--input-folder', str(peptide_folder),
            '--output-folder', str(output_folder),
            '--allele', 'HLA-A*02:01',
        ])
        assert result.exit_code == 1
        assert 'Could not find netMHCpan' in result.output
        # A stack trace would be the wrong way to tell someone their path is wrong.
        assert result.exception.__class__.__name__ == 'SystemExit'

    def test_an_unknown_allele_is_a_clean_message(self, runner, peptide_folder, output_folder,
                                                  fake_netmhcpan):
        result = runner.invoke(cli, [
            'predict',
            '--input-folder', str(peptide_folder),
            '--output-folder', str(output_folder),
            '--netmhcpan', str(fake_netmhcpan),
            '--allele', 'HLA-A*99:99',
        ])
        assert result.exit_code == 1
        assert 'does not know the allele' in result.output

    def test_a_failed_run_exits_non_zero(self, runner, peptide_folder, output_folder,
                                         fake_netmhcpan, monkeypatch):
        monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_MODE', 'error')
        result = runner.invoke(cli, predict_arguments(peptide_folder, output_folder, fake_netmhcpan))
        assert result.exit_code == 1
        assert 'predictions failed' in result.output

    def test_an_unknown_source_names_the_modes_sources(self, runner, peptide_folder, output_folder,
                                                       fake_netmhcpan):
        result = runner.invoke(cli, predict_arguments(
            peptide_folder, output_folder, fake_netmhcpan, '--source', 'immunopeptidome'))
        assert result.exit_code == 1
        assert "Sources for 'proteome'" in result.output

    def test_an_unknown_mode_is_refused_by_click(self, runner, peptide_folder, output_folder,
                                                 fake_netmhcpan):
        result = runner.invoke(cli, predict_arguments(
            peptide_folder, output_folder, fake_netmhcpan, '--mode', 'display_screen'))
        assert result.exit_code != 0
        assert 'display_screen' in result.output


class TestConfigFile:
    def test_a_file_can_drive_the_whole_run(self, runner, tmp_path, peptide_folder, output_folder,
                                            fake_netmhcpan):
        config_file = tmp_path / 'run.toml'
        config_file.write_text(
            f"input_folder = '{peptide_folder}'\n"
            f"output_folder = '{output_folder}'\n"
            f"netmhcpan_path = '{fake_netmhcpan}'\n"
            "alleles = ['HLA-A*02:01']\n"
            "lengths = [9]\n"
        )
        result = runner.invoke(cli, ['predict', '--config', str(config_file)])
        assert result.exit_code == 0, result.output
        assert (output_folder / 'hla_a_02_01' / '9mers' / 'testset.parquet').exists()
        assert not (output_folder / 'hla_a_02_01' / '8mers').exists()

    def test_options_win_over_the_file(self, runner, tmp_path, peptide_folder, output_folder,
                                       fake_netmhcpan):
        config_file = tmp_path / 'run.json'
        config_file.write_text(json.dumps({
            'input_folder': str(peptide_folder),
            'output_folder': str(tmp_path / 'ignored'),
            'netmhcpan_path': str(fake_netmhcpan),
            'alleles': ['HLA-A*02:01'],
        }))
        result = runner.invoke(cli, [
            'predict', '--config', str(config_file), '--output-folder', str(output_folder),
        ])
        assert result.exit_code == 0, result.output
        assert (output_folder / 'hla_a_02_01' / '9mers' / 'testset.parquet').exists()
        assert not (tmp_path / 'ignored').exists()

    def test_a_file_can_set_something_with_a_default(self, runner, tmp_path, peptide_folder,
                                                     output_folder, fake_netmhcpan):
        # Testing for None would not do here: chunk_size has a non-None default,
        # so only click's own parameter source can tell 'unset' from 'set to 5000'.
        config_file = tmp_path / 'run.toml'
        config_file.write_text(
            f"input_folder = '{peptide_folder}'\n"
            f"output_folder = '{output_folder}'\n"
            f"netmhcpan_path = '{fake_netmhcpan}'\n"
            "alleles = ['HLA-A*02:01']\n"
            "chunk_size = 7\n"
        )
        result = runner.invoke(cli, ['predict', '--config', str(config_file)])
        assert result.exit_code == 0, result.output
        # 30 peptides in chunks of 7 is 5 chunks, rather than the default's 1.
        assert '5' in result.output

    def test_one_file_can_hold_both_commands(self, runner, tmp_path, peptide_folder, output_folder,
                                             fake_netmhcpan):
        config_file = tmp_path / 'pipeline.toml'
        config_file.write_text(
            "[prepare]\n"
            f"input_path = '{tmp_path / 'absent.fasta'}'\n"
            "[predict]\n"
            f"input_folder = '{peptide_folder}'\n"
            f"output_folder = '{output_folder}'\n"
            f"netmhcpan_path = '{fake_netmhcpan}'\n"
            "alleles = ['HLA-A*02:01']\n"
        )
        result = runner.invoke(cli, ['predict', '--config', str(config_file)])
        assert result.exit_code == 0, result.output


class TestPrepare:
    def test_it_writes_a_peptide_set_prediction_can_read(self, runner, tmp_path, fake_netmhcpan):
        pytest.importorskip('polars_bio', reason="preparation needs the 'proteome' extra")
        fasta = tmp_path / 'testset.fasta'
        fasta.write_text('>sp|P00001|P_TEST Test OS=Testus testus OX=1\nMESRIWCLVVCVNLTALWTL\n')
        peptides = tmp_path / 'peptides'

        prepared = runner.invoke(cli, [
            'prepare', '--input-path', str(fasta), '--output-folder', str(peptides),
            '--length', '9',
        ])
        assert prepared.exit_code == 0, prepared.output
        assert (peptides / '9mers' / 'testset.parquet').exists()

        # The two halves of the package, joined only by the folder layout.
        predicted = runner.invoke(cli, [
            'predict', '--input-folder', str(peptides),
            '--output-folder', str(tmp_path / 'predictions'),
            '--netmhcpan', str(fake_netmhcpan), '--allele', 'HLA-A*02:01',
        ])
        assert predicted.exit_code == 0, predicted.output
        assert (tmp_path / 'predictions' / 'hla_a_02_01' / '9mers' / 'testset.parquet').exists()

    def test_a_missing_input_is_a_clean_message(self, runner, tmp_path):
        result = runner.invoke(cli, [
            'prepare', '--input-path', str(tmp_path / 'gone.fasta'),
            '--output-folder', str(tmp_path / 'out'),
        ])
        assert result.exit_code == 1
        assert 'No input file' in result.output

    def test_it_needs_no_netmhcpan(self, runner, tmp_path, monkeypatch):
        pytest.importorskip('polars_bio', reason="preparation needs the 'proteome' extra")
        monkeypatch.delenv('NETMHCPAN', raising=False)
        monkeypatch.setattr('shutil.which', lambda program: None)
        fasta = tmp_path / 'testset.fasta'
        fasta.write_text('>sp|P00001|P_TEST Test OS=Testus testus OX=1\nMESRIWCLVVCVNLTALWTL\n')
        result = runner.invoke(cli, [
            'prepare', '--input-path', str(fasta), '--output-folder', str(tmp_path / 'out'),
            '--length', '9',
        ])
        assert result.exit_code == 0, result.output


def test_the_help_lists_the_modes_sources(runner):
    result = runner.invoke(cli, ['predict', '--help'])
    assert 'proteome: proteome' in result.output


def test_the_help_explains_the_score_sets(runner):
    result = runner.invoke(cli, ['predict', '--help'])
    assert 'eluted ligand' in result.output
    assert 'binding affinity' in result.output


def test_a_config_files_netmhcpan_beats_the_environment(runner, tmp_path, peptide_folder,
                                                        output_folder, fake_netmhcpan,
                                                        monkeypatch):
    # The environment is a fallback, not an option someone typed, so it must not
    # win over --config the way a real command line option does.
    monkeypatch.setenv('NETMHCPAN', str(tmp_path / 'stale' / 'netMHCpan'))
    config_file = tmp_path / 'run.toml'
    config_file.write_text(
        f"input_folder = '{peptide_folder}'\n"
        f"output_folder = '{output_folder}'\n"
        f"netmhcpan_path = '{fake_netmhcpan}'\n"
        "alleles = ['HLA-A*02:01']\n"
    )
    result = runner.invoke(cli, ['predict', '--config', str(config_file)])
    assert result.exit_code == 0, result.output


def test_an_option_still_beats_the_environment(runner, peptide_folder, output_folder,
                                               fake_netmhcpan, monkeypatch, tmp_path):
    monkeypatch.setenv('NETMHCPAN', str(tmp_path / 'stale' / 'netMHCpan'))
    result = runner.invoke(cli, predict_arguments(peptide_folder, output_folder, fake_netmhcpan))
    assert result.exit_code == 0, result.output


class TestVerbose:
    def test_it_shows_the_resolved_settings(self, runner, peptide_folder, output_folder,
                                            fake_netmhcpan):
        result = runner.invoke(cli, predict_arguments(
            peptide_folder, output_folder, fake_netmhcpan, '--verbose'))
        assert result.exit_code == 0, result.output
        # The values worth checking are the ones nobody typed: defaults, and
        # anything a config file set.
        assert 'chunk_size' in result.output
        assert 'peptide_sequence' in result.output
        assert 'validate_alleles' in result.output

    def test_it_shows_what_each_allele_resolved_to(self, runner, peptide_folder, output_folder,
                                                   fake_netmhcpan):
        result = runner.invoke(cli, predict_arguments(
            peptide_folder, output_folder, fake_netmhcpan, '--verbose', '--allele', 'hla_b_07_02'))
        assert 'HLA-B07:02' in result.output      # what netMHCpan is asked for
        assert 'hla_b_07_02' in result.output     # what the folder is called

    def test_it_shows_the_input_files_it_found(self, runner, peptide_folder, output_folder,
                                               fake_netmhcpan):
        result = runner.invoke(cli, predict_arguments(
            peptide_folder, output_folder, fake_netmhcpan, '--verbose'))
        assert 'Peptide files' in result.output
        assert 'testset' in result.output

    def test_it_shows_what_netmhcpan_wrote_to_stderr(self, runner, peptide_folder, output_folder,
                                                     fake_netmhcpan, monkeypatch):
        monkeypatch.setenv('FAKE_NETMHCPAN_STDERR_BYTES', '200')
        result = runner.invoke(cli, predict_arguments(
            peptide_folder, output_folder, fake_netmhcpan, '--verbose'))
        assert result.exit_code == 0, result.output
        assert 'netMHCpan said' in result.output
        assert 'grumble' in result.output

    def test_it_says_so_when_netmhcpan_was_quiet(self, runner, peptide_folder, output_folder,
                                                 fake_netmhcpan):
        result = runner.invoke(cli, predict_arguments(
            peptide_folder, output_folder, fake_netmhcpan, '--verbose'))
        assert 'wrote nothing to stderr' in result.output

    def test_it_stops_before_running_if_the_input_is_wrong(self, runner, tmp_path, output_folder,
                                                           fake_netmhcpan):
        empty = tmp_path / 'empty'
        empty.mkdir()
        result = runner.invoke(cli, predict_arguments(
            empty, output_folder, fake_netmhcpan, '--verbose'))
        assert result.exit_code == 1
        assert 'No peptide parquet files' in result.output

    def test_quiet_and_verbose_together_is_refused(self, runner, peptide_folder, output_folder,
                                                   fake_netmhcpan):
        result = runner.invoke(cli, predict_arguments(
            peptide_folder, output_folder, fake_netmhcpan, '--verbose', '--quiet'))
        assert result.exit_code != 0
        assert 'opposite things' in result.output

    def test_prepare_shows_its_settings(self, runner, tmp_path):
        pytest.importorskip('polars_bio', reason="preparation needs the 'proteome' extra")
        fasta = tmp_path / 'testset.fasta'
        fasta.write_text('>sp|P00001|P_TEST Test OS=Testus testus OX=1\nMESRIWCLVVCVNLTALWTL\n')
        result = runner.invoke(cli, [
            'prepare', '--input-path', str(fasta), '--output-folder', str(tmp_path / 'out'),
            '--length', '9', '--verbose',
        ])
        assert result.exit_code == 0, result.output
        assert 'keep_nonstandard' in result.output
