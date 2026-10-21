"""Prediction end to end, against the stand-in netMHCpan."""

import importlib

import polars as pl
import pytest

from netmhcpan_multiprocessor import predict
from netmhcpan_multiprocessor.config import PredictionConfig
from netmhcpan_multiprocessor.errors import (
    ConfigError,
    PeptideInputError,
    PredictionFailedError,
    UnknownAlleleError,
)

from conftest import make_peptides, write_peptide_parquet


def read_predictions(output_folder, allele_slug, length, dataset='testset'):
    return pl.read_parquet(output_folder / allele_slug / f"{length}mers" / f"{dataset}.parquet")


class TestOutput:
    def test_one_parquet_per_length_per_allele(self, prediction_settings, output_folder):
        result = predict({**prediction_settings, 'alleles': ['HLA-A*02:01', 'HLA-B*07:02']})
        assert result.succeeded
        written = sorted(
            path.relative_to(output_folder).as_posix()
            for path in output_folder.rglob('*.parquet')
        )
        assert written == [
            'hla_a_02_01/8mers/testset.parquet',
            'hla_a_02_01/9mers/testset.parquet',
            'hla_b_07_02/8mers/testset.parquet',
            'hla_b_07_02/9mers/testset.parquet',
        ]

    def test_the_output_holds_one_row_per_input_peptide(self, prediction_settings, output_folder):
        predict(prediction_settings)
        assert read_predictions(output_folder, 'hla_a_02_01', 9).height == 30

    def test_the_output_keeps_the_input_row_order(self, prediction_settings, output_folder):
        # The reason no identity columns are carried through: order is preserved,
        # so the two frames join positionally.
        predict(prediction_settings)
        predictions = read_predictions(output_folder, 'hla_a_02_01', 9)
        assert predictions['peptide_sequence'].to_list() == make_peptides(9, 30)

    def test_order_survives_chunks_finishing_out_of_order(self, prediction_settings, output_folder):
        # Small chunks and many workers so the chunks genuinely race each other.
        predict({**prediction_settings, 'chunk_size': 3, 'workers': 8})
        predictions = read_predictions(output_folder, 'hla_a_02_01', 9)
        assert predictions['peptide_sequence'].to_list() == make_peptides(9, 30)

    def test_the_allele_slug_column_matches_its_folder(self, prediction_settings, output_folder):
        predict({**prediction_settings, 'alleles': ['hla_b_07_02']})
        predictions = read_predictions(output_folder, 'hla_b_07_02', 9)
        assert predictions['allele_slug'].unique().to_list() == ['hla_b_07_02']

    def test_the_output_joins_back_to_the_peptide_parquet(self, prediction_settings, peptide_folder,
                                                          output_folder):
        predict(prediction_settings)
        peptides = pl.read_parquet(peptide_folder / '9mers' / 'testset.parquet')
        predictions = read_predictions(output_folder, 'hla_a_02_01', 9)
        joined = peptides.join(predictions, on='peptide_sequence', how='inner')
        assert joined.height == peptides.height

    def test_binder_counts_are_reported(self, prediction_settings):
        result = predict(prediction_settings)
        nine_mers = next(run for run in result.runs if run.length == 9)
        expected_strong = sum(1 for peptide in make_peptides(9, 30) if peptide.startswith('S'))
        expected_weak = sum(1 for peptide in make_peptides(9, 30) if peptide.startswith('W'))
        assert (nine_mers.strong_binders, nine_mers.weak_binders) == (expected_strong, expected_weak)
        assert expected_strong and expected_weak  # the fixture must exercise both

    def test_the_summary_is_available_as_a_frame(self, prediction_settings):
        summary = predict(prediction_settings).frame()
        assert summary.height == 2
        assert summary['peptide_count'].sum() == 60


class TestAlleleForms:
    @pytest.mark.parametrize('allele', ['HLA-A*02:01', 'HLA-A02:01', 'hla_a_02_01'])
    def test_any_spelling_writes_the_same_folder(self, prediction_settings, output_folder, allele):
        predict({**prediction_settings, 'alleles': [allele]})
        assert (output_folder / 'hla_a_02_01' / '9mers' / 'testset.parquet').exists()

    def test_numbers_and_slugs_can_be_mixed(self, prediction_settings, output_folder):
        result = predict({**prediction_settings, 'alleles': ['HLA-A*02:01', 'hla_b_07_02']})
        assert {run.allele for run in result.runs} == {'HLA-A*02:01', 'HLA-B*07:02'}

    def test_the_same_allele_twice_is_one_run(self, prediction_settings):
        result = predict({**prediction_settings, 'alleles': ['HLA-A*02:01', 'hla_a_02_01']})
        assert len(result.runs) == 2  # two lengths, one allele

    def test_an_unknown_allele_stops_before_any_work(self, prediction_settings, output_folder):
        with pytest.raises(UnknownAlleleError):
            predict({**prediction_settings, 'alleles': ['HLA-A*99:99']})
        assert not output_folder.exists()


class TestSelection:
    def test_lengths_can_be_restricted(self, prediction_settings, output_folder):
        predict({**prediction_settings, 'lengths': [9]})
        assert not (output_folder / 'hla_a_02_01' / '8mers').exists()

    def test_several_datasets_are_predicted_separately(self, prediction_settings, peptide_folder,
                                                       output_folder):
        write_peptide_parquet(peptide_folder, length=9, count=10, dataset='second')
        result = predict({**prediction_settings, 'lengths': [9]})
        assert {run.dataset for run in result.runs} == {'testset', 'second'}
        assert (output_folder / 'hla_a_02_01' / '9mers' / 'second.parquet').exists()

    def test_no_overwrite_skips_what_is_already_there(self, prediction_settings, output_folder):
        predict(prediction_settings)
        stamps = {path: path.stat().st_mtime_ns for path in output_folder.rglob('*.parquet')}

        result = predict({**prediction_settings, 'overwrite': False})
        assert all(run.skipped for run in result.runs)
        assert all(path.stat().st_mtime_ns == stamp for path, stamp in stamps.items())

    def test_overwrite_rewrites(self, prediction_settings, output_folder):
        predict(prediction_settings)
        result = predict(prediction_settings)
        assert not any(run.skipped for run in result.runs)


class TestRawOutput:
    def test_no_raw_text_by_default(self, prediction_settings, tmp_path):
        predict(prediction_settings)
        assert not list(tmp_path.rglob('*.txt'))

    def test_raw_text_when_asked_for(self, prediction_settings, tmp_path):
        raw_folder = tmp_path / 'raw'
        result = predict({**prediction_settings, 'raw_output_folder': raw_folder})
        raw_file = raw_folder / 'testset__9mers_hla_a_02_01.txt'
        assert raw_file.exists()
        assert len(raw_file.read_text().splitlines()) == 30
        assert next(run for run in result.runs if run.length == 9).raw_path == raw_file

    def test_raw_positions_are_renumbered_across_chunks(self, prediction_settings, tmp_path):
        raw_folder = tmp_path / 'raw'
        predict({**prediction_settings, 'raw_output_folder': raw_folder, 'chunk_size': 7})
        raw_file = raw_folder / 'testset__9mers_hla_a_02_01.txt'
        positions = [int(line.split()[0]) for line in raw_file.read_text().splitlines()]
        assert positions == list(range(1, 31))


class TestFailures:
    def test_one_bad_allele_does_not_lose_the_others(self, prediction_settings, output_folder,
                                                     monkeypatch):
        monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_MODE', 'error')
        monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_ALLELE', 'HLA-B07:02')

        result = predict({**prediction_settings, 'alleles': ['HLA-A*02:01', 'HLA-B*07:02']})

        assert [run.allele for run in result.failures] == ['HLA-B*07:02'] * 2
        assert (output_folder / 'hla_a_02_01' / '9mers' / 'testset.parquet').exists()
        assert not (output_folder / 'hla_b_07_02').exists()

    def test_an_error_printed_with_a_zero_exit_is_still_a_failure(self, prediction_settings,
                                                                  monkeypatch):
        # netMHCpan reports an unusable allele on stdout and then exits 0, so the
        # return code alone never sees it.
        monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_MODE', 'error')
        result = predict(prediction_settings)
        assert not result.succeeded
        assert 'not found' in result.failures[0].error

    def test_a_non_zero_exit_is_a_failure(self, prediction_settings, monkeypatch):
        monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_MODE', 'exit')
        result = predict(prediction_settings)
        assert all('exit 1' in run.error for run in result.runs)

    def test_a_short_result_is_a_failure(self, prediction_settings, monkeypatch):
        # The backstop for any way netMHCpan quietly returns less than asked:
        # one prediction per peptide is the whole contract of -p mode.
        monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_MODE', 'truncate')
        result = predict({**prediction_settings, 'chunk_size': 30})
        assert all('expected 30 predictions, got 29' in run.error for run in result.runs)

    def test_a_failed_run_writes_nothing(self, prediction_settings, output_folder, monkeypatch):
        monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_MODE', 'error')
        predict(prediction_settings)
        assert not list(output_folder.rglob('*.parquet'))

    def test_strict_raises(self, prediction_settings, monkeypatch):
        monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_MODE', 'error')
        with pytest.raises(PredictionFailedError, match='predictions failed'):
            predict({**prediction_settings, 'strict': True})

    def test_strict_reports_every_run(self, prediction_settings, monkeypatch):
        # Runs cancelled by the strict stop still have to appear, rather than
        # vanishing from the summary as if they had never been planned.
        monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_MODE', 'error')
        monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_ALLELE', 'HLA-A02:01')
        with pytest.raises(PredictionFailedError) as raised:
            predict({**prediction_settings, 'alleles': ['HLA-A*02:01'], 'strict': True})
        assert 'not found' in str(raised.value)

    def test_a_bad_peptide_file_stops_the_run(self, prediction_settings, peptide_folder):
        pl.DataFrame({'peptide_sequence': ['TOOSHORT']}).write_parquet(
            peptide_folder / '9mers' / 'testset.parquet'
        )
        with pytest.raises(PeptideInputError, match='length 8'):
            predict(prediction_settings)


class TestCallingConventions:
    def test_a_config_object(self, prediction_settings):
        assert predict(PredictionConfig.from_dict(prediction_settings)).succeeded

    def test_keyword_arguments(self, prediction_settings):
        settings = dict(prediction_settings)
        alleles = settings.pop('alleles')
        assert predict(settings, alleles=alleles).succeeded

    def test_keyword_arguments_override_a_config_object(self, prediction_settings, output_folder):
        config = PredictionConfig.from_dict(prediction_settings)
        predict(config, lengths=[9])
        assert not (output_folder / 'hla_a_02_01' / '8mers').exists()

    def test_progress_is_reported(self, prediction_settings):
        class RecordingReporter:
            def __init__(self):
                self.runs, self.advanced, self.failed = [], 0, []

            def add_run(self, label, total):
                self.runs.append((label, total))
                return len(self.runs)

            def advance(self, handle, amount=1):
                self.advanced += amount

            def fail(self, handle, message):
                self.failed.append(message)

        reporter = RecordingReporter()
        predict(prediction_settings, progress=reporter)
        assert sorted(reporter.runs) == [('8mers HLA-A*02:01', 30), ('9mers HLA-A*02:01', 30)]
        assert reporter.advanced == 60
        assert reporter.failed == []

    def test_nothing_is_printed_by_default(self, prediction_settings, capsys):
        # A pipeline importing this package should get no terminal output.
        predict(prediction_settings)
        assert capsys.readouterr().out == ''


class TestScoreSets:
    def test_eluted_ligand_by_default(self, prediction_settings, output_folder):
        predict(prediction_settings)
        predictions = read_predictions(output_folder, 'hla_a_02_01', 9)
        assert 'score_el' in predictions.columns
        assert 'affinity_nm' not in predictions.columns

    def test_affinity_adds_its_columns(self, prediction_settings, output_folder):
        predict({**prediction_settings, 'scores': 'ba'})
        predictions = read_predictions(output_folder, 'hla_a_02_01', 9)
        # -BA adds affinity to the eluted ligand output rather than replacing it.
        assert predictions['score_el'].null_count() == 0
        assert predictions['affinity_nm'].null_count() == 0
        assert predictions['rank_ba'].null_count() == 0

    def test_affinity_keeps_the_row_order_and_count(self, prediction_settings, output_folder):
        predict({**prediction_settings, 'scores': 'ba', 'chunk_size': 4})
        predictions = read_predictions(output_folder, 'hla_a_02_01', 9)
        assert predictions['peptide_sequence'].to_list() == make_peptides(9, 30)

    def test_binder_levels_still_parse_under_affinity(self, prediction_settings):
        result = predict({**prediction_settings, 'scores': 'ba'})
        nine_mers = next(run for run in result.runs if run.length == 9)
        assert nine_mers.strong_binders and nine_mers.weak_binders

    def test_an_unknown_score_set_is_rejected(self, prediction_settings):
        with pytest.raises(ConfigError, match='Unknown scores'):
            predict({**prediction_settings, 'scores': 'affinity'})


class TestFailureIsolationIsThorough:
    def test_an_unexpected_exception_does_not_lose_the_other_runs(self, prediction_settings,
                                                                  output_folder, monkeypatch):
        # Not every way a chunk can fail is a netMHCpan error this package knows
        # about. Whatever the cause, the promise when strict is off is that the
        # other runs still land.
        # importlib, because the package exports a `predict` function that
        # shadows the module of the same name as an attribute.
        predict_module = importlib.import_module('netmhcpan_multiprocessor.predict')

        real_run_chunk = predict_module.run_chunk

        def sometimes_broken(*arguments, alleles, **keywords):
            if 'HLA-B07:02' in alleles:
                raise RuntimeError('something nobody anticipated')
            return real_run_chunk(*arguments, alleles=alleles, **keywords)

        monkeypatch.setattr(predict_module, 'run_chunk', sometimes_broken)

        result = predict({**prediction_settings, 'alleles': ['HLA-A*02:01', 'HLA-B*07:02']})
        assert {run.allele for run in result.failures} == {'HLA-B*07:02'}
        assert 'RuntimeError: something nobody anticipated' in result.failures[0].error
        # The other allele knew nothing about any of this and must have written.
        assert (output_folder / 'hla_a_02_01' / '9mers' / 'testset.parquet').exists()

    def test_a_reporter_that_raises_cannot_fail_a_run(self, prediction_settings, output_folder):
        # Progress is cosmetic and the reporter belongs to the caller. A logger
        # or queue that falls over must not lose a long run's results.
        class BrokenReporter:
            def add_run(self, label, total):
                raise RuntimeError('reporter fell over on add_run')

            def advance(self, handle, amount=1):
                raise RuntimeError('reporter fell over on advance')

            def fail(self, handle, message):
                raise RuntimeError('reporter fell over on fail')

        result = predict(prediction_settings, progress=BrokenReporter())
        assert result.succeeded, [run.error for run in result.runs]
        assert (output_folder / 'hla_a_02_01' / '9mers' / 'testset.parquet').exists()

    def test_a_failed_run_does_not_leave_the_previous_output_behind(self, prediction_settings,
                                                                    output_folder, monkeypatch):
        # Otherwise a rerun with overwrite off would skip the path as though it
        # held this run's results, and the documented positional join would be
        # silently joining the previous peptide set's predictions.
        predict(prediction_settings)
        parquet = output_folder / 'hla_a_02_01' / '9mers' / 'testset.parquet'
        assert parquet.exists()

        monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_MODE', 'error')
        result = predict(prediction_settings)
        assert not result.succeeded
        assert not parquet.exists()

    def test_a_failed_run_clears_its_raw_output_too(self, prediction_settings, tmp_path,
                                                    monkeypatch):
        raw_folder = tmp_path / 'raw'
        predict({**prediction_settings, 'raw_output_folder': raw_folder})
        raw_file = raw_folder / 'testset__9mers_hla_a_02_01.txt'
        assert raw_file.exists()

        monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_MODE', 'error')
        predict({**prediction_settings, 'raw_output_folder': raw_folder})
        assert not raw_file.exists()

    def test_no_partial_parquet_is_left_behind(self, prediction_settings, output_folder):
        # The parquet is written beside its target and moved into place, so an
        # interrupted write can never be read as a short peptide set.
        predict(prediction_settings)
        assert not list(output_folder.rglob('*.partial'))

    def test_every_run_is_reported_even_when_strict_stops_the_pool(self, prediction_settings,
                                                                   monkeypatch):
        # A run abandoned by the strict stop still has to appear in the count,
        # rather than vanishing as though it had never been planned.
        monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_MODE', 'error')
        settings = {**prediction_settings, 'alleles': ['HLA-A*02:01', 'HLA-B*07:02'],
                    'strict': True, 'workers': 1}
        with pytest.raises(PredictionFailedError) as raised:
            predict(settings)
        # Two lengths by two alleles, however few of them got to run.
        assert 'of 4 predictions failed' in str(raised.value)

    def test_netmhcpan_writing_a_lot_to_stderr_does_not_hang(self, prediction_settings,
                                                             monkeypatch):
        # stdout and stderr are both pipes. Reading one to EOF before touching
        # the other deadlocks as soon as netMHCpan writes more than the stderr
        # pipe buffer holds, and a hung chunk hangs the whole run with it.
        monkeypatch.setenv('FAKE_NETMHCPAN_STDERR_BYTES', str(512 * 1024))
        result = predict(prediction_settings)
        assert result.succeeded, [run.error for run in result.runs]

    def test_out_of_order_predictions_are_caught(self, prediction_settings, monkeypatch):
        # netMHCpan answering in a different order than it was asked would be
        # self-consistent - the peptide column comes from its own output - and
        # would silently break the positional join, so it is checked explicitly.
        monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_MODE', 'shuffle')
        result = predict({**prediction_settings, 'chunk_size': 30})
        assert not result.succeeded
        assert 'different order' in result.failures[0].error


class TestWarnings:
    def test_netmhcpan_warnings_are_kept_on_a_successful_run(self, prediction_settings,
                                                             monkeypatch):
        # Discarded on success until now, which is the one case where you have
        # results in hand and no way to see what it grumbled about.
        monkeypatch.setenv('FAKE_NETMHCPAN_STDERR_BYTES', '200')
        result = predict(prediction_settings)
        assert result.succeeded
        assert all('grumble' in ' '.join(run.warnings) for run in result.runs)

    def test_repeated_warnings_are_reported_once(self, prediction_settings, monkeypatch):
        # Every chunk repeats the same complaint; the point is what was said.
        monkeypatch.setenv('FAKE_NETMHCPAN_STDERR_BYTES', '4000')
        result = predict({**prediction_settings, 'chunk_size': 3})
        nine_mers = next(run for run in result.runs if run.length == 9)
        assert nine_mers.chunk_count == 10
        assert len(nine_mers.warnings) == 1

    def test_a_quiet_netmhcpan_leaves_no_warnings(self, prediction_settings):
        assert all(run.warnings == [] for run in predict(prediction_settings).runs)

    def test_a_failed_run_puts_it_in_the_error_instead(self, prediction_settings, monkeypatch):
        # Not in both places: on a failure the same output is the error message.
        monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_MODE', 'exit')
        result = predict(prediction_settings)
        assert not result.succeeded
        assert all(run.warnings == [] for run in result.runs)


class TestAllelesPerProcess:
    """Scoring several alleles in one netMHCpan process, which is most of what a long run costs."""

    PANEL = ['HLA-A*02:01', 'HLA-B*07:02', 'HLA-C*07:02']

    def record_batches(self, monkeypatch):
        """list: the alleles each netMHCpan process was given, as the run goes."""
        predict_module = importlib.import_module('netmhcpan_multiprocessor.predict')
        real_run_chunk = predict_module.run_chunk
        batches = []

        def recording(*arguments, alleles, **keywords):
            batches.append(list(alleles))
            return real_run_chunk(*arguments, alleles=alleles, **keywords)

        monkeypatch.setattr(predict_module, 'run_chunk', recording)
        return batches

    def test_one_allele_per_process_by_default(self, prediction_settings, monkeypatch):
        batches = self.record_batches(monkeypatch)
        predict({**prediction_settings, 'alleles': self.PANEL})
        assert all(len(batch) == 1 for batch in batches)

    def test_a_batch_scores_its_alleles_in_one_process(self, prediction_settings, monkeypatch):
        batches = self.record_batches(monkeypatch)
        predict({**prediction_settings, 'alleles': self.PANEL, 'alleles_per_process': 3})
        assert all(len(batch) == 3 for batch in batches)
        # Three alleles over two lengths, in chunks of ten: six processes rather
        # than the eighteen the same work takes one allele at a time.
        assert len(batches) == 6

    def test_zero_means_every_allele_in_one_batch(self, prediction_settings, monkeypatch):
        batches = self.record_batches(monkeypatch)
        predict({**prediction_settings, 'alleles': self.PANEL, 'alleles_per_process': 0})
        assert all(len(batch) == len(self.PANEL) for batch in batches)

    def test_a_batch_larger_than_the_panel_is_the_panel(self, prediction_settings, monkeypatch):
        batches = self.record_batches(monkeypatch)
        predict({**prediction_settings, 'alleles': self.PANEL, 'alleles_per_process': 99})
        assert all(len(batch) == len(self.PANEL) for batch in batches)

    def test_batching_changes_nothing_about_the_output(self, prediction_settings, tmp_path):
        one_at_a_time = tmp_path / 'separate'
        batched = tmp_path / 'batched'
        predict({**prediction_settings, 'alleles': self.PANEL, 'output_folder': one_at_a_time})
        predict({
            **prediction_settings, 'alleles': self.PANEL,
            'output_folder': batched, 'alleles_per_process': 0,
        })
        written = sorted(path.relative_to(batched) for path in batched.rglob('*.parquet'))
        assert written == sorted(path.relative_to(one_at_a_time)
                                 for path in one_at_a_time.rglob('*.parquet'))
        for path in written:
            assert pl.read_parquet(batched / path).equals(pl.read_parquet(one_at_a_time / path))

    def test_every_allele_is_reported_as_its_own_run(self, prediction_settings):
        result = predict({
            **prediction_settings, 'alleles': self.PANEL, 'alleles_per_process': 0,
        })
        # Three alleles over two lengths, whatever they shared a process with.
        assert len(result.runs) == 6
        assert {run.allele for run in result.runs} == set(self.PANEL)
        assert all(run.peptide_count == 30 for run in result.runs)

    def test_a_resumed_run_batches_only_what_is_missing(self, prediction_settings,
                                                        output_folder, monkeypatch):
        # The point of batching after the skip: a run that stopped part way
        # fills its processes with the alleles still to do.
        predict({**prediction_settings, 'alleles': self.PANEL[:1]})
        batches = self.record_batches(monkeypatch)
        result = predict({
            **prediction_settings, 'alleles': self.PANEL,
            'alleles_per_process': 0, 'overwrite': False,
        })
        assert all(batch == ['HLA-B07:02', 'HLA-C07:02'] for batch in batches)
        assert {run.allele for run in result.runs if run.skipped} == {'HLA-A*02:01'}

    def test_a_failed_batch_fails_every_allele_in_it(self, prediction_settings,
                                                     output_folder, monkeypatch):
        # The trade-off of sharing a process: one allele netMHCpan cannot do
        # takes its batch down with it. Alleles are checked against netMHCpan's
        # own list before any work starts, so this is the rare case.
        monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_MODE', 'exit')
        monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_ALLELE', 'HLA-C07:02')
        result = predict({
            **prediction_settings, 'alleles': self.PANEL, 'alleles_per_process': 0,
        })
        assert {run.allele for run in result.failures} == set(self.PANEL)
        assert not list(output_folder.rglob('*.parquet'))

    def test_one_allele_per_process_still_isolates_a_failure(self, prediction_settings,
                                                             output_folder, monkeypatch):
        monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_MODE', 'exit')
        monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_ALLELE', 'HLA-C07:02')
        result = predict({**prediction_settings, 'alleles': self.PANEL})
        assert {run.allele for run in result.failures} == {'HLA-C*07:02'}
        assert (output_folder / 'hla_a_02_01' / '9mers' / 'testset.parquet').exists()

    def test_a_batch_reports_progress_for_every_allele(self, prediction_settings):
        class Reporter:
            def __init__(self):
                self.totals = []
                self.advanced = 0

            def add_run(self, label, total):
                self.totals.append(total)
                return label

            def advance(self, handle, amount=1):
                self.advanced += amount

            def fail(self, handle, message):
                raise AssertionError(f"unexpected failure: {message}")

        reporter = Reporter()
        predict(
            {**prediction_settings, 'alleles': self.PANEL, 'alleles_per_process': 0},
            progress=reporter,
        )
        # One bar per length, each covering its three alleles' 30 peptides.
        assert reporter.totals == [90, 90]
        assert reporter.advanced == 180

    def test_a_negative_batch_is_refused(self, prediction_settings):
        with pytest.raises(ConfigError, match='alleles_per_process cannot be negative'):
            predict({**prediction_settings, 'alleles_per_process': -1})
