"""
The prediction run itself.

Peptides in as parquet, predictions out as parquet. Everything between is about
using more than one core well: each netMHCpan process is single threaded and the
lengths are very uneven - 11mers take nearly ten times as long as 9mers - so the
work is cut into chunks, pooled flat, and handed out longest first.
"""

import os
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from netmhcpan_multiprocessor.alleles import ResolvedAllele, resolve_alleles
from netmhcpan_multiprocessor.chunking import PeptideChunk, write_chunk_files
from netmhcpan_multiprocessor.config import PredictionConfig, normalise_config_keys
from netmhcpan_multiprocessor.errors import (
    NetMHCpanMultiprocessorError,
    PredictionFailedError,
)
from netmhcpan_multiprocessor.modes.base import PeptideFile
from netmhcpan_multiprocessor.parsing import build_predictions_frame, write_raw_predictions
from netmhcpan_multiprocessor.progress import (
    GuardedProgressReporter,
    NullProgressReporter,
    ProgressReporter,
)
from netmhcpan_multiprocessor.runner import run_chunk


# How many distinct netMHCpan warnings to keep per run. Enough to see what it is
# complaining about, few enough that a program stuck in a warning loop cannot
# fill memory with it.
MAX_WARNINGS_PER_RUN = 20


@dataclass
class RunResult:
    """
    What one allele-by-length prediction produced.

    Attributes:
        length (int): the peptide length predicted.
        allele (str): the allele number, e.g. 'HLA-A*02:01'.
        allele_slug (str): the slug naming its folder and its column value.
        dataset (str): the peptide set's name, from the input file name.
        parquet_path (Path): where the predictions were written, or would have been.
        raw_path (Path | None): the netMHCpan text output, if it was kept.
        peptide_count (int): how many peptides were predicted.
        chunk_count (int): how many netMHCpan processes it took.
        strong_binders (int): predictions marked SB.
        weak_binders (int): predictions marked WB.
        error (str | None): why the run failed, or None if it did not.
        skipped (bool): True if the output already existed and overwrite was off.
        warnings (list[str]): distinct lines netMHCpan wrote to stderr while
            succeeding. Usually empty, and never fatal, but the first place to
            look when a run finishes and the numbers seem wrong.
    """
    length: int
    allele: str
    allele_slug: str
    dataset: str
    parquet_path: Path
    raw_path: Path | None = None
    peptide_count: int = 0
    chunk_count: int = 0
    strong_binders: int = 0
    weak_binders: int = 0
    error: str | None = None
    skipped: bool = False
    warnings: list[str] = field(default_factory=list)

    @property
    def succeeded(self) -> bool:
        """bool: True if the run produced predictions, or was deliberately skipped."""
        return self.error is None

    @property
    def label(self) -> str:
        """str: a short human-readable name for the run."""
        return f"{self.dataset} {self.length}mers {self.allele}"


@dataclass
class PredictionResult:
    """
    What a whole prediction run produced.

    Attributes:
        runs (list[RunResult]): one entry per allele-by-length prediction, in
            the order they were planned.
        elapsed_seconds (float): how long the run took.
    """
    runs: list[RunResult] = field(default_factory=list)
    elapsed_seconds: float = 0.0

    @property
    def failures(self) -> list[RunResult]:
        """list[RunResult]: the runs that failed."""
        return [run for run in self.runs if run.error is not None]

    @property
    def succeeded(self) -> bool:
        """bool: True if no run failed."""
        return not self.failures

    @property
    def peptide_count(self) -> int:
        """int: predictions written across every run."""
        return sum(run.peptide_count for run in self.runs)

    def frame(self) -> pl.DataFrame:
        """
        The summary as data, for a pipeline to log or store.

        Returns:
            pl.DataFrame: one row per run, with the paths as strings so the
                frame can be written straight out.
        """
        return pl.DataFrame([
            {
                'dataset': run.dataset,
                'length': run.length,
                'allele': run.allele,
                'allele_slug': run.allele_slug,
                'parquet_path': str(run.parquet_path),
                'raw_path': str(run.raw_path) if run.raw_path else None,
                'peptide_count': run.peptide_count,
                'chunk_count': run.chunk_count,
                'strong_binders': run.strong_binders,
                'weak_binders': run.weak_binders,
                'skipped': run.skipped,
                'error': run.error,
            }
            for run in self.runs
        ], schema={
            'dataset': pl.String, 'length': pl.Int64, 'allele': pl.String,
            'allele_slug': pl.String, 'parquet_path': pl.String, 'raw_path': pl.String,
            'peptide_count': pl.Int64, 'chunk_count': pl.Int64,
            'strong_binders': pl.Int64, 'weak_binders': pl.Int64,
            'skipped': pl.Boolean, 'error': pl.String,
        })


@dataclass(frozen=True)
class _PlannedRun:
    """One allele-by-length prediction, and everything needed to carry it out."""
    peptide_file: PeptideFile
    allele: ResolvedAllele
    parquet_path: Path
    raw_path: Path | None

    @property
    def key(self) -> tuple[int, str, str]:
        return (self.peptide_file.length, self.peptide_file.dataset, self.allele.slug)

    @property
    def label(self) -> str:
        return f"{self.peptide_file.length}mers {self.allele.number}"


def predict(
    config: PredictionConfig | dict | None = None,
    *,
    progress: ProgressReporter | None = None,
    **settings,
) -> PredictionResult:
    """
    Run netMHCpan over a peptide set for one or more alleles.

    Args:
        config (PredictionConfig | dict | None): the run's settings, as a config
            object or the dictionary a pipeline is likely to already hold.
        progress (ProgressReporter | None): where to report progress, or None
            for silence.
        **settings: config keys given as keyword arguments, merged over config.

    Returns:
        PredictionResult: one RunResult per allele-by-length prediction, whether
            it succeeded or failed. Nothing is raised for a failed run unless
            the config asks for it, so one bad allele does not lose the rest of
            a long run's work.

    Raises:
        ConfigError: if the settings are unusable.
        UnknownAlleleError: if an allele cannot be resolved.
        PeptideInputError: if the peptide input is missing or malformed.
        PredictionFailedError: if a run failed and strict is set.
    """
    config = _as_config(config, settings)
    # Guarded, because the reporter comes from the caller and a display that
    # raises must not be able to lose the results of a long run.
    reporter = GuardedProgressReporter(progress or NullProgressReporter())

    alleles = resolve_alleles(config.alleles, config.netmhcpan_path, config.validate_alleles)
    source = config.peptide_source
    peptide_files = source.discover(config.input_folder, config.lengths)

    planned, results = _plan_runs(config, peptide_files, alleles)
    if not planned:
        return PredictionResult(runs=_ordered_results(results, peptide_files, alleles))

    started = time.perf_counter()
    # Chunks live only for the run: the peptide parquet is the input of record,
    # and a half-written chunk folder left behind would be a trap for the next one.
    with tempfile.TemporaryDirectory(prefix='netmhcpan_chunks_') as chunk_folder:
        chunks_by_file = _chunk_peptide_files(config, source, planned, Path(chunk_folder))
        results.update(_run_predictions(config, planned, chunks_by_file, reporter))
    elapsed = time.perf_counter() - started

    result = PredictionResult(
        runs=_ordered_results(results, peptide_files, alleles),
        elapsed_seconds=elapsed,
    )
    if config.strict and result.failures:
        raise PredictionFailedError(
            f"{len(result.failures)} of {len(result.runs)} predictions failed: "
            + '; '.join(f"{run.label}: {run.error}" for run in result.failures)
        )
    return result


def _as_config(config, settings: dict) -> PredictionConfig:
    if isinstance(config, PredictionConfig) and not settings:
        return config
    # Each mapping's aliases are resolved before the merge, so 'allele' given as
    # a keyword overrides an 'alleles' already on the config rather than looking
    # like two keys fighting over one field.
    settings = normalise_config_keys(settings, PredictionConfig)
    if isinstance(config, PredictionConfig):
        # Keyword arguments over a config object: rebuilt through from_dict so
        # the merged values get the same validation as any other route in.
        base = {name: getattr(config, name) for name in config.__dataclass_fields__}
        return PredictionConfig.from_dict({**base, **settings})
    return PredictionConfig.from_dict({**normalise_config_keys(config or {}, PredictionConfig), **settings})


def _plan_runs(
    config: PredictionConfig,
    peptide_files: list[PeptideFile],
    alleles: list[ResolvedAllele],
) -> tuple[list[_PlannedRun], dict[tuple, RunResult]]:
    # Returns what still has to be run, and the results of everything that does
    # not - today only the runs skipped because their output is already there.
    planned: list[_PlannedRun] = []
    skipped: dict[tuple, RunResult] = {}
    for peptide_file in peptide_files:
        for allele in alleles:
            parquet_path = (
                config.output_folder / allele.slug
                / f"{peptide_file.length}mers" / f"{peptide_file.dataset}.parquet"
            )
            raw_path = None
            if config.raw_output_folder is not None:
                raw_path = config.raw_output_folder / (
                    f"{peptide_file.dataset}__{peptide_file.length}mers_{allele.slug}.txt"
                )
            run = _PlannedRun(peptide_file, allele, parquet_path, raw_path)
            if not config.overwrite and parquet_path.exists():
                skipped[run.key] = RunResult(
                    length=peptide_file.length,
                    allele=allele.number,
                    allele_slug=allele.slug,
                    dataset=peptide_file.dataset,
                    parquet_path=parquet_path,
                    skipped=True,
                )
                continue
            planned.append(run)
    return planned, skipped


def _chunk_peptide_files(
    config: PredictionConfig,
    source,
    planned: list[_PlannedRun],
    chunk_folder: Path,
) -> dict[Path, list[PeptideChunk]]:
    # Chunked once per peptide file rather than once per allele: the chunks are
    # the same peptides whichever allele is being predicted, so every extra
    # allele would otherwise rewrite the whole peptide set to disk again.
    chunks_by_file: dict[Path, list[PeptideChunk]] = {}
    for run in planned:
        peptide_file = run.peptide_file
        if peptide_file.path in chunks_by_file:
            continue
        peptides = source.read_peptides(peptide_file, config.peptide_column)
        chunks_by_file[peptide_file.path] = write_chunk_files(
            peptides,
            chunk_size=config.chunk_size,
            chunk_folder=chunk_folder,
            name=f"{peptide_file.dataset}__{peptide_file.length}mers",
        )
    return chunks_by_file


def _run_predictions(
    config: PredictionConfig,
    planned: list[_PlannedRun],
    chunks_by_file: dict[Path, list[PeptideChunk]],
    reporter: ProgressReporter,
) -> dict[tuple, RunResult]:
    handles = {}
    for run in planned:
        peptide_count = sum(chunk.peptide_count for chunk in chunks_by_file[run.peptide_file.path])
        handles[run.key] = reporter.add_run(run.label, peptide_count)

    results: dict[tuple, RunResult] = {}
    rows_by_run: dict[tuple, dict[int, list[str]]] = defaultdict(dict)
    # Appended to from the worker threads. list.append is atomic under the GIL,
    # and the same warning from every chunk is deduplicated when the run is
    # assembled, so a bounded number of distinct lines survives however many
    # chunks repeat them.
    warnings_by_run: dict[tuple, list[str]] = defaultdict(list)
    outstanding = {run.key: len(chunks_by_file[run.peptide_file.path]) for run in planned}
    failures: dict[tuple, str] = {}
    runs_by_key = {run.key: run for run in planned}

    # Longest job first: submit the slowest lengths before the quick ones, so the
    # long tail starts early and the short chunks fill in around it. Submitting
    # 8mers first would leave the 11mer chunks queued last, which is what decides
    # the finishing time.
    longest_first = sorted(planned, key=lambda run: run.peptide_file.length, reverse=True)

    # Threads rather than processes because the work happens in a subprocess, so
    # the GIL is released for the whole of the wait. One flat pool of chunk jobs
    # rather than a pool per length, so a worker that finishes a quick 9mer chunk
    # picks up an 11mer chunk straight away instead of waiting on its own length.
    executor = ThreadPoolExecutor(max_workers=config.workers)
    try:
        futures = {}
        for run in longest_first:
            for chunk_index, chunk in enumerate(chunks_by_file[run.peptide_file.path]):
                future = executor.submit(
                    run_chunk,
                    netmhcpan_path=config.netmhcpan_path,
                    chunk=chunk,
                    allele=run.allele.netmhcpan_name,
                    advance=lambda amount, key=run.key: reporter.advance(handles[key], amount),
                    arguments=config.score_set.arguments,
                    on_warning=lambda line, key=run.key: warnings_by_run[key].append(line),
                )
                futures[future] = (run.key, chunk_index)

        for future in as_completed(futures):
            run_key, chunk_index = futures[future]
            # Recording failures per run means one bad allele cannot pass for a
            # successful run, while the other runs still finish. Every exception
            # is caught, not just this package's own: a chunk that fails for a
            # reason nobody anticipated - a reporter that raises, undecodable
            # output - must still leave the other runs their results, which is
            # what a caller is promised when strict is off.
            try:
                rows_by_run[run_key][chunk_index] = future.result()
            except Exception as chunk_error:
                if run_key not in failures:
                    failures[run_key] = _describe(chunk_error)
                    reporter.fail(handles[run_key], failures[run_key])

            outstanding[run_key] -= 1
            if outstanding[run_key]:
                continue
            # Assembled and written as soon as its last chunk lands, rather than
            # after the whole pool drains, so a finished length's rows can be
            # released instead of every run's rows piling up until the end.
            results[run_key] = _finalise_run(
                runs_by_key[run_key],
                rows_by_run.pop(run_key, {}),
                chunk_count=len(chunks_by_file[runs_by_key[run_key].peptide_file.path]),
                score_set=config.score_set,
                error=failures.get(run_key),
                warnings=warnings_by_run.pop(run_key, []),
            )
            if config.strict and results[run_key].error:
                # Nothing queued is worth running once the result is known to be
                # a failure, so give back the workers rather than the wall clock.
                executor.shutdown(wait=False, cancel_futures=True)
                break
    except BaseException:
        # Ctrl-C, most likely. Without cancelling the queue, shutting down would
        # quietly work through every chunk still waiting - minutes or hours of
        # what looks like a hang - and then discard the results anyway.
        executor.shutdown(wait=False, cancel_futures=True)
        raise
    finally:
        executor.shutdown(wait=True)

    # Anything still outstanding was abandoned by a strict stop. Some of it may
    # have finished, since a run is only assembled once its last chunk is
    # collected, so this says the run was abandoned rather than that it was
    # never started.
    for run_key in outstanding:
        if run_key not in results:
            run = runs_by_key[run_key]
            error = failures.get(run_key, 'abandoned: an earlier failure stopped the run')
            if run_key not in failures:
                reporter.fail(handles[run_key], error)
            # Abandoned counts as failed, so its path must not keep an earlier
            # run's parquet either.
            _clear_stale_output(run)
            results[run_key] = _run_result(
                run,
                peptide_count=0,
                chunk_count=len(chunks_by_file[run.peptide_file.path]),
                error=error,
            )
    return results


def _describe(error: Exception) -> str:
    # An unanticipated failure needs its type in the message, since something
    # like a bare KeyError says nothing on its own.
    if isinstance(error, NetMHCpanMultiprocessorError):
        return str(error)
    return f"{type(error).__name__}: {error}"


def _finalise_run(
    run: _PlannedRun,
    rows_by_chunk: dict[int, list[str]],
    chunk_count: int,
    score_set,
    error: str | None,
    warnings: list[str] | None = None,
) -> RunResult:
    # Deduplicated and capped: every chunk of a long run tends to repeat the same
    # complaint, and the point is to see what netMHCpan said, not how often.
    warnings = list(dict.fromkeys(warnings or []))[:MAX_WARNINGS_PER_RUN]
    if error is not None:
        _clear_stale_output(run)
        return _run_result(
            run, peptide_count=0, chunk_count=chunk_count, error=error, warnings=warnings,
        )

    # Chunks finish out of order, so reassemble by chunk index to put the
    # predictions back into peptide input order. Nothing sorts them after this,
    # so the frame and the parquet keep the input peptide file's order, and join
    # back to the peptide parquet row for row.
    merged_rows = [row for chunk_index in sorted(rows_by_chunk) for row in rows_by_chunk[chunk_index]]

    try:
        if run.raw_path is not None:
            run.raw_path.parent.mkdir(parents=True, exist_ok=True)
            write_raw_predictions(merged_rows, run.raw_path)

        predictions = build_predictions_frame(merged_rows, run.allele.slug, score_set)
        _write_parquet(predictions, run.parquet_path)
    except Exception as writing_error:
        # Broad, because parsing goes through polars, which raises its own
        # exceptions for a value that will not cast. Letting one of those out
        # would lose every other run's results on the way past.
        _clear_stale_output(run)
        return _run_result(
            run, peptide_count=0, chunk_count=chunk_count,
            error=_describe(writing_error), warnings=warnings,
        )

    binder_counts = predictions['bind_level'].value_counts().to_dict(as_series=False)
    binder_counts = dict(zip(binder_counts['bind_level'], binder_counts['count']))
    return _run_result(
        run,
        peptide_count=predictions.height,
        chunk_count=chunk_count,
        error=None,
        strong_binders=binder_counts.get('SB', 0),
        weak_binders=binder_counts.get('WB', 0),
        wrote_raw=run.raw_path is not None,
        warnings=warnings,
    )


def _write_parquet(predictions: pl.DataFrame, parquet_path: Path) -> None:
    # Written beside the target and moved into place, because os.replace is
    # atomic within a filesystem: an interrupted run leaves either the previous
    # file or none, never a half written one that reads as a short peptide set.
    parquet_path.parent.mkdir(parents=True, exist_ok=True)
    # Named for this process, so that two runs writing the same output at once
    # cannot publish or delete each other's half-written file.
    partial_path = parquet_path.with_name(f"{parquet_path.name}.{os.getpid()}.partial")
    try:
        predictions.write_parquet(partial_path)
        os.replace(partial_path, parquet_path)
    finally:
        partial_path.unlink(missing_ok=True)


def _clear_stale_output(run: _PlannedRun) -> None:
    # A failed run must not leave the previous run's parquet sitting at the path
    # it was asked to replace. It no longer corresponds to anything: the peptide
    # set may have changed underneath it, and the documented positional join
    # would then be silently wrong. Worse, a retry with overwrite off would skip
    # the path as though it held this run's results.
    # Best effort: a read-only output tree, or a directory sitting where the
    # parquet should go, must not turn one run's failure into the loss of every
    # other run's results on the way out.
    for path in (run.parquet_path, run.raw_path):
        if path is None:
            continue
        try:
            path.unlink(missing_ok=True)
        except OSError:
            pass


def _run_result(
    run: _PlannedRun,
    peptide_count: int,
    chunk_count: int,
    error: str | None,
    strong_binders: int = 0,
    weak_binders: int = 0,
    wrote_raw: bool = False,
    warnings: list[str] | None = None,
) -> RunResult:
    return RunResult(
        length=run.peptide_file.length,
        allele=run.allele.number,
        allele_slug=run.allele.slug,
        dataset=run.peptide_file.dataset,
        parquet_path=run.parquet_path,
        raw_path=run.raw_path if wrote_raw else None,
        peptide_count=peptide_count,
        chunk_count=chunk_count,
        strong_binders=strong_binders,
        weak_binders=weak_binders,
        error=error,
        warnings=warnings or [],
    )


def _ordered_results(
    results: dict[tuple, RunResult],
    peptide_files: list[PeptideFile],
    alleles: list[ResolvedAllele],
) -> list[RunResult]:
    # Reported in the order the runs were planned - by length, then dataset,
    # then the order the alleles were given - rather than the order they
    # happened to finish in, so two runs of the same job read the same.
    ordered = []
    for peptide_file in peptide_files:
        for allele in alleles:
            key = (peptide_file.length, peptide_file.dataset, allele.slug)
            if key in results:
                ordered.append(results[key])
    return ordered
