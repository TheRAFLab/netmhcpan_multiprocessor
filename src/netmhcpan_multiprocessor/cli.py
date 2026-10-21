"""
The command line interface.

A thin layer over the same functions a pipeline calls: every command builds a
config and hands it straight to predict() or prepare(). The rich progress bar
and summary table live here rather than in the library, so importing the
package into a pipeline brings no terminal output with it.
"""

import click
from click.core import ParameterSource
from rich.console import Console
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
)
from rich.table import Table

from netmhcpan_multiprocessor.alleles import resolve_alleles
from netmhcpan_multiprocessor.config import (
    NETMHCPAN_ENVIRONMENT_VARIABLE,
    PredictionConfig,
    PrepareConfig,
    normalise_config_keys,
    read_config_file,
)
from netmhcpan_multiprocessor.errors import NetMHCpanMultiprocessorError
from netmhcpan_multiprocessor.modes import DEFAULT_MODE, MODES, get_preparer, registered_modes
from netmhcpan_multiprocessor.parsing import DEFAULT_SCORES, SCORE_SETS
from netmhcpan_multiprocessor.predict import predict as run_prediction


def _mode_source_help() -> str:
    # Built from the registry rather than written out, so --help cannot drift
    # from what is actually registered.
    return '; '.join(f"{name}: {', '.join(mode.sources)}" for name, mode in MODES.items())


def _score_set_help() -> str:
    return '; '.join(f"{name}: {each.description}" for name, each in SCORE_SETS.items())


class RichProgressReporter:
    """Shows one progress bar per run, advancing as predictions stream back."""

    def __init__(self, console: Console):
        self.progress = Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            MofNCompleteColumn(),
            TextColumn("peptides"),
            TimeElapsedColumn(),
            console=console,
        )

    def __enter__(self):
        self.progress.__enter__()
        return self

    def __exit__(self, *exception):
        return self.progress.__exit__(*exception)

    def add_run(self, label: str, total: int):
        """
        Add a progress bar for a run.

        Args:
            label (str): the run's name.
            total (int): how many peptides it will predict.

        Returns:
            TaskID: the rich task to advance.
        """
        return self.progress.add_task(label, total=total)

    def advance(self, handle, amount: int = 1) -> None:
        """
        Advance a run's bar.

        Args:
            handle (TaskID): the run's task.
            amount (int): how many peptides have been predicted since the last call.

        Returns:
            None
        """
        self.progress.advance(handle, amount)

    def fail(self, handle, message: str) -> None:
        """
        Stop a failed run's bar, so it does not sit spinning at the end.

        Args:
            handle (TaskID): the run's task.
            message (str): why it failed; shown in the summary table instead.

        Returns:
            None
        """
        self.progress.stop_task(handle)


def _explicit_options(context: click.Context, names: tuple[str, ...]) -> dict:
    # Only options actually typed are merged over the config file. Using click's
    # own parameter source rather than testing for None is what lets a file set
    # a value that has a non-None default, such as workers or chunk_size.
    return {
        name: context.params[name]
        for name in names
        if context.get_parameter_source(name) is not ParameterSource.DEFAULT
    }


def _build_settings(
    context: click.Context,
    names: tuple[str, ...],
    config_file,
    section: str,
    config_class,
) -> dict:
    # The file's aliases are resolved before the merge, so a file using 'allele'
    # and a command line using --allele are the same key, and the option wins
    # rather than the two looking like a conflict.
    from_file = read_config_file(config_file, section) if config_file else {}
    return {
        **normalise_config_keys(from_file, config_class),
        **_explicit_options(context, names),
    }


@click.group()
@click.version_option(package_name='netmhcpan-multiprocessor')
def cli() -> None:
    """Run netMHCpan over peptide parquet sets, in parallel."""


@cli.command()
@click.option('--input-folder', 'input_folder', type=click.Path(file_okay=False),
              help="Folder holding the peptide set, laid out as the source expects.")
@click.option('--output-folder', 'output_folder', type=click.Path(file_okay=False),
              help="Written as <output>/<allele_slug>/<N>mers/<dataset>.parquet.")
@click.option('--allele', 'alleles', multiple=True,
              help="Allele to predict, as a number ('HLA-A*02:01') or a slug ('hla_a_02_01'). Repeat for several.")
@click.option('--length', 'lengths', multiple=True, type=int,
              help="Peptide length to predict. Repeat for several. Default: every length found.")
@click.option('--mode', 'mode', default=DEFAULT_MODE, show_default=True,
              type=click.Choice(registered_modes()),
              help="The kind of peptide collection.")
@click.option('--source', 'source', default=None,
              help=f"Input handler for the mode. Default: the mode's own. [{_mode_source_help()}]")
@click.option('--scores', 'scores', default=DEFAULT_SCORES, show_default=True,
              type=click.Choice(tuple(SCORE_SETS)),
              help=f"Which predictions to make. [{_score_set_help()}]")
# No envvar= here: click would report an environment value as an explicitly given
# option, which would then beat a netmhcpan_path set in --config. The environment
# is a fallback, and resolve_netmhcpan_path already reads it at that precedence.
@click.option('--netmhcpan', 'netmhcpan_path', type=click.Path(dir_okay=False),
              help=f"The netMHCpan executable. Falls back to "
                   f"${NETMHCPAN_ENVIRONMENT_VARIABLE}, then PATH.")
@click.option('--chunk-size', 'chunk_size', type=int, default=5000, show_default=True,
              help="Peptides per netMHCpan process. Smaller chunks balance the workers better.")
@click.option('--alleles-per-process', 'alleles_per_process', type=int,
              default=PredictionConfig.alleles_per_process, show_default=True,
              help="Alleles scored in one netMHCpan process, or 0 for all of them. A batch "
                   "costs far less than a process each, but multiplies what one chunk costs, "
                   "so raise it alongside a smaller --chunk-size.")
@click.option('--workers', 'workers', type=int, default=PredictionConfig.workers, show_default=True,
              help="Concurrent netMHCpan processes.")
@click.option('--peptide-column', 'peptide_column', default='peptide_sequence', show_default=True,
              help="Column holding the peptide sequences.")
@click.option('--raw-output-folder', 'raw_output_folder', type=click.Path(file_okay=False),
              help="Also write netMHCpan's own text output here. Off by default.")
@click.option('--validate-alleles/--no-validate-alleles', 'validate_alleles', default=True,
              show_default=True, help="Check the alleles against netMHCpan's allele list first.")
@click.option('--overwrite/--no-overwrite', 'overwrite', default=True, show_default=True,
              help="--no-overwrite skips any run whose output parquet already exists.")
@click.option('--strict', 'strict', is_flag=True, default=False,
              help="Stop and raise on the first failed run, rather than finishing the rest.")
@click.option('--config', 'config_file', type=click.Path(dir_okay=False),
              help="TOML or JSON file of the same settings. Options given here win.")
@click.option('--quiet', is_flag=True, default=False, help="No progress bar and no summary table.")
@click.option('--verbose', '-v', is_flag=True, default=False,
              help="Show the resolved settings, alleles and input files before running, "
                   "and anything netMHCpan wrote to stderr afterwards.")
@click.pass_context
def predict(context: click.Context, config_file, quiet, verbose, **_options) -> None:
    """
    Predict a peptide set against one or more alleles.

    Reads one peptide parquet per length and writes one prediction parquet per
    length per allele, in the same row order, so the two join directly.
    """
    _check_output_flags(quiet, verbose)
    console = Console(quiet=quiet)

    try:
        settings = _build_settings(
            context, PREDICT_OPTIONS, config_file, 'predict', PredictionConfig)
        config = PredictionConfig.from_dict(settings)
    except NetMHCpanMultiprocessorError as configuration_error:
        raise click.ClickException(str(configuration_error)) from configuration_error

    console.print(
        f"Predicting [bold]{config.input_folder}[/bold] "
        f"({config.mode}/{config.source_name}, {config.scores.upper()} scores) for "
        f"[bold]{len(config.alleles)}[/bold] allele(s) over "
        f"[bold]{config.workers}[/bold] workers"
    )
    if verbose:
        # Before the run rather than after: the point of asking is to check what
        # is about to happen, while stopping it is still cheap.
        try:
            _describe_prediction(console, config)
        except NetMHCpanMultiprocessorError as inspection_error:
            raise click.ClickException(str(inspection_error)) from inspection_error

    try:
        if quiet:
            result = run_prediction(config)
        else:
            with RichProgressReporter(console) as reporter:
                result = run_prediction(config, progress=reporter)
    except NetMHCpanMultiprocessorError as prediction_error:
        raise click.ClickException(str(prediction_error)) from prediction_error

    if not quiet:
        console.print(_prediction_summary(result))
    if verbose:
        _report_warnings(console, result)
    if result.failures:
        raise click.ClickException(
            f"{len(result.failures)} of {len(result.runs)} predictions failed"
        )


PREDICT_OPTIONS = (
    'input_folder', 'output_folder', 'alleles', 'lengths', 'mode', 'source',
    'scores', 'netmhcpan_path', 'chunk_size', 'alleles_per_process', 'workers',
    'peptide_column', 'raw_output_folder', 'validate_alleles', 'overwrite', 'strict',
)


@cli.command()
@click.option('--input-path', 'input_path', type=click.Path(dir_okay=False),
              help="The file to build a peptide set from; a FASTA, for the proteome mode.")
@click.option('--output-folder', 'output_folder', type=click.Path(file_okay=False),
              help="Written as <output>/<N>mers/<dataset>.parquet, ready to predict from.")
@click.option('--length', 'lengths', multiple=True, type=int,
              help="Peptide length to generate. Repeat for several. [default: 8 9 10 11]")
@click.option('--mode', 'mode', default=DEFAULT_MODE, show_default=True,
              type=click.Choice(registered_modes()),
              help="The kind of peptide collection to prepare.")
@click.option('--dataset', 'dataset', default=None,
              help="Name for the peptide set. Default: the input file's stem.")
@click.option('--keep-nonstandard', 'keep_nonstandard', is_flag=True, default=False,
              help="Keep peptides with residues outside the twenty standard amino acids.")
@click.option('--config', 'config_file', type=click.Path(dir_okay=False),
              help="TOML or JSON file of the same settings. Options given here win.")
@click.option('--quiet', is_flag=True, default=False, help="No summary table.")
@click.option('--verbose', '-v', is_flag=True, default=False,
              help="Show the resolved settings before building.")
@click.pass_context
def prepare(context: click.Context, config_file, quiet, verbose, **_options) -> None:
    """
    Build a peptide set for prediction.

    For the proteome mode: every k-mer of every sequence in a FASTA, written as
    one parquet per length.
    """
    _check_output_flags(quiet, verbose)
    console = Console(quiet=quiet)

    try:
        settings = _build_settings(
            context, PREPARE_OPTIONS, config_file, 'prepare', PrepareConfig)
        config = PrepareConfig.from_dict(settings)
        console.print(
            f"Preparing [bold]{config.dataset}[/bold] ({config.mode}) from {config.input_path}"
        )
        if verbose:
            console.print(_settings_table(config, 'Settings'))
        result = get_preparer(config.mode).prepare(config)
    except NetMHCpanMultiprocessorError as preparation_error:
        raise click.ClickException(str(preparation_error)) from preparation_error

    if not quiet:
        console.print(_prepare_summary(result))


PREPARE_OPTIONS = ('input_path', 'output_folder', 'lengths', 'mode', 'dataset', 'keep_nonstandard')


def _check_output_flags(quiet: bool, verbose: bool) -> None:
    if quiet and verbose:
        raise click.UsageError("--quiet and --verbose ask for opposite things; use one")


def _settings_table(config, title: str) -> Table:
    # Every field of the config, so what is shown is what the run will use -
    # including the values that came from a file, an alias, or a default, which
    # are exactly the ones worth checking.
    settings = Table(title=title, show_header=False, title_justify='left')
    settings.add_column('setting', style='bold')
    settings.add_column('value')
    for name in config.__dataclass_fields__:
        settings.add_row(name, str(getattr(config, name)))
    return settings


def _describe_prediction(console: Console, config) -> None:
    console.print(_settings_table(config, 'Settings'))

    alleles = Table(title='Alleles', title_justify='left')
    alleles.add_column('given')
    alleles.add_column('netMHCpan -a')
    alleles.add_column('slug, and output folder')
    for supplied, allele in zip(
        config.alleles,
        resolve_alleles(config.alleles, config.netmhcpan_path, config.validate_alleles),
    ):
        alleles.add_row(supplied, allele.netmhcpan_name, allele.slug)
    console.print(alleles)

    found = Table(title='Peptide files', title_justify='left')
    found.add_column('length')
    found.add_column('dataset')
    found.add_column('path')
    for peptide_file in config.peptide_source.discover(config.input_folder, config.lengths):
        found.add_row(f"{peptide_file.length}mers", peptide_file.dataset, str(peptide_file.path))
    console.print(found)


def _report_warnings(console: Console, result) -> None:
    reported = [run for run in result.runs if run.warnings]
    if not reported:
        console.print("netMHCpan wrote nothing to stderr.")
        return
    for run in reported:
        console.print(f"[yellow]netMHCpan said, during {run.label}:[/yellow]")
        for warning in run.warnings:
            console.print(f"  {warning}")


def _prediction_summary(result) -> Table:
    summary = Table(title=f"netMHCpan finished in {result.elapsed_seconds:.1f}s")
    summary.add_column("run")
    summary.add_column("chunks", justify="right")
    summary.add_column("predictions", justify="right")
    summary.add_column("SB", justify="right")
    summary.add_column("WB", justify="right")
    summary.add_column("parquet")
    for run in result.runs:
        if run.error:
            summary.add_row(run.label, str(run.chunk_count), "[red]failed[/red]", "", "", run.error)
        elif run.skipped:
            summary.add_row(run.label, "", "[yellow]skipped[/yellow]", "", "", str(run.parquet_path))
        else:
            summary.add_row(
                run.label,
                str(run.chunk_count),
                f"{run.peptide_count:,}",
                f"{run.strong_binders:,}",
                f"{run.weak_binders:,}",
                str(run.parquet_path),
            )
    return summary


def _prepare_summary(result) -> Table:
    summary = Table(title=f"Prepared {result.dataset} in {result.elapsed_seconds:.1f}s")
    summary.add_column("length")
    summary.add_column("peptides", justify="right")
    summary.add_column("parquet")
    for prepared in result.files:
        summary.add_row(f"{prepared.length}mers", f"{prepared.peptide_count:,}", str(prepared.path))
    return summary


if __name__ == '__main__':
    cli()
