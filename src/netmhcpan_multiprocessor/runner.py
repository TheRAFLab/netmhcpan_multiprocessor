"""
Running netMHCpan over a single chunk of peptides.

One process per chunk, its stdout read as it arrives so predictions can be
counted while they stream. Everything about detecting a failed run lives here,
because netMHCpan has more than one way of not doing what was asked.
"""

import os
import signal
import subprocess
import threading
from collections.abc import Callable
from pathlib import Path

from netmhcpan_multiprocessor.chunking import PeptideChunk
from netmhcpan_multiprocessor.errors import NetMHCpanRunError
from netmhcpan_multiprocessor.parsing import is_prediction_row


# How long to give a netMHCpan that has closed its pipes but not yet exited.
# Generous, because this is only ever reached when something is already wrong,
# and stingy enough that one stuck process cannot hold a worker indefinitely.
SHUTDOWN_TIMEOUT_SECONDS = 30

# How long to keep collecting stderr after netMHCpan has exited. Short, because
# by then everything it wrote is already in the pipe: this only expires when
# something else it left running is still holding the write end open.
STDERR_GRACE_SECONDS = 5


def run_chunk(
    netmhcpan_path: Path,
    chunk: PeptideChunk,
    allele: str,
    advance: Callable[[int], None] | None = None,
    arguments: tuple[str, ...] = (),
    on_warning: Callable[[str], None] | None = None,
) -> list[str]:
    """
    Run netMHCpan over one chunk of peptides.

    Args:
        netmhcpan_path (Path): the netMHCpan executable.
        chunk (PeptideChunk): the chunk file and how many peptides it holds.
        allele (str): the allele in netMHCpan's own naming, e.g. 'HLA-A02:01'.
        advance (Callable[[int], None] | None): called with 1 for each prediction
            as it arrives, so a caller can show progress within the chunk.
        arguments (tuple[str, ...]): extra netMHCpan arguments, such as the -BA
            that asks for binding affinity as well as eluted ligand likelihood.
        on_warning (Callable[[str], None] | None): called with each line netMHCpan
            wrote to stderr on a run that otherwise succeeded. On a failed run the
            same output goes into the error instead, so it is never lost twice.

    Returns:
        list[str]: the chunk's prediction rows, banner and summary stripped, in
            the order the peptides were given.

    Raises:
        NetMHCpanRunError: if netMHCpan could not be started, exited non-zero,
            reported an error, or returned anything other than one prediction
            per peptide in the order they were supplied.
    """
    # No shell, and stdout read here rather than redirected with '>', so a failed
    # run raises instead of silently leaving a zero byte results file behind.
    # Popen rather than run() so predictions can be counted as they stream, which
    # turns a spinner that only proves the process is alive into a real bar.
    try:
        process = subprocess.Popen(
            [str(netmhcpan_path), '-p', str(chunk.path), '-a', allele, *arguments],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            # Undecodable bytes are replaced rather than raised on: a stray byte
            # in a warning should not lose a chunk's worth of real predictions.
            errors='replace',
            # Its own process group, so that killing it kills what it started.
            # netMHCpan is a shell wrapper around the binary, and signalling only
            # the wrapper would leave a grandchild holding the pipes open.
            start_new_session=True,
        )
    except OSError as start_error:
        raise NetMHCpanRunError(f"Could not run {netmhcpan_path}: {start_error}") from start_error

    prediction_rows: list[str] = []
    error_lines: list[str] = []
    stderr_lines: list[str] = []

    # stderr is drained by its own thread rather than read after stdout. Both are
    # pipes, so a netMHCpan that writes more than the pipe buffer holds - warnings
    # over a 5000 peptide chunk get there easily - would block on its stderr write,
    # stop producing stdout, and hang this read for good.
    try:
        stderr_reader = threading.Thread(
            target=stderr_lines.extend, args=(process.stderr,), daemon=True,
        )
        stderr_reader.start()
    except BaseException:
        _terminate(process)
        raise

    # One try, with both clauses: the except has to run before the finally, so
    # that the process is dead before anything waits on the pipes it holds.
    try:
        for line in process.stdout:
            if is_prediction_row(line):
                prediction_rows.append(line)
                if advance is not None:
                    advance(1)
            elif line.startswith('Error'):
                # netMHCpan reports an unusable allele by printing 'Error.' to
                # stdout and then exiting 0, so the return code alone never
                # sees it.
                error_lines.append(line.strip())
    except BaseException:
        # Killed rather than left running: without this the process keeps its
        # pipes open and the thread pool's shutdown waits on it for ever.
        _terminate(process)
        raise
    finally:
        process.stdout.close()
        _reap(process, stderr_reader)

    stderr = ''.join(stderr_lines).strip()
    if process.returncode != 0:
        raise NetMHCpanRunError(
            f"exit {process.returncode}: {stderr or ' '.join(error_lines) or 'no output'}"
        )
    if error_lines:
        raise NetMHCpanRunError(' '.join(error_lines))
    # Backstop for any other way netMHCpan might quietly return less than asked:
    # one prediction per peptide is the whole contract of -p mode.
    if len(prediction_rows) != chunk.peptide_count:
        raise NetMHCpanRunError(
            f"expected {chunk.peptide_count} predictions, got {len(prediction_rows)}"
        )
    _check_peptide_order(prediction_rows, chunk)

    # Only on the way out of a successful run: a chunk that failed has already
    # put this same output into its error message.
    if on_warning is not None:
        for line in stderr.splitlines():
            if line.strip():
                on_warning(line.strip())
    return prediction_rows


def _reap(process: subprocess.Popen, stderr_reader: threading.Thread) -> None:
    # The process first, since that is what the result depends on, and only then
    # the stderr drain. Waiting on the drain first would mean anything else
    # holding the pipe open - a helper the wrapper script started and left
    # running - could keep a worker busy long after netMHCpan itself had
    # finished and its output was complete.
    try:
        process.wait(timeout=SHUTDOWN_TIMEOUT_SECONDS)
    except subprocess.TimeoutExpired:
        # stdout is at EOF and it still has not exited, so nothing more is
        # coming. Waiting longer would only hold a worker hostage.
        _terminate(process)
        process.wait()

    # It has exited, so everything it wrote to stderr is already in the pipe;
    # this is a short grace to collect it, not a wait on anything.
    stderr_reader.join(timeout=STDERR_GRACE_SECONDS)
    if not stderr_reader.is_alive():
        process.stderr.close()


def _terminate(process: subprocess.Popen) -> None:
    # The whole process group, not just the process: netMHCpan is a shell wrapper
    # that starts the real binary, and killing only the wrapper would leave the
    # grandchild running with the pipes still open - which is the hang this is
    # here to prevent.
    try:
        os.killpg(os.getpgid(process.pid), signal.SIGKILL)
    except (OSError, AttributeError):
        # No process group to kill, or a platform without them.
        process.kill()


def _check_peptide_order(prediction_rows: list[str], chunk: PeptideChunk) -> None:
    # The output parquet is written in this order and joined back to the peptide
    # parquet positionally, so 'netMHCpan answers in the order it was asked' is
    # load bearing. It is also the one assumption nothing else would catch: a
    # reordered result is self-consistent, since the peptide column comes from
    # netMHCpan's own output, and would silently mismatch the input instead.
    expected = [line.strip() for line in chunk.path.read_text().splitlines() if line.strip()]
    returned = [row.split(maxsplit=3)[2] for row in prediction_rows]
    if returned != expected:
        difference = next(
            (index for index, (was, got) in enumerate(zip(expected, returned)) if was != got),
            min(len(expected), len(returned)),
        )
        raise NetMHCpanRunError(
            f"netMHCpan returned peptides in a different order to the input: "
            f"position {difference + 1} was "
            f"{expected[difference] if difference < len(expected) else 'nothing'!r}, "
            f"got {returned[difference] if difference < len(returned) else 'nothing'!r}"
        )
