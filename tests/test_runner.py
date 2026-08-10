"""
Running one chunk: the subprocess handling.

These are the failure modes that hang rather than fail, so several of them run
in a child process with a timeout. A test for a deadlock cannot be allowed to
deadlock the suite.
"""

import stat
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from netmhcpan_multiprocessor.chunking import write_chunk_files
from netmhcpan_multiprocessor.errors import NetMHCpanRunError
from netmhcpan_multiprocessor.runner import run_chunk

from conftest import FIXTURES, make_peptides


@pytest.fixture
def chunk(tmp_path):
    return write_chunk_files(make_peptides(9, 5), 5, tmp_path, 'testset')[0]


@pytest.fixture
def wrapped_netmhcpan(fake_netmhcpan):
    """Path: the stand-in behind a shell wrapper, as netMHCpan itself is shipped."""
    wrapper = FIXTURES / 'wrapped_netmhcpan.sh'
    wrapper.chmod(wrapper.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return wrapper


def run_in_a_child(source: str, environment: dict, timeout: int = 30):
    """
    Run a snippet in its own interpreter, so a hang is a failure rather than a hang.

    Args:
        source (str): the code to run.
        environment (dict): environment variables to add.
        timeout (int): seconds to allow before calling it hung.

    Returns:
        subprocess.CompletedProcess: what the child did.
    """
    import os

    return subprocess.run(
        [sys.executable, '-c', textwrap.dedent(source)],
        capture_output=True,
        text=True,
        timeout=timeout,
        env={**os.environ, **environment},
        cwd=Path(__file__).parent.parent,
    )


def test_a_chunk_comes_back_in_order(fake_netmhcpan, chunk):
    rows = run_chunk(fake_netmhcpan, chunk, 'HLA-A02:01')
    assert [row.split()[2] for row in rows] == make_peptides(9, 5)


def test_progress_is_advanced_once_per_prediction(fake_netmhcpan, chunk):
    advanced = []
    run_chunk(fake_netmhcpan, chunk, 'HLA-A02:01', advance=advanced.append)
    assert sum(advanced) == 5


def test_a_missing_executable_is_reported(tmp_path, chunk):
    with pytest.raises(NetMHCpanRunError, match='Could not run'):
        run_chunk(tmp_path / 'not-there', chunk, 'HLA-A02:01')


def test_a_shell_wrapped_netmhcpan_works(wrapped_netmhcpan, chunk, fake_netmhcpan, monkeypatch):
    monkeypatch.setenv('FAKE_NETMHCPAN', str(fake_netmhcpan))
    rows = run_chunk(wrapped_netmhcpan, chunk, 'HLA-A02:01')
    assert len(rows) == 5


def test_killing_a_shell_wrapped_netmhcpan_takes_its_children_too(wrapped_netmhcpan, chunk,
                                                                  fake_netmhcpan):
    # A failure part way through the output has to leave nothing behind holding
    # the pipes open. Killing only the wrapper would leave its sleeping child
    # with the stderr pipe, and the drain thread would never see EOF - so this
    # would hang for ever rather than fail.
    child = run_in_a_child(
        f"""
        import netmhcpan_multiprocessor.runner as runner
        from netmhcpan_multiprocessor.chunking import PeptideChunk
        from pathlib import Path

        def explode(line):
            raise RuntimeError('something went wrong part way through')

        runner.is_prediction_row = explode
        chunk = PeptideChunk(path=Path({str(chunk.path)!r}), peptide_count=5)
        try:
            runner.run_chunk(Path({str(wrapped_netmhcpan)!r}), chunk, 'HLA-A02:01')
        except RuntimeError:
            print('raised and returned')
        """,
        environment={'FAKE_NETMHCPAN': str(fake_netmhcpan)},
    )
    assert 'raised and returned' in child.stdout, child.stderr


def test_a_lot_of_stderr_does_not_deadlock(fake_netmhcpan, chunk, monkeypatch):
    # stdout and stderr are both pipes: reading one to EOF before touching the
    # other blocks as soon as the child fills the stderr buffer.
    monkeypatch.setenv('FAKE_NETMHCPAN_STDERR_BYTES', str(512 * 1024))
    assert len(run_chunk(fake_netmhcpan, chunk, 'HLA-A02:01')) == 5


def test_undecodable_output_does_not_lose_the_chunk(fake_netmhcpan, chunk, monkeypatch):
    # A stray byte in a warning should cost the warning, not the predictions.
    monkeypatch.setenv('FAKE_NETMHCPAN_STDERR_BYTES', '64')
    monkeypatch.setenv('FAKE_NETMHCPAN_STDERR_RAW', '1')
    assert len(run_chunk(fake_netmhcpan, chunk, 'HLA-A02:01')) == 5


@pytest.mark.parametrize('fail_mode, expected', [
    ('exit', 'exit 1'),
    ('error', 'not found'),
    ('truncate', 'expected 5 predictions, got 4'),
    ('shuffle', 'different order'),
])
def test_every_way_netmhcpan_fails_is_caught(fake_netmhcpan, chunk, monkeypatch,
                                             fail_mode, expected):
    monkeypatch.setenv('FAKE_NETMHCPAN_FAIL_MODE', fail_mode)
    with pytest.raises(NetMHCpanRunError, match=expected):
        run_chunk(fake_netmhcpan, chunk, 'HLA-A02:01')


def test_affinity_arguments_are_passed_through(fake_netmhcpan, chunk):
    rows = run_chunk(fake_netmhcpan, chunk, 'HLA-A02:01', arguments=('-BA',))
    # -BA widens every row from 13 fields to 16.
    assert len(rows[0].split()) in (16, 18)
