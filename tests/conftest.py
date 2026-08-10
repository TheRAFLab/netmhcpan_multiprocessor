"""
Shared test fixtures.

The suite runs against the stand-in netMHCpan in fixtures/, so it needs no
netMHCpan install. Tests that do want the real thing are marked `netmhcpan` and
skipped unless $NETMHCPAN points at one.
"""

import os
import stat
from pathlib import Path

import polars as pl
import pytest


AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'

FIXTURES = Path(__file__).parent / 'fixtures'


def peptide_at(index: int, length: int) -> str:
    """
    Build the index-th distinct peptide of a given length.

    Args:
        index (int): which peptide to build.
        length (int): how long it should be.

    Returns:
        str: a peptide of standard amino acids, distinct for each index. The
            first residue varies fastest, so a short run still covers the S and
            W that the stand-in netMHCpan treats as strong and weak binders.
    """
    residues = []
    remaining = index
    for _ in range(length):
        residues.append(AMINO_ACIDS[remaining % len(AMINO_ACIDS)])
        remaining //= len(AMINO_ACIDS)
    return ''.join(residues)


def make_peptides(length: int, count: int) -> list[str]:
    """
    Build a list of distinct peptides.

    Args:
        length (int): the peptide length.
        count (int): how many to build.

    Returns:
        list[str]: the peptides, in a stable order.
    """
    return [peptide_at(index, length) for index in range(count)]


def write_peptide_parquet(folder: Path, length: int, count: int, dataset: str = 'testset') -> Path:
    """
    Write one peptide parquet in the layout the proteome source expects.

    Args:
        folder (Path): the input folder to write a `<length>mers` folder into.
        length (int): the peptide length.
        count (int): how many peptides to write.
        dataset (str): the file name, which becomes the dataset name.

    Returns:
        Path: the parquet written.
    """
    peptides = make_peptides(length, count)
    length_folder = folder / f"{length}mers"
    length_folder.mkdir(parents=True, exist_ok=True)
    path = length_folder / f"{dataset}.parquet"
    pl.DataFrame({
        'peptide_sequence': peptides,
        'id': [f"P{index:05d}" for index in range(count)],
        'start': list(range(1, count + 1)),
        'end': list(range(length, count + length)),
    }).write_parquet(path)
    return path


@pytest.fixture(scope='session')
def fake_netmhcpan() -> Path:
    """
    Path: an executable stand-in for netMHCpan.

    The executable bit is set here rather than relied on from the checkout, so
    the suite works however the files arrived.
    """
    stub = FIXTURES / 'fake_netmhcpan.py'
    stub.chmod(stub.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return stub


@pytest.fixture
def peptide_folder(tmp_path: Path) -> Path:
    """Path: an input folder holding 8mers and 9mers for one dataset."""
    folder = tmp_path / 'peptides'
    write_peptide_parquet(folder, length=8, count=30)
    write_peptide_parquet(folder, length=9, count=30)
    return folder


@pytest.fixture
def output_folder(tmp_path: Path) -> Path:
    """Path: where a prediction run should write."""
    return tmp_path / 'predictions'


@pytest.fixture
def prediction_settings(peptide_folder, output_folder, fake_netmhcpan) -> dict:
    """dict: a config dictionary that predicts the whole fixture set."""
    return {
        'input_folder': peptide_folder,
        'output_folder': output_folder,
        'alleles': ['HLA-A*02:01'],
        'netmhcpan_path': fake_netmhcpan,
        'chunk_size': 10,
        'workers': 4,
    }


@pytest.fixture
def real_netmhcpan() -> Path:
    """Path: a real netMHCpan, or skip the test."""
    configured = os.environ.get('NETMHCPAN')
    if not configured or not Path(configured).is_file():
        pytest.skip("set NETMHCPAN to a real netMHCpan executable to run this test")
    return Path(configured)
