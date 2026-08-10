"""Splitting a peptide set into chunk files."""

import pytest

from netmhcpan_multiprocessor.chunking import write_chunk_files
from netmhcpan_multiprocessor.errors import ConfigError

from conftest import make_peptides


def read_chunks(chunks):
    return [chunk.path.read_text().splitlines() for chunk in chunks]


def test_an_exact_multiple(tmp_path):
    peptides = make_peptides(9, 20)
    chunks = write_chunk_files(peptides, chunk_size=10, chunk_folder=tmp_path, name='testset')
    assert [chunk.peptide_count for chunk in chunks] == [10, 10]


def test_a_remainder(tmp_path):
    peptides = make_peptides(9, 25)
    chunks = write_chunk_files(peptides, chunk_size=10, chunk_folder=tmp_path, name='testset')
    assert [chunk.peptide_count for chunk in chunks] == [10, 10, 5]


def test_one_chunk_when_the_set_is_small(tmp_path):
    chunks = write_chunk_files(make_peptides(9, 3), 5000, tmp_path, 'testset')
    assert len(chunks) == 1
    assert chunks[0].peptide_count == 3


def test_chunks_reassemble_into_the_input_order(tmp_path):
    # This is the whole contract: reassembling in chunk order has to give the
    # peptide file back exactly, or the output parquet will not line up with it.
    peptides = make_peptides(9, 25)
    chunks = write_chunk_files(peptides, chunk_size=7, chunk_folder=tmp_path, name='testset')
    assert [peptide for chunk in read_chunks(chunks) for peptide in chunk] == peptides


def test_chunk_files_are_one_peptide_per_line(tmp_path):
    chunks = write_chunk_files(make_peptides(9, 4), 2, tmp_path, 'testset')
    assert read_chunks(chunks)[0] == make_peptides(9, 4)[:2]


def test_chunk_files_are_named_for_their_set(tmp_path):
    chunks = write_chunk_files(make_peptides(9, 4), 2, tmp_path, 'testset__9mers')
    assert chunks[0].path.name == 'testset__9mers__chunk0000.txt'
    assert chunks[1].path.name == 'testset__9mers__chunk0001.txt'


def test_an_empty_peptide_list_makes_no_chunks(tmp_path):
    assert write_chunk_files([], 10, tmp_path, 'testset') == []


def test_a_useless_chunk_size_is_rejected(tmp_path):
    with pytest.raises(ConfigError, match='at least 1'):
        write_chunk_files(make_peptides(9, 4), 0, tmp_path, 'testset')
