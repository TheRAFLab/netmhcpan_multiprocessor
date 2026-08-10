"""Finding and reading a proteome peptide set."""

import polars as pl
import pytest

from netmhcpan_multiprocessor.errors import PeptideInputError
from netmhcpan_multiprocessor.modes.base import PeptideFile
from netmhcpan_multiprocessor.modes.proteome.source import ProteomeSource

from conftest import make_peptides, write_peptide_parquet


@pytest.fixture
def source():
    return ProteomeSource()


class TestDiscovery:
    def test_every_length_is_found_by_default(self, source, peptide_folder):
        found = source.discover(peptide_folder, None)
        assert [(each.length, each.dataset) for each in found] == [(8, 'testset'), (9, 'testset')]

    def test_lengths_can_be_restricted(self, source, peptide_folder):
        assert [each.length for each in source.discover(peptide_folder, (9,))] == [9]

    def test_results_are_ordered_by_length(self, source, tmp_path):
        for length in (11, 8, 10):
            write_peptide_parquet(tmp_path, length=length, count=3)
        assert [each.length for each in source.discover(tmp_path, None)] == [8, 10, 11]

    def test_several_datasets_in_one_length_folder(self, source, tmp_path):
        write_peptide_parquet(tmp_path, length=9, count=3, dataset='alpha')
        write_peptide_parquet(tmp_path, length=9, count=3, dataset='beta')
        found = source.discover(tmp_path, None)
        # Each keeps its own name, because that name goes on to name the output.
        assert [each.dataset for each in found] == ['alpha', 'beta']

    def test_a_requested_length_that_is_absent_is_an_error(self, source, peptide_folder):
        # Quietly predicting three of the four lengths asked for is worse than
        # stopping, because nothing downstream would notice.
        with pytest.raises(PeptideInputError) as raised:
            source.discover(peptide_folder, (9, 10))
        assert '10mers' in str(raised.value)
        assert 'Found: 8mers, 9mers' in str(raised.value)

    def test_a_missing_folder_is_an_error(self, source, tmp_path):
        with pytest.raises(PeptideInputError, match='No peptide input folder'):
            source.discover(tmp_path / 'absent', None)

    def test_a_folder_with_no_length_folders_is_an_error(self, source, tmp_path):
        (tmp_path / 'peptides.parquet').write_text('not really')
        with pytest.raises(PeptideInputError, match='<length>mers folders'):
            source.discover(tmp_path, None)

    def test_an_empty_length_folder_is_ignored(self, source, tmp_path):
        write_peptide_parquet(tmp_path, length=9, count=3)
        (tmp_path / '10mers').mkdir()
        assert [each.length for each in source.discover(tmp_path, None)] == [9]

    def test_unrelated_folders_are_ignored(self, source, tmp_path):
        write_peptide_parquet(tmp_path, length=9, count=3)
        (tmp_path / 'scratch').mkdir()
        (tmp_path / 'scratch' / 'other.parquet').write_text('x')
        assert len(source.discover(tmp_path, None)) == 1


class TestReading:
    def test_peptides_come_back_in_file_order(self, source, peptide_folder):
        peptide_file = source.discover(peptide_folder, (9,))[0]
        assert source.read_peptides(peptide_file, 'peptide_sequence') == make_peptides(9, 30)

    def test_a_missing_column_names_what_is_there(self, source, peptide_folder):
        peptide_file = source.discover(peptide_folder, (9,))[0]
        with pytest.raises(PeptideInputError) as raised:
            source.read_peptides(peptide_file, 'sequence')
        assert "no 'sequence' column" in str(raised.value)
        assert 'peptide_sequence' in str(raised.value)

    def test_peptides_of_the_wrong_length_are_caught(self, source, tmp_path):
        # netMHCpan only reports a length mismatch as a short result, which is a
        # long way to travel to find out the input was wrong.
        folder = tmp_path / '9mers'
        folder.mkdir(parents=True)
        path = folder / 'testset.parquet'
        pl.DataFrame({'peptide_sequence': ['MESRIWCLV', 'MESRIWCL']}).write_parquet(path)
        with pytest.raises(PeptideInputError, match='holds peptides of length 8-9'):
            source.read_peptides(PeptideFile(9, path, 'testset'), 'peptide_sequence')

    def test_a_uniformly_wrong_length_is_caught(self, source, tmp_path):
        folder = tmp_path / '9mers'
        folder.mkdir(parents=True)
        path = folder / 'testset.parquet'
        pl.DataFrame({'peptide_sequence': ['MESRIWCL', 'MESRIWCM']}).write_parquet(path)
        with pytest.raises(PeptideInputError, match='holds peptides of length 8'):
            source.read_peptides(PeptideFile(9, path, 'testset'), 'peptide_sequence')

    def test_nulls_are_caught(self, source, tmp_path):
        folder = tmp_path / '9mers'
        folder.mkdir(parents=True)
        path = folder / 'testset.parquet'
        pl.DataFrame({'peptide_sequence': ['MESRIWCLV', None]}).write_parquet(path)
        with pytest.raises(PeptideInputError, match='null peptides'):
            source.read_peptides(PeptideFile(9, path, 'testset'), 'peptide_sequence')

    def test_an_empty_file_is_caught(self, source, tmp_path):
        folder = tmp_path / '9mers'
        folder.mkdir(parents=True)
        path = folder / 'testset.parquet'
        pl.DataFrame({'peptide_sequence': []}, schema={'peptide_sequence': pl.String}).write_parquet(path)
        with pytest.raises(PeptideInputError, match='holds no peptides'):
            source.read_peptides(PeptideFile(9, path, 'testset'), 'peptide_sequence')
