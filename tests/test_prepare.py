"""Building a proteome peptide set from a FASTA."""

import polars as pl
import pytest

from netmhcpan_multiprocessor import prepare
from netmhcpan_multiprocessor.config import PrepareConfig
from netmhcpan_multiprocessor.errors import PeptideInputError
from netmhcpan_multiprocessor.modes.proteome.source import ProteomeSource

pytest.importorskip('polars_bio', reason="preparation needs the 'proteome' extra")
pytest.importorskip('shared_helpers', reason="preparation needs the 'proteome' extra")


LONG_SEQUENCE = 'MESRIWCLVVCVNLTALWTL'          # 20 residues
SHORT_SEQUENCE = 'MKV'                          # too short for any length asked for
NONSTANDARD_SEQUENCE = 'MKVXQTPRLLIYDTSKLA'     # 18 residues, one of them an X


def write_fasta(tmp_path, records, name='testset.fasta'):
    """
    Write a UniProt-style FASTA.

    Args:
        tmp_path (Path): where to write it.
        records (list[tuple[str, str]]): accession and sequence pairs.
        name (str): the file name, which becomes the default dataset name.

    Returns:
        Path: the FASTA written.
    """
    fasta = tmp_path / name
    fasta.write_text(''.join(
        f">sp|{accession}|{accession}_TEST Test protein OS=Testus testus OX=1 GN=T PE=1 SV=1\n"
        f"{sequence}\n"
        for accession, sequence in records
    ))
    return fasta


@pytest.fixture
def fasta(tmp_path):
    return write_fasta(tmp_path, [('P00001', LONG_SEQUENCE), ('P00002', SHORT_SEQUENCE)])


def test_one_parquet_per_length(fasta, tmp_path):
    result = prepare({'input_path': fasta, 'output_folder': tmp_path / 'peptides', 'lengths': [8, 9]})
    assert [prepared.length for prepared in result.files] == [8, 9]
    assert all(prepared.path.exists() for prepared in result.files)


def test_the_layout_is_the_one_prediction_reads(fasta, tmp_path):
    # The seam between the two halves of the package: what prepare writes has to
    # be what the proteome source discovers, or they drift apart silently.
    output_folder = tmp_path / 'peptides'
    prepare({'input_path': fasta, 'output_folder': output_folder, 'lengths': [8, 9]})

    found = ProteomeSource().discover(output_folder, None)
    assert [(each.length, each.dataset) for each in found] == [(8, 'testset'), (9, 'testset')]
    assert ProteomeSource().read_peptides(found[0], 'peptide_sequence')


def test_every_kmer_of_every_long_enough_sequence(fasta, tmp_path):
    result = prepare({'input_path': fasta, 'output_folder': tmp_path / 'out', 'lengths': [9]})
    # 20 residues gives 12 nonamers; the 3 residue sequence gives none.
    assert result.files[0].peptide_count == len(LONG_SEQUENCE) - 9 + 1


def test_the_peptide_columns_are_what_prediction_expects(fasta, tmp_path):
    result = prepare({'input_path': fasta, 'output_folder': tmp_path / 'out', 'lengths': [9]})
    peptides = pl.read_parquet(result.files[0].path)
    assert peptides.columns == ['peptide_sequence', 'id', 'start', 'end']
    assert peptides['peptide_sequence'][0] == LONG_SEQUENCE[:9]
    assert peptides['id'][0] == 'P00001'


def test_coordinates_are_one_based_and_inclusive(fasta, tmp_path):
    result = prepare({'input_path': fasta, 'output_folder': tmp_path / 'out', 'lengths': [9]})
    peptides = pl.read_parquet(result.files[0].path)
    assert peptides['start'][0] == 1
    assert peptides['end'][0] == 9
    assert (peptides['end'] - peptides['start'] + 1).unique().to_list() == [9]


def test_nonstandard_residues_are_dropped(tmp_path):
    # netMHCpan cannot score a peptide containing an X, so leaving them in would
    # only fail a prediction run much later.
    fasta = write_fasta(tmp_path, [('P00003', NONSTANDARD_SEQUENCE)])
    result = prepare({'input_path': fasta, 'output_folder': tmp_path / 'out', 'lengths': [9]})
    peptides = pl.read_parquet(result.files[0].path)
    assert peptides.height < len(NONSTANDARD_SEQUENCE) - 9 + 1
    assert not peptides.filter(pl.col('peptide_sequence').str.contains('X')).height


def test_nonstandard_residues_can_be_kept(tmp_path):
    fasta = write_fasta(tmp_path, [('P00003', NONSTANDARD_SEQUENCE)])
    result = prepare({
        'input_path': fasta, 'output_folder': tmp_path / 'out',
        'lengths': [9], 'keep_nonstandard': True,
    })
    assert result.files[0].peptide_count == len(NONSTANDARD_SEQUENCE) - 9 + 1


def test_the_dataset_names_the_output_file(fasta, tmp_path):
    result = prepare({
        'input_path': fasta, 'output_folder': tmp_path / 'out',
        'lengths': [9], 'dataset': 'UP000008992',
    })
    assert result.files[0].path.name == 'UP000008992.parquet'
    assert result.dataset == 'UP000008992'


def test_a_length_nothing_is_long_enough_for_is_an_error(tmp_path):
    fasta = write_fasta(tmp_path, [('P00002', SHORT_SEQUENCE)])
    with pytest.raises(PeptideInputError, match='yielded no 9mers'):
        prepare({'input_path': fasta, 'output_folder': tmp_path / 'out', 'lengths': [9]})


def test_a_config_object_works_too(fasta, tmp_path):
    config = PrepareConfig.from_dict({
        'input_path': fasta, 'output_folder': tmp_path / 'out', 'lengths': [9],
    })
    assert prepare(config).files[0].peptide_count


def test_the_total_is_reported(fasta, tmp_path):
    result = prepare({'input_path': fasta, 'output_folder': tmp_path / 'out', 'lengths': [8, 9]})
    assert result.peptide_count == (len(LONG_SEQUENCE) - 8 + 1) + (len(LONG_SEQUENCE) - 9 + 1)


def test_nothing_is_printed_to_stdout(fasta, tmp_path, capfd):
    # stderr is not asserted on: polars-bio draws a tqdm progress bar there while
    # it reads, and tqdm fixes whether it is disabled when it is imported, so
    # this package cannot turn it off after the fact. The README says how.
    prepare({'input_path': fasta, 'output_folder': tmp_path / 'out', 'lengths': [9]})
    assert capfd.readouterr().out == ''
