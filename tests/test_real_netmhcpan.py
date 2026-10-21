"""
The same pipeline against a real netMHCpan.

Skipped unless $NETMHCPAN points at one, so the suite never depends on a licensed
install, but the assumptions the stand-in encodes can still be checked against
the real program on a machine that has it.
"""

import polars as pl
import pytest

from netmhcpan_multiprocessor import predict
from netmhcpan_multiprocessor.alleles import build_allele_index, resolve_alleles

pytestmark = pytest.mark.netmhcpan


# Well-known epitopes, so a real run has something to find: SLYNTVATL is an
# HLA-A*02:01 strong binder and NLVPMVATV is one of the best characterised.
PEPTIDES = ['SLYNTVATL', 'NLVPMVATV', 'MESRIWCLV', 'AAAAAAAAA']


@pytest.fixture
def peptide_set(tmp_path):
    folder = tmp_path / 'peptides' / '9mers'
    folder.mkdir(parents=True)
    pl.DataFrame({'peptide_sequence': PEPTIDES}).write_parquet(folder / 'testset.parquet')
    return tmp_path / 'peptides'


@pytest.fixture
def settings(peptide_set, tmp_path, real_netmhcpan):
    return {
        'input_folder': peptide_set,
        'output_folder': tmp_path / 'predictions',
        'alleles': ['HLA-A*02:01'],
        'netmhcpan_path': real_netmhcpan,
    }


def test_the_allele_list_covers_both_spellings(real_netmhcpan):
    index = build_allele_index(real_netmhcpan)
    assert index['hla_a_02_01'] == 'HLA-A02:01'
    assert index['h_2_db'] == 'H-2-Db'


def test_alleles_resolve_from_numbers_and_slugs(real_netmhcpan):
    resolved = resolve_alleles(['HLA-A*02:01', 'hla_b_07_02'], real_netmhcpan)
    assert [allele.netmhcpan_name for allele in resolved] == ['HLA-A02:01', 'HLA-B07:02']


def test_a_real_run_writes_predictions_in_input_order(settings, tmp_path):
    result = predict(settings)
    assert result.succeeded, [run.error for run in result.runs]

    predictions = pl.read_parquet(
        tmp_path / 'predictions' / 'hla_a_02_01' / '9mers' / 'testset.parquet')
    assert predictions['peptide_sequence'].to_list() == PEPTIDES


def test_a_real_run_finds_the_known_binders(settings, tmp_path):
    predict(settings)
    predictions = pl.read_parquet(
        tmp_path / 'predictions' / 'hla_a_02_01' / '9mers' / 'testset.parquet')
    binders = predictions.filter(pl.col('bind_level').is_not_null())['peptide_sequence'].to_list()
    assert 'SLYNTVATL' in binders
    assert 'NLVPMVATV' in binders


def test_chunking_a_real_run_preserves_order(settings, tmp_path):
    predict({**settings, 'chunk_size': 1, 'workers': 4})
    predictions = pl.read_parquet(
        tmp_path / 'predictions' / 'hla_a_02_01' / '9mers' / 'testset.parquet')
    assert predictions['peptide_sequence'].to_list() == PEPTIDES


def test_a_real_affinity_run_parses(settings, tmp_path):
    # The one assumption most likely to rot: how many columns -BA adds, and where.
    result = predict({**settings, 'scores': 'ba'})
    assert result.succeeded, [run.error for run in result.runs]

    predictions = pl.read_parquet(
        tmp_path / 'predictions' / 'hla_a_02_01' / '9mers' / 'testset.parquet')
    assert predictions['affinity_nm'].null_count() == 0
    assert predictions['score_el'].null_count() == 0
    strong = predictions.filter(pl.col('peptide_sequence') == 'SLYNTVATL')
    assert strong['bind_level'][0] == 'SB'
    assert strong['affinity_nm'][0] < 500


def test_a_real_unusable_allele_is_caught(settings):
    # netMHCpan prints 'Error.' and exits 0 for an allele it cannot use, so
    # validation is turned off here to let the request reach it.
    result = predict({**settings, 'alleles': ['HLA-A*99:99'], 'validate_alleles': False})
    assert not result.succeeded


PANEL = ['HLA-A*02:01', 'HLA-B*07:02', 'HLA-C*07:02']


def test_a_real_batch_matches_one_allele_per_process(settings, tmp_path):
    # The assumption batching rests on: netMHCpan given several alleles answers
    # with a block per allele, and those blocks hold exactly what the same
    # alleles return one process at a time.
    separate = tmp_path / 'separate'
    batched = tmp_path / 'batched'
    assert predict({**settings, 'alleles': PANEL, 'output_folder': separate}).succeeded
    assert predict({
        **settings, 'alleles': PANEL, 'output_folder': batched, 'alleles_per_process': 0,
    }).succeeded

    for allele_slug in ('hla_a_02_01', 'hla_b_07_02', 'hla_c_07_02'):
        path = allele_slug + '/9mers/testset.parquet'
        assert pl.read_parquet(batched / path).equals(pl.read_parquet(separate / path))


def test_a_real_batch_keeps_each_allele_in_input_order(settings, tmp_path):
    # Order is what the positional join back to the peptide parquet rests on, and
    # a batch is where it could go wrong: the peptides come back once per allele.
    assert predict({
        **settings, 'alleles': PANEL, 'alleles_per_process': 0, 'chunk_size': 2,
    }).succeeded
    for allele_slug in ('hla_a_02_01', 'hla_b_07_02', 'hla_c_07_02'):
        predictions = pl.read_parquet(
            tmp_path / 'predictions' / allele_slug / '9mers' / 'testset.parquet')
        assert predictions['peptide_sequence'].to_list() == PEPTIDES
        assert predictions['allele_slug'].to_list() == [allele_slug] * len(PEPTIDES)
