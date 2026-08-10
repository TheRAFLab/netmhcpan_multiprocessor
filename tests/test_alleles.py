"""Allele naming: slugs, numbers, netMHCpan names, and resolving between them."""

import pytest

from netmhcpan_multiprocessor.alleles import (
    build_allele_index,
    list_netmhcpan_alleles,
    resolve_allele,
    resolve_alleles,
    slugify,
    to_allele_number,
    to_netmhcpan_name,
)
from netmhcpan_multiprocessor.errors import UnknownAlleleError


@pytest.fixture
def allele_index(fake_netmhcpan):
    return build_allele_index(fake_netmhcpan)


@pytest.mark.parametrize('text, expected', [
    ('HLA-A*02:01', 'hla_a_02_01'),
    ('HLA-A02:01', 'hla_a02_01'),
    ('hla_a_02_01', 'hla_a_02_01'),
    ('H-2-Db', 'h_2_db'),
    ('Mamu-A01', 'mamu_a01'),
    ('BoLA-1:00901', 'bola_1_00901'),
])
def test_slugify(text, expected):
    assert slugify(text) == expected


def test_slugify_is_idempotent():
    # The whole resolution scheme leans on this: a number and its own slug have
    # to land on the same lookup key.
    assert slugify(slugify('HLA-A*02:01')) == slugify('HLA-A*02:01')


@pytest.mark.parametrize('number, netmhcpan_name', [
    ('HLA-A*02:01', 'HLA-A02:01'),
    ('H-2-Db', 'H-2-Db'),
])
def test_to_netmhcpan_name(number, netmhcpan_name):
    assert to_netmhcpan_name(number) == netmhcpan_name


@pytest.mark.parametrize('netmhcpan_name, number', [
    ('HLA-A02:01', 'HLA-A*02:01'),
    ('HLA-B07:02', 'HLA-B*07:02'),
    ('HLA-A*02:01', 'HLA-A*02:01'),
    # No starred form exists for these, so inventing one would break the lookup.
    ('H-2-Db', 'H-2-Db'),
    ('Mamu-A01', 'Mamu-A01'),
    ('HLA-A0101', 'HLA-A0101'),
])
def test_to_allele_number(netmhcpan_name, number):
    assert to_allele_number(netmhcpan_name) == number


def test_listing_alleles_is_cached(fake_netmhcpan):
    first = list_netmhcpan_alleles(fake_netmhcpan)
    assert first is list_netmhcpan_alleles(fake_netmhcpan)
    assert 'HLA-A02:01' in first


@pytest.mark.parametrize('supplied', ['HLA-A*02:01', 'HLA-A02:01', 'hla_a_02_01', 'HLA-A*02:01  '])
def test_every_spelling_resolves_to_one_allele(supplied, allele_index):
    resolved = resolve_allele(supplied, allele_index)
    assert resolved.number == 'HLA-A*02:01'
    assert resolved.netmhcpan_name == 'HLA-A02:01'
    assert resolved.slug == 'hla_a_02_01'


@pytest.mark.parametrize('supplied', ['H-2-Db', 'h_2_db'])
def test_non_hla_alleles_resolve_through_the_index(supplied, allele_index):
    resolved = resolve_allele(supplied, allele_index)
    assert resolved.netmhcpan_name == 'H-2-Db'
    assert resolved.slug == 'h_2_db'


def test_unknown_allele_names_close_matches(allele_index):
    with pytest.raises(UnknownAlleleError) as raised:
        resolve_allele('HLA-A*02:02', allele_index)
    message = str(raised.value)
    assert 'HLA-A*02:02' in message
    assert 'HLA-A*02:01' in message


def test_unknown_allele_without_close_matches(allele_index):
    with pytest.raises(UnknownAlleleError, match='does not know'):
        resolve_allele('not-an-allele', allele_index)


def test_empty_allele_is_rejected(allele_index):
    with pytest.raises(UnknownAlleleError):
        resolve_allele('   ', allele_index)


def test_resolving_without_an_index_covers_hla_class_one():
    resolved = resolve_allele('hla_a_02_01', None)
    assert resolved.number == 'HLA-A*02:01'
    assert resolved.netmhcpan_name == 'HLA-A02:01'


def test_resolving_without_an_index_normalises_the_netmhcpan_spelling():
    # Both spellings have to produce the same slug with and without the index,
    # or the same job would write to two different folders.
    assert resolve_allele('HLA-A02:01', None).slug == 'hla_a_02_01'


def test_resolving_without_an_index_rejects_what_it_cannot_read():
    with pytest.raises(UnknownAlleleError, match='HLA class I only'):
        resolve_allele('h_2_db', None)


def test_resolve_alleles_keeps_order_and_drops_duplicates(fake_netmhcpan):
    resolved = resolve_alleles(
        ['HLA-B*07:02', 'HLA-A*02:01', 'hla_b_07_02'],
        netmhcpan_path=fake_netmhcpan,
    )
    assert [allele.slug for allele in resolved] == ['hla_b_07_02', 'hla_a_02_01']


def test_resolve_alleles_without_validation_needs_no_netmhcpan():
    resolved = resolve_alleles(['hla_a_02_01'], netmhcpan_path=None, validate=False)
    assert resolved[0].netmhcpan_name == 'HLA-A02:01'


def test_resolve_alleles_rejects_an_empty_list(fake_netmhcpan):
    with pytest.raises(UnknownAlleleError, match='No alleles'):
        resolve_alleles([], netmhcpan_path=fake_netmhcpan)


def test_the_index_holds_both_spellings(allele_index):
    # 'HLA-A01:01' and 'HLA-A0101' are separate entries in netMHCpan's own list,
    # and both have to remain reachable.
    assert allele_index['hla_a_01_01'] == 'HLA-A01:01'
    assert allele_index['hla_a01_01'] == 'HLA-A01:01'
    assert allele_index['hla_a0101'] == 'HLA-A0101'
