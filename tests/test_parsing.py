"""Reading netMHCpan's output rows."""

import polars as pl
import pytest

from netmhcpan_multiprocessor.errors import ConfigError, NetMHCpanRunError
from netmhcpan_multiprocessor.parsing import (
    DEFAULT_SCORES,
    build_predictions_frame,
    get_score_set,
    split_rows_by_allele,
    is_prediction_row,
    write_raw_predictions,
)


NON_BINDER = (
    "   1    HLA-A*02:01         MESRIWCLV MESRIWCLV  0  0  0  0  0"
    "         MESRIWCLV         PEPLIST 0.000160   47.644\n"
)
STRONG_BINDER = (
    "   2    HLA-A*02:01         SLYNTVATL SLYNTVATL  0  0  0  0  0"
    "         SLYNTVATL         PEPLIST 0.912345    0.041 <= SB\n"
)
WEAK_BINDER = (
    "   3    HLA-A*02:01         WMDSTPPPG WMDSTPPPG  0  0  0  0  0"
    "         WMDSTPPPG         PEPLIST 0.123456    1.500 <= WB\n"
)


@pytest.mark.parametrize('line', [NON_BINDER, STRONG_BINDER, WEAK_BINDER])
def test_prediction_rows_are_recognised(line):
    assert is_prediction_row(line)


def test_a_position_past_a_thousand_is_still_a_prediction():
    # The position column is right aligned in four characters, so it loses its
    # leading space at 1000. Anything testing for the indent breaks here.
    assert is_prediction_row("1000    HLA-A*02:01         MESRIWCLV" + NON_BINDER[35:])


@pytest.mark.parametrize('line', [
    '# NetMHCpan version 4.2\n',
    '\n',
    '-----------------------------------------\n',
    ' Pos          MHC        Peptide  Core Of Gp Gl Ip Il  Icore  Identity  Score_EL %Rank_EL\n',
    'Protein PEPLIST. Allele HLA-A02:01. Number of high binders 1. Number of peptides 3\n',
    'Error. MHC molecule HLA-A99:99 not found.\n',
])
def test_everything_else_is_not_a_prediction(line):
    assert not is_prediction_row(line)


def test_the_frame_holds_one_row_per_prediction_in_order():
    predictions = build_predictions_frame([NON_BINDER, STRONG_BINDER, WEAK_BINDER], 'hla_a_02_01')
    assert predictions.height == 3
    assert predictions['peptide_sequence'].to_list() == ['MESRIWCLV', 'SLYNTVATL', 'WMDSTPPPG']


def test_the_allele_is_stored_as_the_slug():
    predictions = build_predictions_frame([NON_BINDER], 'hla_a_02_01')
    assert predictions['allele_slug'].to_list() == ['hla_a_02_01']
    # netMHCpan's own allele column and its PEPLIST identity carry nothing the
    # slug and the file name do not, so neither is kept.
    assert 'allele' not in predictions.columns
    assert 'identity' not in predictions.columns


def test_bind_level_is_null_for_non_binders():
    predictions = build_predictions_frame([NON_BINDER, STRONG_BINDER, WEAK_BINDER], 'hla_a_02_01')
    assert predictions['bind_level'].to_list() == [None, 'SB', 'WB']


def test_the_schema_is_typed():
    predictions = build_predictions_frame([STRONG_BINDER], 'hla_a_02_01')
    assert predictions.schema['score_el'] == pl.Float64
    assert predictions.schema['rank_el'] == pl.Float64
    assert predictions.schema['of'] == pl.Int64
    assert predictions.schema['peptide_sequence'] == pl.String
    assert predictions.columns == [
        'allele_slug', 'peptide_sequence', 'score_el', 'rank_el', 'bind_level',
        'core', 'of', 'gp', 'gl', 'ip', 'il', 'icore',
    ]


def test_an_unexpected_row_shape_is_rejected():
    with pytest.raises(NetMHCpanRunError, match='Unexpected netMHCpan row'):
        build_predictions_frame([NON_BINDER.strip() + ' extra\n'], 'hla_a_02_01')


def test_a_binder_row_is_told_apart_by_its_marker_not_its_width():
    # Any 15 field row that is not a binder row has to be refused rather than
    # read, because the only thing separating the two is the '<='.
    fifteen_fields = NON_BINDER.strip() + ' 12345.67 42.000\n'
    assert len(fifteen_fields.split()) == len(STRONG_BINDER.split())
    with pytest.raises(NetMHCpanRunError, match='Unexpected netMHCpan row'):
        build_predictions_frame([fifteen_fields], 'hla_a_02_01')


def test_raw_output_renumbers_across_chunks(tmp_path):
    # Every chunk numbers its own peptides from 1, so a merged file would repeat
    # positions if they were not rewritten.
    raw_file = tmp_path / 'merged.txt'
    write_raw_predictions([NON_BINDER, NON_BINDER, NON_BINDER], raw_file)
    positions = [line.split()[0] for line in raw_file.read_text().splitlines()]
    assert positions == ['1', '2', '3']


def test_raw_output_keeps_the_rest_of_the_row(tmp_path):
    raw_file = tmp_path / 'merged.txt'
    write_raw_predictions([STRONG_BINDER], raw_file)
    assert raw_file.read_text().split(maxsplit=1)[1] == STRONG_BINDER.split(maxsplit=1)[1]


AFFINITY_NON_BINDER = (
    "   1    HLA-A*02:01         MESRIWCLV MESRIWCLV  0  0  0  0  0"
    "         MESRIWCLV         PEPLIST 0.000087   41.151 0.185539    9.331  6716.24\n"
)
AFFINITY_STRONG_BINDER = (
    "   2    HLA-A*02:01         SLYNTVATL SLYNTVATL  0  0  0  0  0"
    "         SLYNTVATL         PEPLIST 0.852778    0.083 0.599970    0.655    75.81 <= SB\n"
)


class TestScoreSets:
    def test_the_default_is_eluted_ligand(self):
        assert get_score_set(DEFAULT_SCORES).name == 'el'
        assert get_score_set('el').arguments == ()

    def test_affinity_asks_netmhcpan_for_it(self):
        # -BA adds affinity to the eluted ligand output rather than replacing it,
        # which is why the ba columns are a superset of the el ones.
        assert get_score_set('ba').arguments == ('-BA',)
        assert set(get_score_set('el').columns) < set(get_score_set('ba').columns)

    def test_an_unknown_score_set_is_rejected(self):
        with pytest.raises(ConfigError, match='Unknown scores'):
            get_score_set('affinity')

    def test_affinity_rows_parse_into_their_own_columns(self):
        predictions = build_predictions_frame(
            [AFFINITY_NON_BINDER, AFFINITY_STRONG_BINDER], 'hla_a_02_01', 'ba',
        )
        assert predictions.columns == list(get_score_set('ba').columns)
        assert predictions['affinity_nm'].to_list() == [6716.24, 75.81]
        assert predictions['rank_ba'].to_list() == [9.331, 0.655]
        assert predictions['bind_level'].to_list() == [None, 'SB']
        assert predictions.schema['affinity_nm'] == pl.Float64

    def test_affinity_rows_keep_the_eluted_ligand_columns_too(self):
        predictions = build_predictions_frame([AFFINITY_NON_BINDER], 'hla_a_02_01', 'ba')
        assert predictions['score_el'].to_list() == [0.000087]
        assert predictions['rank_el'].to_list() == [41.151]

    def test_reading_affinity_rows_as_eluted_ligand_is_refused(self):
        # The failure this guards against is silent: a 16 field affinity row read
        # as el would take '0.185539' for the '<=' marker.
        with pytest.raises(NetMHCpanRunError, match="scores='ba'"):
            build_predictions_frame([AFFINITY_NON_BINDER], 'hla_a_02_01', 'el')

    def test_reading_eluted_ligand_rows_as_affinity_is_refused(self):
        with pytest.raises(NetMHCpanRunError, match="scores='el'"):
            build_predictions_frame([NON_BINDER], 'hla_a_02_01', 'ba')


class TestSplittingABatch:
    def test_rows_are_grouped_by_allele_slug(self):
        other_allele = NON_BINDER.replace('HLA-A*02:01', 'HLA-B*07:02')
        rows_by_allele = split_rows_by_allele([NON_BINDER, STRONG_BINDER, other_allele])
        assert set(rows_by_allele) == {'hla_a_02_01', 'hla_b_07_02'}
        assert rows_by_allele['hla_a_02_01'] == [NON_BINDER, STRONG_BINDER]
        assert rows_by_allele['hla_b_07_02'] == [other_allele]

    def test_an_allele_keeps_the_order_its_rows_arrived_in(self):
        rows = [STRONG_BINDER, NON_BINDER, WEAK_BINDER]
        assert split_rows_by_allele(rows)['hla_a_02_01'] == rows

    def test_a_non_hla_allele_groups_under_its_own_slug(self):
        # H-2-Db and the like have no starred form, so the slug is the name.
        row = NON_BINDER.replace('HLA-A*02:01', 'H-2-Db')
        assert list(split_rows_by_allele([row])) == ['h_2_db']

    def test_a_row_with_no_allele_column_is_an_error(self):
        with pytest.raises(NetMHCpanRunError, match='no allele column'):
            split_rows_by_allele(['   1\n'])

    def test_a_split_block_parses_as_its_own_frame(self):
        other_allele = STRONG_BINDER.replace('HLA-A*02:01', 'HLA-B*07:02')
        rows_by_allele = split_rows_by_allele([STRONG_BINDER, other_allele])
        frame = build_predictions_frame(rows_by_allele['hla_b_07_02'], 'hla_b_07_02')
        assert frame['allele_slug'].to_list() == ['hla_b_07_02']
        assert frame['peptide_sequence'].to_list() == ['SLYNTVATL']
