"""
Reading netMHCpan's output.

netMHCpan writes a fixed-width table wrapped in a banner and a summary, so
everything here is about telling the prediction rows apart from the decoration
and turning them into a typed frame.

Two score sets are supported. By default netMHCpan predicts eluted ligand
likelihood only; -BA adds binding affinity, which widens every row by three
columns. Which one is in use has to be known before a row can be read, since a
16 field affinity row and a 15 field eluted-ligand binder row are told apart by
what is in the columns, not by how many there are.
"""

from dataclasses import dataclass

import polars as pl

from netmhcpan_multiprocessor.alleles import slugify
from netmhcpan_multiprocessor.errors import ConfigError, NetMHCpanRunError


# netMHCpan's EL output columns, in the order they appear on a prediction row.
# BindLevel is the trailing '<= SB' / '<= WB' marker, absent on non-binders.
# This describes what is parsed, not what is written: the written frame swaps
# the allele column for the slug and drops identity. Of, Gp, Gl, Ip and Il give
# the peptide's alignment to the 9mer binding core, so they are all zero for
# 9mers but carry real information at every other length.
EL_ROW_SCHEMA = {
    'pos': pl.Int64,
    'allele': pl.String,
    'peptide': pl.String,
    'core': pl.String,
    'of': pl.Int64,
    'gp': pl.Int64,
    'gl': pl.Int64,
    'ip': pl.Int64,
    'il': pl.Int64,
    'icore': pl.String,
    'identity': pl.String,
    'score_el': pl.Float64,
    'rank_el': pl.Float64,
    'bind_level': pl.String,
}

# -BA inserts its three columns between %Rank_EL and BindLevel: the affinity
# score, its percentile rank, and the affinity itself in nanomolar.
BA_ROW_SCHEMA = {
    **{name: dtype for name, dtype in EL_ROW_SCHEMA.items() if name != 'bind_level'},
    'score_ba': pl.Float64,
    'rank_ba': pl.Float64,
    'affinity_nm': pl.Float64,
    'bind_level': pl.String,
}

# The columns written to parquet, in order.
EL_COLUMNS = (
    'allele_slug', 'peptide_sequence',
    'score_el', 'rank_el', 'bind_level',
    'core', 'of', 'gp', 'gl', 'ip', 'il', 'icore',
)
BA_COLUMNS = (
    'allele_slug', 'peptide_sequence',
    'score_el', 'rank_el', 'score_ba', 'rank_ba', 'affinity_nm', 'bind_level',
    'core', 'of', 'gp', 'gl', 'ip', 'il', 'icore',
)


@dataclass(frozen=True)
class ScoreSet:
    """
    One of netMHCpan's prediction outputs, and how to read it.

    Attributes:
        name (str): the name used in configuration, 'el' or 'ba'.
        arguments (tuple[str, ...]): the netMHCpan arguments that ask for it.
        row_schema (dict): the parsed row's columns and types, in row order.
        columns (tuple[str, ...]): the columns written to parquet, in order.
        field_count (int): how many whitespace-separated fields a row has before
            the optional '<= SB' marker.
        description (str): a one-line description, for --help.
    """
    name: str
    arguments: tuple[str, ...]
    row_schema: dict
    columns: tuple[str, ...]
    field_count: int
    description: str


EL = ScoreSet(
    name='el',
    arguments=(),
    row_schema=EL_ROW_SCHEMA,
    columns=EL_COLUMNS,
    field_count=13,
    description='eluted ligand likelihood only (netMHCpan default)',
)

BA = ScoreSet(
    name='ba',
    arguments=('-BA',),
    row_schema=BA_ROW_SCHEMA,
    columns=BA_COLUMNS,
    field_count=16,
    # netMHCpan does not predict affinity instead of eluted ligand likelihood,
    # it predicts it as well, so this is the wider output rather than the other
    # one. Said plainly here because 'BA' on its own suggests a swap.
    description='binding affinity as well as eluted ligand likelihood (netMHCpan -BA)',
)

SCORE_SETS = {score_set.name: score_set for score_set in (EL, BA)}

DEFAULT_SCORES = EL.name

# Kept for callers that want the default output's columns without resolving a
# score set first.
PREDICTION_COLUMNS = EL_COLUMNS


def get_score_set(scores: str | ScoreSet) -> ScoreSet:
    """
    Look a score set up by name.

    Args:
        scores (str | ScoreSet): 'el', 'ba', or a ScoreSet, returned as it is.

    Returns:
        ScoreSet: the score set to predict and parse with.

    Raises:
        ConfigError: if there is no score set of that name.
    """
    if isinstance(scores, ScoreSet):
        return scores
    try:
        return SCORE_SETS[str(scores).lower()]
    except KeyError:
        raise ConfigError(
            f"Unknown scores {scores!r}. Available: "
            + ', '.join(f"{name} ({each.description})" for name, each in SCORE_SETS.items())
        ) from None


def is_prediction_row(line: str) -> bool:
    """
    Test whether a line of netMHCpan output is a prediction, rather than part of
    the banner, the column headings or the trailing summary.

    Args:
        line (str): a single line of netMHCpan output.

    Returns:
        bool: True if the line is a prediction row.
    """
    # Prediction rows start with the position number and nothing else does. Do
    # not test for a leading space: the position column is right aligned in four
    # characters, so it loses its indent once the position reaches 1000.
    fields = line.split(maxsplit=1)
    return bool(fields) and fields[0].isdigit()


def split_rows_by_allele(prediction_rows: list[str]) -> dict[str, list[str]]:
    """
    Group prediction rows by the allele netMHCpan names in each one.

    Args:
        prediction_rows (list[str]): prediction rows, as a run against one or
            more alleles returns them.

    Returns:
        dict[str, list[str]]: the rows for each allele, keyed by its slug and
            holding that allele's rows in the order they were returned, which is
            the order the peptides were given. Keyed by slug rather than by the
            printed name because the name is netMHCpan's own spelling, and the
            slug is what the rest of this package identifies an allele by.

    Raises:
        NetMHCpanRunError: if a row has no allele column to group on.
    """
    rows_by_allele: dict[str, list[str]] = {}
    for row in prediction_rows:
        fields = row.split(maxsplit=2)
        if len(fields) < 2:
            raise NetMHCpanRunError(f"netMHCpan row has no allele column: {row.strip()!r}")
        # The rows print the allele in its starred form, which is the form the
        # slug is built from, so one slugify lands on the same slug a resolved
        # allele carries.
        rows_by_allele.setdefault(slugify(fields[1]), []).append(row)
    return rows_by_allele


def build_predictions_frame(
    prediction_rows: list[str],
    allele_slug: str,
    scores: str | ScoreSet = DEFAULT_SCORES,
) -> pl.DataFrame:
    """
    Parse merged prediction rows into a typed frame, preserving row order.

    Args:
        prediction_rows (list[str]): prediction rows in peptide input order.
        allele_slug (str): the slugified allele, stored in place of netMHCpan's
            own allele column.
        scores (str | ScoreSet): the score set the rows were produced with.

    Returns:
        pl.DataFrame: one row per prediction, in the order the peptides were
            given to netMHCpan, so it lines up row for row with the input
            peptide parquet.

    Raises:
        ConfigError: if the score set is not one this package knows.
        NetMHCpanRunError: if a row has an unexpected number of fields, or is
            not in the score set's shape.
    """
    score_set = get_score_set(scores)
    width = score_set.field_count

    records = []
    for position, row in enumerate(prediction_rows, start=1):
        fields = row.split()
        # The score set's own width, plus '<=' and the binder level on strong
        # and weak binders. The '<=' is checked rather than assumed, because an
        # affinity row is the same width as an eluted-ligand binder row and
        # would otherwise be read with a binder level of '12.3'.
        if len(fields) == width:
            bind_level = None
        elif len(fields) == width + 2 and fields[width] == '<=':
            bind_level = fields[width + 1]
        else:
            raise NetMHCpanRunError(
                f"Unexpected netMHCpan row with {len(fields)} fields, "
                f"reading {score_set.name.upper()} output: {row.strip()!r}. "
                f"Predictions made with scores='{'ba' if score_set is EL else 'el'}' "
                "have to be read back the same way."
            )
        records.append((position, *fields[1:width], bind_level))

    # orient='row' because records are positional tuples matching the schema.
    predictions = pl.DataFrame(records, schema=score_set.row_schema, orient='row')
    return predictions.with_columns(
        # identity is always 'PEPLIST' in -p mode, so it carries nothing, and pos
        # is just the row's own index. The allele is better held as the slug that
        # names the file it lives in. peptide_sequence rather than peptide to
        # match what the peptide parquet holds, so the two frames join directly.
        pl.lit(allele_slug, dtype=pl.String).alias('allele_slug'),
        pl.col('peptide').alias('peptide_sequence'),
    ).select(score_set.columns)


def write_raw_predictions(prediction_rows: list[str], output_file) -> None:
    """
    Write every chunk's prediction rows out as one file of predictions only.

    Args:
        prediction_rows (list[str]): prediction rows in peptide input order.
        output_file (Path | str): the path to write.

    Returns:
        None
    """
    with open(output_file, 'w') as output_handle:
        for position, row in enumerate(prediction_rows, start=1):
            # Every chunk numbers its own peptides from 1, so renumber to keep
            # the position column unique across the merged file. It is only ever
            # the index into the peptide list, so renumbering loses nothing.
            _, rest_of_row = row.split(maxsplit=1)
            output_handle.write(f"{position:>4}    {rest_of_row}")
