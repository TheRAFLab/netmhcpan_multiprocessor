"""
Allele naming.

Three forms of the same allele are in play, and each has exactly one job:

    allele number    HLA-A*02:01     what a human writes, and what is reported
    netMHCpan name   HLA-A02:01      what is passed to netMHCpan's -a option
    allele slug      hla_a_02_01     folder names, and the allele_slug column

A caller may supply either the number or the slug, mixed freely within one run,
so everything here is about resolving whatever was supplied onto all three.

The authoritative mapping is netMHCpan's own allele list, read once per run from
`-listMHC`. Resolving against the list rather than against string rules is what
makes non-HLA names such as H-2-Db and Mamu-A1:00101 work, and it turns a typo
into an error in the first second of a run rather than a failure forty minutes
into one. A rule-based fallback covers HLA class I when the list is unavailable.
"""

import difflib
import re
import subprocess
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from netmhcpan_multiprocessor.errors import NetMHCpanRunError, UnknownAlleleError


# The allele list writes class I HLA names without the field separator, e.g.
# 'HLA-A01:01'. Putting the '*' back gives the form people write and, through
# slugify(), the slug this package has always used for folder names. Anything
# that does not match this shape is left alone: H-2-Db and BoLA-1:00901 have no
# starred form, and inventing one would only break the lookup.
STARRABLE_ALLELE = re.compile(r'^([A-Za-z]+-[A-Z]+[0-9]?)\*?([0-9]{2,3}):([0-9]{2,5})$')

# The inverse, for the fallback path: hla_a_02_01 -> HLA-A*02:01. Deliberately
# limited to the HLA class I shape, because that is the only slug that can be
# read back unambiguously without consulting the allele list.
HLA_SLUG = re.compile(r'^hla_([a-z]{1,4}[0-9]?)_([0-9]{2,3})_([0-9]{2,5})$')


@dataclass(frozen=True)
class ResolvedAllele:
    """
    One allele in all three of the forms this package needs.

    Attributes:
        number (str): the allele number, e.g. 'HLA-A*02:01'.
        netmhcpan_name (str): the name netMHCpan's -a option takes, e.g. 'HLA-A02:01'.
        slug (str): the slug used for folder names and the allele_slug column.
    """
    number: str
    netmhcpan_name: str
    slug: str


def slugify(text: str) -> str:
    """
    Reduce a name to a lowercase, underscore-separated slug safe for a path.

    Args:
        text (str): the text to slugify.

    Returns:
        str: the slug, e.g. 'HLA-A*02:01' becomes 'hla_a_02_01'.
    """
    slug_char = '_'
    to_replace = [' ', '-', '.', ',', '[', ']', '{', '}', '(', ')', '/', '\\', '*', ':']
    for replace_me in to_replace:
        if replace_me in text:
            text = text.replace(replace_me, slug_char)
    while slug_char * 2 in text:
        text = text.replace(slug_char * 2, slug_char)
    text = text.strip(slug_char)
    return text.lower()


def to_netmhcpan_name(allele_number: str) -> str:
    """
    Convert an allele number to the name netMHCpan's -a option expects.

    Args:
        allele_number (str): the allele number, e.g. 'HLA-A*02:01'.

    Returns:
        str: the netMHCpan name, e.g. 'HLA-A02:01'.
    """
    return allele_number.replace('*', '')


def to_allele_number(netmhcpan_name: str) -> str:
    """
    Convert a netMHCpan allele name to the number form a human would write.

    Args:
        netmhcpan_name (str): a name as netMHCpan lists it, e.g. 'HLA-A02:01'.

    Returns:
        str: the allele number, e.g. 'HLA-A*02:01'. Names with no starred form,
            such as 'H-2-Db', are returned unchanged.
    """
    match = STARRABLE_ALLELE.match(netmhcpan_name)
    if not match:
        return netmhcpan_name
    gene, field_one, field_two = match.groups()
    return f"{gene}*{field_one}:{field_two}"


def list_netmhcpan_alleles(netmhcpan_path: Path | str) -> tuple[str, ...]:
    """
    Read the alleles a netMHCpan install knows about.

    Args:
        netmhcpan_path (Path | str): path to the netMHCpan executable.

    Returns:
        tuple[str, ...]: every allele name, in the order netMHCpan lists them.

    Raises:
        NetMHCpanRunError: if netMHCpan could not be run, or listed nothing.
    """
    return _list_netmhcpan_alleles_cached(str(netmhcpan_path))


@lru_cache(maxsize=8)
def _list_netmhcpan_alleles_cached(netmhcpan_path: str) -> tuple[str, ...]:
    # Cached on the executable path: the list runs to ~13,000 names and never
    # changes for a given install, so one subprocess per process is plenty.
    try:
        listing = subprocess.run(
            [netmhcpan_path, '-listMHC'],
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as listing_error:
        raise NetMHCpanRunError(f"Could not run {netmhcpan_path} -listMHC: {listing_error}") from listing_error

    if listing.returncode != 0:
        raise NetMHCpanRunError(
            f"{netmhcpan_path} -listMHC exited {listing.returncode}: {listing.stderr.strip()}"
        )
    # The listing is prefixed with a comment block describing the run, which
    # every netMHCpan invocation prints regardless of what was asked for.
    alleles = tuple(
        line.strip() for line in listing.stdout.splitlines()
        if line.strip() and not line.startswith('#')
    )
    if not alleles:
        raise NetMHCpanRunError(f"{netmhcpan_path} -listMHC listed no alleles")
    return alleles


def build_allele_index(netmhcpan_path: Path | str) -> dict[str, str]:
    """
    Build a lookup from allele slug to netMHCpan allele name.

    Args:
        netmhcpan_path (Path | str): path to the netMHCpan executable.

    Returns:
        dict[str, str]: every allele the install knows, keyed by slug. Each is
            indexed under the slug of its own listed name and under the slug of
            its starred form, so 'hla_a02_01' and 'hla_a_02_01' both resolve.

    Raises:
        NetMHCpanRunError: if the allele list could not be read.
    """
    index: dict[str, str] = {}
    for name in list_netmhcpan_alleles(netmhcpan_path):
        # setdefault, not assignment: where two listed names collapse onto one
        # slug the first in netMHCpan's own order wins, which keeps the mapping
        # stable from run to run.
        index.setdefault(slugify(name), name)
        index.setdefault(slugify(to_allele_number(name)), name)
    return index


def resolve_allele(value: str, allele_index: dict[str, str] | None = None) -> ResolvedAllele:
    """
    Resolve one supplied allele, given as a number or a slug, onto all three forms.

    Args:
        value (str): an allele number such as 'HLA-A*02:01', or a slug such as
            'hla_a_02_01'.
        allele_index (dict[str, str] | None): the lookup from build_allele_index.
            None resolves by string rules alone, which only covers HLA class I.

    Returns:
        ResolvedAllele: the allele in number, netMHCpan and slug form.

    Raises:
        UnknownAlleleError: if the value is empty, or names no allele netMHCpan knows.
    """
    value = value.strip()
    if not value:
        raise UnknownAlleleError("An allele was supplied as an empty string")

    # slugify() is idempotent and maps both input forms onto the same key, so one
    # lookup covers 'HLA-A*02:01', 'HLA-A02:01' and 'hla_a_02_01' alike.
    slug = slugify(value)
    if allele_index is not None:
        netmhcpan_name = allele_index.get(slug)
        if netmhcpan_name is None:
            raise UnknownAlleleError(_unknown_allele_message(value, slug, allele_index))
        return _resolved_from_netmhcpan_name(netmhcpan_name)

    return _resolve_without_index(value, slug)


def resolve_alleles(
    values: list[str] | tuple[str, ...],
    netmhcpan_path: Path | str | None = None,
    validate: bool = True,
) -> list[ResolvedAllele]:
    """
    Resolve every supplied allele, dropping duplicates and keeping the given order.

    Args:
        values (list[str] | tuple[str, ...]): allele numbers and/or slugs, mixed freely.
        netmhcpan_path (Path | str | None): path to the netMHCpan executable, needed
            when validating.
        validate (bool): check each allele against netMHCpan's own allele list.
            False resolves by string rules alone, which only covers HLA class I but
            needs no netMHCpan install.

    Returns:
        list[ResolvedAllele]: one entry per distinct allele, in the order supplied.

    Raises:
        UnknownAlleleError: if no alleles were supplied, or one cannot be resolved.
        NetMHCpanRunError: if validating and the allele list could not be read.
    """
    if not values:
        raise UnknownAlleleError("No alleles were supplied")

    allele_index = build_allele_index(netmhcpan_path) if validate and netmhcpan_path else None

    resolved: dict[str, ResolvedAllele] = {}
    for value in values:
        allele = resolve_allele(value, allele_index)
        # Keyed by slug so 'HLA-A*02:01' and 'hla_a_02_01' in one call are the
        # single run they obviously mean, rather than two runs writing the same
        # output file twice.
        resolved.setdefault(allele.slug, allele)
    return list(resolved.values())


def _resolved_from_netmhcpan_name(netmhcpan_name: str) -> ResolvedAllele:
    number = to_allele_number(netmhcpan_name)
    # The slug always comes from the number, never from the listed name, so an
    # allele reachable by two spellings still writes to exactly one folder.
    return ResolvedAllele(number=number, netmhcpan_name=netmhcpan_name, slug=slugify(number))


def _resolve_without_index(value: str, slug: str) -> ResolvedAllele:
    match = HLA_SLUG.match(slug)
    if match:
        gene, field_one, field_two = match.groups()
        number = f"HLA-{gene.upper()}*{field_one}:{field_two}"
        return ResolvedAllele(number=number, netmhcpan_name=to_netmhcpan_name(number), slug=slug)
    if STARRABLE_ALLELE.match(value):
        number = to_allele_number(value)
        return ResolvedAllele(number=number, netmhcpan_name=to_netmhcpan_name(number), slug=slugify(number))
    raise UnknownAlleleError(
        f"Cannot resolve allele {value!r} without netMHCpan's allele list. "
        "Rule-based resolution covers HLA class I only, e.g. 'HLA-A*02:01' or 'hla_a_02_01'. "
        "Supply a netMHCpan path with validation enabled to use any other allele."
    )


def _unknown_allele_message(value: str, slug: str, allele_index: dict[str, str]) -> str:
    suggestions = difflib.get_close_matches(slug, allele_index, n=5, cutoff=0.7)
    message = f"netMHCpan does not know the allele {value!r}"
    if suggestions:
        # Suggest in the form the user is most likely to want to type back.
        close = ', '.join(sorted({to_allele_number(allele_index[each]) for each in suggestions}))
        message = f"{message}. Did you mean: {close}?"
    return message
