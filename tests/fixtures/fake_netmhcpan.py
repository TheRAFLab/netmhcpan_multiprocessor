#!/usr/bin/env python3
"""
A stand-in for netMHCpan, so the whole pipeline can be tested without one.

It writes the same shape of output as the real program - the same banner, the
same column header, the same fixed-width prediction rows and the same trailing
summary - and it can be told, through the environment, to fail in each of the
ways netMHCpan actually fails. Only the standard library is used, so it runs
under whatever python3 the shebang finds.

It takes one or more alleles on -a, comma separated, and writes a block of rows
per allele as the real program does.

Environment:
    FAKE_NETMHCPAN_FAIL_ALLELE  only misbehave for this allele, so a run with
                                several alleles - whether as separate processes
                                or batched into one - can have exactly one go
                                wrong. Unset means every allele misbehaves.
    FAKE_NETMHCPAN_FAIL_MODE    'exit'     exit non-zero with a message on stderr
                                'error'    print 'Error.' to stdout and exit 0,
                                           which is what an unusable allele does
                                'truncate' return one prediction fewer than asked
                                'shuffle'  return every prediction, out of order
    FAKE_NETMHCPAN_SLEEP        seconds to sleep before writing anything.
    FAKE_NETMHCPAN_STDERR_BYTES roughly this many bytes of warnings on stderr,
                                written before stdout, as a real netMHCpan may.
                                More than a pipe buffer's worth deadlocks any
                                caller that reads stdout to EOF first.

Scores are derived from the peptide with crc32, not hash(), so they are the same
in every process and every run. Binder levels follow a rule a test can predict:
a peptide starting with S is a strong binder, one starting with W is a weak one.
"""

import os
import sys
import time
import zlib


ALLELE_LIST = (
    'HLA-A01:01',
    'HLA-A0101',
    'HLA-A02:01',
    'HLA-A02:06',
    'HLA-B07:02',
    'HLA-C07:02',
    'H-2-Db',
    'Mamu-A01',
)

HEADER_RULE = '-' * 100
COLUMN_HEADER = (
    ' Pos          MHC        Peptide   Core Of Gp Gl Ip Il        Icore '
    '       Identity  Score_EL %Rank_EL'
)
AFFINITY_HEADER = ' Score_BA %Rank_BA  Aff(nM)'


def parse_arguments(argv):
    """
    Read the -p, -a, -BA and -listMHC arguments, ignoring everything else.

    Args:
        argv (list[str]): the command line, without the program name.

    Returns:
        tuple[str | None, str | None, bool, bool]: peptide file, the -a value
            (which may name several alleles, comma separated), list mode, and
            whether binding affinity was asked for.
    """
    peptide_file, allele, list_alleles, affinity = None, None, False, False
    remaining = list(argv)
    while remaining:
        argument = remaining.pop(0)
        if argument == '-p':
            peptide_file = remaining.pop(0)
        elif argument == '-a':
            allele = remaining.pop(0)
        elif argument == '-listMHC':
            list_alleles = True
        elif argument == '-BA':
            affinity = True
    return peptide_file, allele, list_alleles, affinity


def banner(peptide_file, allele):
    """
    The comment block every netMHCpan run prints before its output.

    Args:
        peptide_file (str | None): the -p file.
        allele (str | None): the -a allele.

    Returns:
        list[str]: the banner lines.
    """
    return [
        '# Fake NetMHCpan version 4.2',
        f'# Input is a PEPTIDE list: {peptide_file}',
        f'# Allele: {allele}',
        '',
    ]


def starred(allele):
    """
    Put the field separator back, as netMHCpan does in its output rows.

    Args:
        allele (str): the allele as given to -a, e.g. 'HLA-A02:01'.

    Returns:
        str: the allele as printed, e.g. 'HLA-A*02:01'.
    """
    if allele.startswith('HLA-') and ':' in allele and '*' not in allele:
        gene, fields = allele[4:].split(':', 1)
        return f"HLA-{gene[:-2]}*{gene[-2:]}:{fields}"
    return allele


def prediction_row(position, allele, peptide, affinity):
    """
    One prediction row, in netMHCpan's own column layout.

    Args:
        position (int): the peptide's 1-based position in the input file.
        allele (str): the allele as given to -a.
        peptide (str): the peptide sequence.
        affinity (bool): whether -BA was asked for, which adds three columns
            between %Rank_EL and BindLevel.

    Returns:
        str: the row, without its newline.
    """
    fingerprint = zlib.crc32(peptide.encode())
    score = (fingerprint % 1_000_000) / 1_000_000
    rank = (fingerprint % 10_000) / 100

    # A rule a test can predict, rather than a threshold on the fake score.
    if peptide.startswith('S'):
        bind_level = '<= SB'
    elif peptide.startswith('W'):
        bind_level = '<= WB'
    else:
        bind_level = ''

    affinity_columns = ''
    if affinity:
        affinity_columns = (
            f"{score / 2:>9.6f} {rank * 1.5:>8.3f} {(1 - score) * 20000:>8.2f} "
        )

    core = peptide[:9]
    return (
        f"{position:>4} {starred(allele):>16} {peptide:>14} {core:>14} "
        f" 0  0  0  0  0 {peptide:>14} {'PEPLIST':>15} "
        f"{score:>9.6f} {rank:>8.3f} {affinity_columns}{bind_level}"
    ).rstrip()


def main(argv):
    """
    Run the stand-in.

    Args:
        argv (list[str]): the command line, without the program name.

    Returns:
        int: the exit code.
    """
    peptide_file, allele, list_alleles, affinity = parse_arguments(argv)

    if list_alleles:
        print('# Fake NetMHCpan version 4.2')
        for name in ALLELE_LIST:
            print(name)
        return 0

    sleep_for = float(os.environ.get('FAKE_NETMHCPAN_SLEEP', 0))
    if sleep_for:
        time.sleep(sleep_for)

    stderr_bytes = int(os.environ.get('FAKE_NETMHCPAN_STDERR_BYTES', 0))
    if stderr_bytes:
        warning = 'Warning: something a real netMHCpan might grumble about\n'
        if os.environ.get('FAKE_NETMHCPAN_STDERR_RAW'):
            # Bytes that are not valid UTF-8, as a program written in C and
            # printing whatever it was given may well emit.
            sys.stderr.buffer.write(b'Warning: \xff\xfe not text\n' * (stderr_bytes // 24 + 1))
            sys.stderr.buffer.flush()
        else:
            sys.stderr.write(warning * (stderr_bytes // len(warning) + 1))
            sys.stderr.flush()

    # -a takes a comma separated list, and the real program then writes one block
    # of rows per allele, in the order given, each block numbered from 1.
    alleles = [name for name in (allele or '').split(',') if name]

    fail_mode = os.environ.get('FAKE_NETMHCPAN_FAIL_MODE')
    fail_allele = os.environ.get('FAKE_NETMHCPAN_FAIL_ALLELE')
    misbehaving = bool(fail_mode) and (fail_allele in (None, '') or fail_allele in alleles)

    if misbehaving and fail_mode == 'exit':
        print(f"could not open {peptide_file}", file=sys.stderr)
        return 1
    if misbehaving and fail_mode == 'error':
        # An unusable allele: netMHCpan says so on stdout and then exits 0.
        for line in banner(peptide_file, allele):
            print(line)
        print(f'Error. MHC molecule {fail_allele or allele} not found.')
        return 0

    with open(peptide_file) as peptide_handle:
        peptides = [line.strip() for line in peptide_handle if line.strip()]

    for line in banner(peptide_file, allele):
        print(line)
    print(HEADER_RULE)
    print(COLUMN_HEADER + (AFFINITY_HEADER if affinity else '') + ' BindLevel')
    print(HEADER_RULE)
    for name in alleles:
        # Only the allele asked to misbehave does, so a run of several alleles
        # can have exactly one of them come back wrong.
        wrong = bool(fail_mode) and (fail_allele in (None, '', name))
        answered = list(peptides)
        if wrong and fail_mode == 'truncate':
            answered = answered[:-1]
        if wrong and fail_mode == 'shuffle':
            # Every peptide, answered in a different order than it was asked. A
            # reversal rather than a random shuffle, so the test is deterministic.
            answered = answered[::-1]
        for position, peptide in enumerate(answered, start=1):
            print(prediction_row(position, name, peptide, affinity))
    print(HEADER_RULE)
    print()
    strong = sum(1 for peptide in peptides if peptide.startswith('S'))
    weak = sum(1 for peptide in peptides if peptide.startswith('W'))
    for name in alleles:
        print(
            f"Protein PEPLIST. Allele {name}. Number of high binders {strong}. "
            f"Number of weak binders {weak}. Number of peptides {len(peptides)}"
        )
    print()
    print(HEADER_RULE)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
