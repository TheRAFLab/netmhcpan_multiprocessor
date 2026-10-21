# netMHCpan multiprocessor

Run netMHCpan over a set of peptides in parallel, from a terminal or from a Python
pipeline.

netMHCpan is single threaded, which is inefficient on modern multiprocessor machines, and
execution time for different peptide lengths is very uneven — 11mers take nearly ten times
as long as 9mers. This package cuts every peptide set into chunks, pools all the chunks
together, and hands them out longest first, so a worker that finishes a quick 9mer chunk
picks up an 11mer chunk instead of leaving a core idle behind the one slowest length.

It also cuts how much work there is to spread. Most of what netMHCpan does to a peptide —
enumerating how a long peptide aligns to the 9mer binding core — does not depend on the
allele, so `--alleles-per-process` scores a batch of alleles in one call. A 12mer chunk
against 32 alleles costs about a quarter of what those alleles cost a process each, and a
whole panel over a proteome comes down from days to hours.

Peptide parquet in, prediction parquet out, one file per length per allele, in the input's
row order so the two join directly.

```bash
uv run netmhcpan-mp predict \
    --input-folder peptides \
    --output-folder predictions \
    --allele 'HLA-A*02:01' --allele hla_b_07_02
```

```python
from netmhcpan_multiprocessor import predict

result = predict({
    'input_folder': 'peptides',
    'output_folder': 'predictions',
    'alleles': ['HLA-A*02:01', 'hla_b_07_02'],
})
```

Developed and tested against **netMHCpan 4.2**, which you install yourself — see
[Installing](#installing).

---

**New here?** [Installing](#installing) → [Running the commands](#running-the-commands) →
[End to end, from a FASTA](#end-to-end-from-a-fasta). If something does not work,
[When something goes wrong](#when-something-goes-wrong).

**Reference:** [Input and output](#input-and-output) · [Alleles](#alleles) ·
[Scores: EL and BA](#scores-el-and-ba) · [Modes and sources](#modes-and-sources) ·
[Preparing a peptide set](#preparing-a-peptide-set) · [Configuration](#configuration) ·
[Using it from a pipeline](#using-it-from-a-pipeline) · [Tests](#tests)

---

## Installing

You need three things: this repository, `uv` to install it with, and your own copy of
netMHCpan. Each step below prints something you can check against, so you can tell which
one went wrong rather than finding out later.

### 1. uv

`uv` installs Python and every dependency for you, so it is the only thing you have to
install by hand. If `uv --version` prints a version, skip this.

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh    # macOS and Linux
uv --version                                       # e.g. uv 0.11.31
```

Restart your terminal if `uv` is still not found after installing it.

### 2. This package

```bash
git clone https://github.com/TheRAFLab/netmhcpan_multiprocessor.git
cd netmhcpan_multiprocessor
uv sync --extra proteome
```

`uv sync` reads `pyproject.toml`, fetches the right Python and every dependency, and puts
them in a `.venv` folder inside the project. Nothing is installed system-wide, and nothing
outside this folder is touched. `--extra proteome` adds what `prepare` needs to read FASTA
files; leave it off if you only ever predict from peptide parquet that somebody else made.

Check it worked:

```bash
uv run netmhcpan-mp --version      # netmhcpan-mp, version 0.1.0
```

If you would rather manage the environment yourself, `python -m pip install -e '.[proteome]'`
inside your own virtual environment does the same job; everything below then works with the
bare `netmhcpan-mp` command.

### 3. netMHCpan

netMHCpan is licensed by DTU and cannot be installed by this package. Download it from
[services.healthtech.dtu.dk](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/),
follow its own installation instructions, then tell this tool where it is.

> **This package is developed and tested against netMHCpan 4.2.** It reads netMHCpan's
> output by column position, so another version whose output differs may not parse. It is
> not version-locked and other versions may well work — the checks below simply tell you
> quickly if one does not, rather than producing quietly wrong numbers.

```bash
export NETMHCPAN=/path/to/netMHCpan-4.2/netMHCpan
```

Put that `export` line in your `~/.bashrc` or `~/.zshrc` so you do not have to type it
again. You can also pass `--netmhcpan /path/to/netMHCpan` on any command, or set
`netmhcpan_path` in a config file. Check it worked:

```bash
$NETMHCPAN -h 2>&1 | head -3     # should print netMHCpan's own usage text
```

If you are on a version other than 4.2, run the test suite against it before trusting a
long run. Those tests check the assumptions this package makes about netMHCpan's output —
in particular how many columns `-BA` adds, and where:

```bash
uv run pytest -m netmhcpan       # needs $NETMHCPAN set, takes a few seconds
```

## Running the commands

The `netmhcpan-mp` command is installed inside the project's `.venv`, not system-wide, so
**`netmhcpan-mp` on its own will say `command not found`.** Use one of these:

```bash
# From the project folder. Nothing to remember, nothing to activate.
uv run netmhcpan-mp predict ...

# Or activate the environment first, then use the bare command for the rest of the session.
source .venv/bin/activate
netmhcpan-mp predict ...
deactivate                       # when you are done

# Or from anywhere at all, if you give the full path.
/path/to/netmhcpan_multiprocessor/.venv/bin/netmhcpan-mp predict ...
```

`uv run` is the one to use if you are not sure: it works from anywhere inside the project
folder, and it re-checks the environment each time, so a half-installed or stale `.venv`
repairs itself instead of failing confusingly.

The rest of this README writes `netmhcpan-mp`. Put `uv run` in front of it, or activate
first — whichever you picked above.

## Input and output

**Input** is a set of peptide parquet files, one per length, in a folder per length:

```
peptides/
    8mers/mydataset.parquet
    9mers/mydataset.parquet
   10mers/mydataset.parquet
   11mers/mydataset.parquet
```

Each file needs a `peptide_sequence` column (`--peptide-column` to use another name).
Everything in it is checked before any prediction starts: the column is there, no peptide
is null, and every peptide is the length its folder claims. netMHCpan only reports a
length mismatch as a short result, which is a long way to travel to find out the input was
wrong.

More than one parquet in a length folder is fine. Each is a separate dataset and keeps its
file name through to the output.

**Output** is one parquet per length per allele:

```
predictions/
    hla_a_02_01/
        8mers/mydataset.parquet
        9mers/mydataset.parquet
    hla_b_07_02/
        8mers/mydataset.parquet
        9mers/mydataset.parquet
```

| column | |
| --- | --- |
| `allele_slug` | the allele, matching the folder it lives in |
| `peptide_sequence` | the peptide, matching the input column |
| `score_el`, `rank_el` | eluted ligand likelihood and its percentile rank |
| `score_ba`, `rank_ba`, `affinity_nm` | binding affinity — only with `--scores ba` |
| `bind_level` | `SB`, `WB`, or null for a non-binder |
| `core`, `of`, `gp`, `gl`, `ip`, `il`, `icore` | the peptide's alignment to the 9mer binding core; all zero for 9mers, meaningful at every other length |

The rows come back in the order the peptides were given, so the predictions line up row for
row with the peptide parquet. That is checked rather than assumed: every chunk's returned
peptides are compared against the peptides it was sent, because a reordered result would
be internally consistent — the peptide column is netMHCpan's own output — and would break
a positional join silently.

Nothing else is carried through. Join back on `peptide_sequence`, or positionally, to pick
up `id`, `start`, `end` or anything else the peptide file holds:

```python
import polars as pl

peptides = pl.read_parquet('peptides/9mers/mydataset.parquet')
predictions = pl.read_parquet('predictions/hla_a_02_01/9mers/mydataset.parquet')
peptides.join(predictions, on='peptide_sequence')
```

## Alleles

Give an allele either as a number or as a slug, mixed freely in one run:

```bash
--allele 'HLA-A*02:01'    --allele hla_a_02_01    --allele H-2-Db
```

Every allele is resolved against netMHCpan's own list before anything runs, so a typo
fails in the first second rather than after a long prediction, and the error suggests what
you probably meant. The output always uses one canonical slug per allele, whichever
spelling was typed, so the same job never writes to two different folders. Pass
`--no-validate-alleles` to skip the check — resolution then falls back to string rules,
which cover HLA class I only.

## Scores: EL and BA

By default netMHCpan predicts eluted ligand likelihood, and that is what `--scores el`
gives you. `--scores ba` adds binding affinity — netMHCpan's `-BA` — which is *additional*
to the eluted ligand columns rather than a replacement for them, so a `ba` run's output is
a superset of an `el` run's.

```bash
netmhcpan-mp predict --scores ba ...
```

Affinity costs more time and three more columns, which is why `el` is the default.

## Modes and sources

Two axes describe where the peptides came from:

- `--mode` — the kind of peptide collection. `proteome` today; immunopeptidome and display
  screen modes are the reason this exists.
- `--source` — the concrete input handler registered for that mode. Optional: it defaults
  to the mode's own.

```
proteome  ->  {proteome}
```

Nothing downstream of the peptide list cares where the peptides came from — netMHCpan sees
a list of equal-length strings either way. What differs between collections is how the set
is laid out on disk and how its lengths are established, and that is all a mode and source
decide. Adding a kind of input means adding a source object in
`src/netmhcpan_multiprocessor/modes/`; the predictor does not change.

## Preparing a peptide set

`prepare` builds a peptide set that `predict` can then read. For the proteome mode, that
is every k-mer of every sequence in a FASTA:

```bash
netmhcpan-mp prepare \
    --input-path proteome.fasta \
    --output-folder peptides \
    --length 8 --length 9 --length 10 --length 11 \
    --dataset UP000008992
```

Peptides containing residues outside the twenty standard amino acids are dropped, since
netMHCpan cannot score them (`--keep-nonstandard` to keep them). The output is the
`<N>mers/<dataset>.parquet` layout above, ready to predict from.

This step needs the `proteome` extra (`polars-bio`, `shared-helpers`); prediction does not.
It is deliberately fenced off inside `modes/proteome/` and imported lazily, so it can move
into its own tool later without disturbing anything else.

## End to end, from a FASTA

The two commands are joined by nothing but the folder layout: `prepare`'s
`--output-folder` is `predict`'s `--input-folder`.

Run this from the project folder, with `uv run` spelled out so it works whether or not you
have activated anything:

```bash
export NETMHCPAN=/path/to/netMHCpan-4.2/netMHCpan

# 1. FASTA -> one peptide parquet per length
uv run netmhcpan-mp prepare \
    --input-path inputs/uniprotkb_proteome_UP000008992_2025_11_10.fasta \
    --output-folder outputs/peptides_parquet \
    --length 8 --length 9 --length 10 --length 11 \
    --dataset UP000008992_2025_11_10

# 2. those peptides -> one prediction parquet per length per allele
uv run netmhcpan-mp predict \
    --input-folder outputs/peptides_parquet \
    --output-folder outputs/netmhcpan_parquet \
    --allele 'HLA-A*01:01' --allele hla_a_02_01
```

which leaves:

```
outputs/
    peptides_parquet/
        8mers/UP000008992_2025_11_10.parquet         <- written by prepare
        9mers/UP000008992_2025_11_10.parquet
       10mers/UP000008992_2025_11_10.parquet
       11mers/UP000008992_2025_11_10.parquet
    netmhcpan_parquet/
        hla_a_01_01/
            8mers/UP000008992_2025_11_10.parquet     <- written by predict
            9mers/UP000008992_2025_11_10.parquet
           10mers/UP000008992_2025_11_10.parquet
           11mers/UP000008992_2025_11_10.parquet
        hla_a_02_01/
            ...
```

`--length` is only needed on `prepare`, which has to be told what to build. `predict`
takes every length it finds unless told otherwise, so the two cannot drift apart by
forgetting to update one of them.

For a real proteome, that is a few seconds for step 1 and a couple of minutes per allele
for step 2 on sixteen cores.

### Watching it more closely

`--verbose` (or `-v`) on either command prints the settings the run actually resolved to —
every field, including the ones that came from a config file, an alias, or a default —
before doing any work. On `predict` it also shows what each allele resolved to and which
input files were found, so a run against the wrong folder or the wrong allele can be
stopped in the first second rather than the fortieth minute:

```bash
netmhcpan-mp predict --input-folder outputs/peptides_parquet \
    --output-folder outputs/netmhcpan_parquet --allele hla_a_01_01 --verbose
```

```
Settings
┌───────────────────┬───────────────────────────────────────────────────┐
│ input_folder      │ outputs/peptides_parquet                          │
│ scores            │ el                                                │
│ netmhcpan_path    │ /home/chris/science/tools/netMHCpan-4.2/netMHCpan  │
│ chunk_size        │ 5000                                              │
│ workers           │ 16                                                │
│ ...               │                                                   │
└───────────────────┴───────────────────────────────────────────────────┘
Alleles
┏━━━━━━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ given       ┃ netMHCpan -a ┃ slug, and output folder ┃
┡━━━━━━━━━━━━━╇━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ hla_a_01_01 │ HLA-A01:01   │ hla_a_01_01             │
└─────────────┴──────────────┴─────────────────────────┘
Peptide files
┏━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ length ┃ dataset                ┃ path                                ┃
┡━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ 8mers  │ UP000008992_2025_11_10 │ .../peptides_parquet/8mers/...      │
│ 9mers  │ UP000008992_2025_11_10 │ .../peptides_parquet/9mers/...      │
└────────┴────────────────────────┴─────────────────────────────────────┘
```

Afterwards it prints anything netMHCpan wrote to stderr, which is otherwise thrown away on
a successful run — or says plainly that it wrote nothing. `--verbose` and `--quiet` ask for
opposite things and are refused together.

## Configuration

Every setting can come from a TOML or JSON file instead of the command line, and options
given on the command line win over the file. One file can configure both steps of a
pipeline: keys at the top level are shared, and a `[prepare]` or `[predict]` table
overrides them for that command.

```toml
# pipeline.toml
lengths = [8, 9, 10, 11]

[prepare]
input_path = "proteome.fasta"
output_folder = "peptides"
dataset = "UP000008992"

[predict]
input_folder = "peptides"
output_folder = "predictions"
alleles = ["HLA-A*02:01", "hla_b_07_02"]
scores = "el"
```

```bash
netmhcpan-mp prepare --config pipeline.toml
netmhcpan-mp predict --config pipeline.toml
```

### predict

| key / option | default | |
| --- | --- | --- |
| `input_folder` / `--input-folder` | — | the peptide set |
| `output_folder` / `--output-folder` | — | `<output>/<allele_slug>/<N>mers/<dataset>.parquet` |
| `alleles` / `--allele` | — | numbers and/or slugs, repeatable |
| `lengths` / `--length` | every length found | repeatable |
| `mode` / `--mode` | `proteome` | kind of peptide collection |
| `source` / `--source` | the mode's own | input handler |
| `scores` / `--scores` | `el` | `el`, or `ba` for affinity as well |
| `netmhcpan_path` / `--netmhcpan` | `$NETMHCPAN`, then `PATH` | the executable |
| `chunk_size` / `--chunk-size` | `5000` | peptides per netMHCpan process |
| `alleles_per_process` / `--alleles-per-process` | `1` | alleles scored in one process, or `0` for all |
| `workers` / `--workers` | CPU count | concurrent processes |
| `peptide_column` / `--peptide-column` | `peptide_sequence` | |
| `raw_output_folder` / `--raw-output-folder` | none | also keep netMHCpan's text output |
| `validate_alleles` / `--no-validate-alleles` | on | check against `-listMHC` |
| `overwrite` / `--no-overwrite` | on | off skips runs whose output exists |
| `strict` / `--strict` | off | stop on the first failed run |
| — / `--verbose`, `-v` | off | show resolved settings, alleles and inputs first, and netMHCpan's stderr after |
| — / `--quiet` | off | no progress bar and no summary table |

### prepare

| key / option | default | |
| --- | --- | --- |
| `input_path` / `--input-path` | — | the FASTA, for `mode=proteome` |
| `output_folder` / `--output-folder` | — | `<output>/<N>mers/<dataset>.parquet` |
| `lengths` / `--length` | `8 9 10 11` | repeatable |
| `mode` / `--mode` | `proteome` | |
| `dataset` / `--dataset` | the input file's stem | names the output files |
| `keep_nonstandard` / `--keep-nonstandard` | off | keep peptides netMHCpan cannot score |
| — / `--verbose`, `-v` | off | show the resolved settings before building |
| — / `--quiet` | off | no summary table |

## Using it from a pipeline

Both steps take a config object, a plain dictionary, or keyword arguments, and return
results rather than printing them. Nothing is written to the terminal unless you ask for it.

### prepare

```python
from netmhcpan_multiprocessor import prepare

prepared = prepare(
    input_path='inputs/uniprotkb_proteome_UP000008992_2025_11_10.fasta',
    output_folder='outputs/peptides_parquet',
    lengths=[8, 9, 10, 11],
    dataset='UP000008992_2025_11_10',
)

print(prepared.dataset, prepared.peptide_count, f"{prepared.elapsed_seconds:.1f}s")
for peptide_file in prepared.files:
    print(peptide_file.length, peptide_file.path, peptide_file.peptide_count)
```

`PrepareResult.files` is one `PreparedFile(length, path, peptide_count)` per length, in
ascending length order, and `peptide_count` is the total across them.

### prepare, then predict

The whole pipeline is two calls, joined by one folder — the same seam the two commands
use, so nothing has to be told twice:

```python
from netmhcpan_multiprocessor import predict, prepare

peptides = 'outputs/peptides_parquet'

prepared = prepare(
    input_path='inputs/uniprotkb_proteome_UP000008992_2025_11_10.fasta',
    output_folder=peptides,
    dataset='UP000008992_2025_11_10',
)

result = predict(
    input_folder=peptides,          # what prepare just wrote
    output_folder='outputs/netmhcpan_parquet',
    alleles=['HLA-A*01:01', 'hla_a_02_01'],
)

if result.failures:
    raise RuntimeError(f"{len(result.failures)} of {len(result.runs)} predictions failed")
```

`lengths` is left off `predict` deliberately: it takes every length it finds, which is
every length `prepare` wrote. Pass it only to predict a subset.

Preparation needs the `proteome` extra. Without it, `prepare()` raises `ConfigError` naming
the extra to install; `predict()` is unaffected either way.

One caveat on the "no terminal output" claim: `prepare()` is silent on stdout, but
polars-bio draws a progress bar on **stderr** while it reads the FASTA. tqdm decides
whether it is disabled when it is imported, so this package cannot switch it off after the
fact. Set it in the environment before the process starts if you want prepare fully quiet:

```bash
TQDM_DISABLE=1 python run_my_pipeline.py
```

`predict()` has no such dependency and is silent on both streams.

### Seeing what a run is actually doing

`RunResult.warnings` holds anything netMHCpan wrote to stderr while succeeding — usually
empty, and never fatal, but the first place to look when a run finishes and the numbers
look wrong. It is deduplicated, since every chunk tends to repeat the same complaint.

```python
for run in result.runs:
    for warning in run.warnings:
        logger.warning('%s: %s', run.label, warning)
```

On a failed run the same output goes into `run.error` instead, so it is never reported in
two places at once.

### predict

`predict()` takes a `PredictionConfig`, a plain dictionary, or keyword arguments, and
returns results rather than printing them.

```python
from netmhcpan_multiprocessor import predict

result = predict(
    input_folder='peptides',
    output_folder='predictions',
    alleles=['HLA-A*02:01'],
    lengths=[9, 10],
    scores='ba',
)

for run in result.runs:
    print(run.dataset, run.length, run.allele, run.peptide_count, run.strong_binders)

result.frame().write_parquet('run_summary.parquet')
```

A failed allele does not lose the rest of a long run: each run records its own outcome, and
`predict()` returns them all.

```python
if result.failures:
    for run in result.failures:
        print(run.label, run.error)
```

Pass `strict=True` to raise `PredictionFailedError` instead, cancelling whatever is still
queued. Runs abandoned that way are still reported, since a run is only assembled once its
last chunk lands and some may have got most of the way.

A failed run writes nothing, and **removes any output already sitting at its path**. That
file was written by an earlier run over a peptide set that may since have changed, so
leaving it would mean a later `--no-overwrite` retry skipped the path as though it held
current results. Successful writes go to a temporary file and are moved into place, so an
interrupted run leaves either the previous file or none — never half of one.

To show progress, pass anything with `add_run`, `advance` and `fail` — the CLI's rich
progress bar is one implementation of that protocol:

```python
class LoggingProgress:
    def add_run(self, label, total): return label
    def advance(self, handle, amount=1): pass
    def fail(self, handle, message): logger.error('%s: %s', handle, message)

predict(settings, progress=LoggingProgress())
```

Every error the package raises deliberately derives from `NetMHCpanMultiprocessorError`, so
one `except` clause catches the lot.

## When something goes wrong

Errors are meant to say what to do about them. These are the ones you are most likely to
meet, with what they actually mean.

**`netmhcpan-mp: command not found`**

The command lives inside the project's `.venv`, not system-wide. Run it as
`uv run netmhcpan-mp ...` from the project folder, or `source .venv/bin/activate` first.
See [Running the commands](#running-the-commands).

If you *have* activated the environment and the command is still not found, the `.venv` is
stale — most often because the project folder was renamed or moved after it was created,
which leaves `activate` pointing at a path that no longer exists. Everything inside it then
vanishes from `PATH` at once, which looks like nothing is installed. Rebuild it:

```bash
rm -rf .venv
uv sync --extra proteome
```

**`Error: Could not find netMHCpan. Set it with the --netmhcpan option, ...`**

netMHCpan is not installed, or this tool has not been told where it is. See
[step 3 of Installing](#3-netmhcpan). `echo $NETMHCPAN` shows whether the environment
variable survived into this terminal — it does not carry across sessions unless you put it
in your shell profile.

**`Error: No netMHCpan executable at /some/path`** — the path is set but wrong. Check for a
typo, and that it points at the `netMHCpan` script itself rather than the folder holding it.

**`Error: netMHCpan does not know the allele 'HLA-A*02:02'. Did you mean: HLA-A*02:01?`**

The allele is checked against netMHCpan's own list before anything runs, so this costs a
second rather than a wasted hour. Either fix the spelling or use one of the suggestions.
Numbers (`HLA-A*02:01`) and slugs (`hla_a_02_01`) are both accepted, so the form is not the
problem — the allele is.

**`Error: No peptide parquet files under peptides. The proteome source expects <length>mers folders`**

`--input-folder` should point at the folder that *contains* the `8mers/`, `9mers/` folders,
not at one of them, and not at a single parquet file. Run with `--verbose` to see exactly
where it looked and what it found.

**`Error: peptides/9mers/x.parquet is in a 9mers folder but holds peptides of length 8-9`**

The peptide file does not hold what its folder claims. This is checked up front because
netMHCpan reports a length mismatch only as a short result, which is a long way to travel
to discover the input was wrong.

**`Error: Preparing a proteome peptide set needs the optional FASTA dependencies`**

`prepare` needs the `proteome` extra, which `uv sync` alone does not install:

```bash
uv sync --extra proteome
```

**`Error: 1 of 4 predictions failed`**

The other three succeeded and their files were written; only the failed one is missing. The
summary table names which run failed and why. A failed run also removes any older output
sitting at its path, so nothing is left behind that looks current but is not.

**It is slower than expected**

Every netMHCpan process is single threaded, so throughput is `--workers` processes at once,
defaulting to your core count. `--verbose` shows what was resolved. If the machine is shared,
lower `--workers`; if a run is dominated by one long length, lower `--chunk-size` so the work
spreads more evenly.

For more than a handful of alleles, `--alleles-per-process 0` is the one setting that
changes the order of magnitude: the alleles share a process instead of repeating the same
allele-independent work. Because a batch multiplies what a chunk costs, pair it with a
smaller `--chunk-size` — aim for a chunk that takes seconds rather than minutes, or one
slow chunk will be the whole run's finishing time. The trade-off is isolation: a chunk
that fails takes every allele in its batch with it, so if one allele is doubtful, predict
it on its own.

**Still stuck**

`--verbose` prints the settings the run actually resolved to — including everything that
came from a config file or a default — plus what each allele resolved to and which input
files were found, all before any prediction starts.

## Tests

```bash
uv run pytest
```

The suite runs against a stand-in netMHCpan in `tests/fixtures/`, so it needs no netMHCpan
install: it writes the same banner, columns and prediction rows as the real program, and
can be told to fail in each of the ways netMHCpan actually fails — a non-zero exit, an
`Error.` printed to stdout followed by a zero exit, a silently short result, predictions
out of order, and enough output on stderr to deadlock a careless reader.

Tests against a real netMHCpan are marked `netmhcpan` and skipped unless `$NETMHCPAN`
points at one:

```bash
NETMHCPAN=/path/to/netMHCpan uv run pytest
```

They are worth running after a netMHCpan upgrade: they check the assumptions the stand-in
encodes, in particular how many columns `-BA` adds and where.
