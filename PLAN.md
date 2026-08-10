# PLAN — `netmhcpan_multiprocessor` as a library

> **Status: done.** Built as planned, with three additions agreed during the work and
> recorded in [§14](#14-what-changed-during-the-build): a `--scores el|ba` switch, the
> `--mode`/`--source` split, and preparation absorbed into the package. Verified against
> the real netMHCpan: output is byte-identical to the old script's, in input order.


Turn the single-file `run_netmhcpan.py` script into an installable Python package with a
stable programmatic entry point and a thin CLI over it, so the same prediction run can be
driven either from a terminal or from a larger pipeline that already holds its settings in
a dictionary.

Nothing about how the prediction is done changes: chunked peptide lists, one netMHCpan
process per chunk, a flat longest-first thread pool, per-run failure isolation. What
changes is where the settings come from, what the inputs and outputs are, and where the
boundaries between the pieces sit — in particular, a boundary between *how a peptide set
arrives* and *how it is predicted*, so new kinds of peptide collection can be added later
without touching the predictor.

---

## 1. Decisions taken

Settled with the user before writing this plan:

| Question | Decision |
| --- | --- |
| Output layout | One parquet **per length per allele**: `<output>/<allele_slug>/<length>mers/<dataset>.parquet`. Unchanged from today. |
| Input layout | Mirrors what `prepare_proteome.py` writes: `<input>/<length>mers/<dataset>.parquet`. |
| Raw netMHCpan text | Optional, **off by default** (`--raw-output-folder` / `raw_output_folder`). |
| Passthrough columns | Predictions only (`allele_slug`, `peptide_sequence`, prediction columns), in the input parquet's row order so a positional or `peptide_sequence` join back is trivial. |
| netMHCpan location | A config / CLI option, never hard-coded. |
| Vocabulary | "proteome" is not the general unit. A proteome is one kind of peptide collection among several — immunopeptidomics, display screens — so the neutral terms throughout are **peptide set** and **dataset**. |
| `--mode` / `--source` | Two axes. `--mode` names the kind of peptide collection; `--source` names the concrete input handler registered **for that mode**, defaulting to the mode's own. Registry keyed by `(mode, source)`. |
| Proteome preparation | Absorbed into the package as `netmhcpan-mp prepare --mode proteome`, fenced inside the proteome mode folder so it can be lifted out later by deleting one directory. |

Taken as author's judgement, flagged here so they can be overruled:

- **Package name** `netmhcpan_multiprocessor`, laid out under `src/`, console script
  `netmhcpan-mp` with two subcommands, `predict` and `prepare`.
- **`rich` stays out of the library core.** Progress is reported through a small protocol;
  the rich display is a CLI-side implementation of it. A pipeline importing the package
  gets no terminal output unless it asks for it.
- **Config files** (`--config x.toml|.json`) are supported by the CLI on top of the config
  dictionary, since the dictionary already has to exist for the Python API. Explicit CLI
  options override file values.
- **`prepare` stops writing the one-per-line `.txt` peptide files.** Parquet is the stated
  input contract; the chunk files `-p` needs are written to a temporary directory at
  predict time and never persisted.

---

## 2. Package layout

```
src/netmhcpan_multiprocessor/
    __init__.py      re-exports: predict, prepare, PredictionConfig, PrepareConfig, errors
    config.py        PredictionConfig / PrepareConfig + from_dict/from_file + validation
    alleles.py       slugify, allele-number <-> slug resolution, -listMHC index
    parsing.py       row schema, is_prediction_row, build_predictions_frame
    runner.py        the netMHCpan subprocess call for one chunk
    chunking.py      peptide list -> temporary chunk files
    predict.py       predict(): planning, thread pool, assembly, writing
    progress.py      ProgressReporter protocol + NullProgressReporter
    errors.py        error hierarchy
    modes/
        __init__.py  the (mode, source) registry: get_mode, resolve_source
        base.py      Mode / PeptideSource protocols, PeptideFile dataclass
        proteome/
            __init__.py   ProteomeMode: name, sources, default_source, preparer
            source.py     reads <input>/<N>mers/*.parquet
            prepare.py    FASTA -> <output>/<N>mers/<dataset>.parquet
    cli.py           click group, rich progress + summary table
tests/
    ...
    fixtures/fake_netmhcpan.py    stand-in binary used by the tests
```

`run_netmhcpan.py` and `prepare_proteome.py` are both deleted at the end — their contents
now live in the package.

---

## 3. Modes and sources

Nothing downstream of the peptide list cares where the peptides came from: netMHCpan sees
a list of equal-length strings either way. What differs between collections is **how the
peptide set is laid out on disk, how its lengths are established, and what has to happen
before it exists at all**. That is the whole job of this layer, and it is the only part
that a new kind of input touches.

```python
@dataclass(frozen=True)
class PeptideFile:
    length: int
    path: Path
    dataset: str        # file name without suffix; names the output file

class PeptideSource(Protocol):
    name: str
    def discover(self, input_folder: Path, lengths: tuple[int, ...] | None) -> list[PeptideFile]: ...
    def read_peptides(self, peptide_file: PeptideFile, peptide_column: str) -> list[str]: ...

class Mode(Protocol):
    name: str
    sources: dict[str, PeptideSource]
    default_source: str
    preparer: Preparer | None      # None = this mode has no preparation step
```

The registry is keyed by `(mode, source)`:

```
proteome         -> {proteome}            <- built now, default source 'proteome'
immunopeptidome  -> {immunopeptidome}     <- later
display_screen   -> {display_screen}      <- later
```

- `--source` is optional; omitted, it takes the mode's `default_source`, so
  `--mode proteome` alone is enough today.
- `get_mode(name)` raises `UnknownModeError` listing registered modes.
  `resolve_source(mode, name)` raises `UnknownSourceError` listing **that mode's** sources,
  so asking for a source that exists under a different mode gives an accurate message
  rather than a bare "not found".
- The intended growth path is that a source object can be registered under several modes.
  A future pre-processed-parquet input becomes one `PeptideParquetSource` registered into
  every mode, giving `proteome -> {proteome, peptide_parquet}` without a second copy of
  anything. Nothing in the design assumes one source per mode.

**`proteome`**, the source built now: the k-merised layout `<input>/<N>mers/*.parquet`.
Length comes from the folder name, and every peptide in the file is checked against it —
netMHCpan only reports a length mismatch as a short result, so it is worth catching before
a long run rather than after one.

## 4. Preparation

`prepare` is a **separate command and a separate function** from `predict`. It produces
peptide parquet; predict consumes peptide parquet; neither knows about the other. That
keeps the stated contract — predict's input is always a set of peptide parquet files —
while letting the proteome k-merisation live here for now.

`modes/proteome/prepare.py` is `prepare_proteome.py` reworked into a function:

```python
def prepare(config: PrepareConfig | dict) -> PrepareResult
```

FASTA in, `<output_folder>/<N>mers/<dataset>.parquet` out, one file per requested length,
carrying `peptide_sequence`, `id`, `start`, `end` as today, non-standard amino acids
filtered out, streamed with `sink_parquet`. `dataset` defaults to the FASTA file's stem and
is overridable, since the current file name (`uniprotkb_proteome_UP000008992_2025_11_10`)
is not the dataset name anyone wants on the output.

**Fenced for removal**: `polars-bio` and `shared-helpers` are needed only here, so they
move to an optional `proteome` extra and are imported inside the function, with a clear
message naming the extra if they are missing. `predict` imports nothing from this folder.
Lifting proteome preparation out later means deleting `modes/proteome/prepare.py`, dropping
the extra, and setting `preparer = None` — and a mode with `preparer is None` already has
to be handled, because the future modes will have one before they have the other.

## 5. Configuration

Two frozen dataclasses, one per command; each is the single source of truth for its
command, and the CLI builds one and calls the same function a pipeline would.

```python
@dataclass(frozen=True)
class PredictionConfig:
    input_folder: Path                  # <input>/<length>mers/<dataset>.parquet
    output_folder: Path                 # <output>/<slug>/<length>mers/<dataset>.parquet
    alleles: tuple[str, ...]            # numbers and/or slugs, mixed freely
    lengths: tuple[int, ...] | None = None   # None -> discover from input folder
    mode: str = 'proteome'
    source: str | None = None           # None -> the mode's default source
    netmhcpan_path: Path = ...          # option -> $NETMHCPAN -> which('netMHCpan')
    chunk_size: int = 5000
    workers: int = os.cpu_count()
    peptide_column: str = 'peptide_sequence'
    raw_output_folder: Path | None = None    # None -> no .txt written
    validate_alleles: bool = True       # check against netMHCpan -listMHC
    overwrite: bool = True              # False -> skip runs whose parquet exists
    strict: bool = False                # True -> raise on first failed run

@dataclass(frozen=True)
class PrepareConfig:
    input_path: Path                    # the FASTA
    output_folder: Path
    lengths: tuple[int, ...] = (8, 9, 10, 11)
    mode: str = 'proteome'
    dataset: str | None = None          # None -> input file stem
    keep_nonstandard: bool = False      # True -> keep peptides with non-standard residues
```

- `from_dict(mapping)` — unknown keys raise, paths coerced from `str`, `allele`/`length`
  accepted as singular aliases for ergonomics.
- `from_file(path)` — TOML (stdlib `tomllib`) or JSON. Keys may sit at the top level or
  under a `[predict]` / `[prepare]` table; the command's own table wins, so one file can
  configure both steps of a pipeline.
- Validation at construction: folders exist, lengths positive, `chunk_size`/`workers` >= 1,
  at least one allele, mode and source resolve.
- **netMHCpan resolution order**: explicit `netmhcpan_path` -> `NETMHCPAN` environment
  variable -> `shutil.which('netMHCpan')`. If none resolves to an executable file, raise
  `NetMHCpanNotFoundError` naming all three routes.

## 6. Alleles: numbers and slugs

Three forms are in play, and each has one job:

| Form | Example | Used for |
| --- | --- | --- |
| Allele number | `HLA-A*02:01` | what a human writes; what the summary shows |
| netMHCpan name | `HLA-A02:01` | what is passed to `-a` |
| Allele slug | `hla_a_02_01` | folder names, the `allele_slug` column |

`alleles.resolve_alleles(values, netmhcpan_path, validate)` maps every supplied value —
number *or* slug, mixed within one run — onto a `ResolvedAllele(number, netmhcpan_name,
slug)`.

- The authoritative index is built once per run from `netMHCpan -listMHC` (~13k entries
  here, one cheap subprocess, cached per path). Each listed name is indexed under its own
  slug **and** under the slug of its starred reconstruction, so both `hla_a_02_01` and the
  list's own `hla_a02_01` resolve, along with non-HLA names like `H-2-Db` and
  `Mamu-A1*001:01`.
- A pure-string fallback (regex on `hla_<gene>_<field>_<field>`) covers the case where the
  index is unavailable or `validate_alleles=False`, so the library is usable — and
  testable — without the real binary.
- An unresolvable allele raises `UnknownAlleleError` listing near matches. This catches
  typos before a 40-minute run rather than after it, which the current script does not.

## 7. Output

`<output_folder>/<allele_slug>/<length>mers/<dataset>.parquet`, written from
`build_predictions_frame` exactly as today: chunks reassembled by index so the frame keeps
the input parquet's row order, allowing a positional join back to the peptide file.
Optional raw text goes to `<raw_output_folder>/<dataset>__<N>mers_<slug>.txt`, with
positions renumbered across chunk boundaries as now.

## 8. The programmatic entry point

```python
from netmhcpan_multiprocessor import predict

result = predict({
    'input_folder': 'outputs/peptides_parquet',
    'output_folder': 'outputs/netmhcpan_parquet',
    'alleles': ['HLA-A*02:01', 'hla_a_01_01'],
    'lengths': [8, 9, 10, 11],
    'mode': 'proteome',
    'netmhcpan_path': '/home/chris/science/tools/netMHCpan-4.2/netMHCpan',
})
```

`predict(config, *, progress=None)` accepts a `PredictionConfig`, a plain dict, or keyword
arguments, and returns a `PredictionResult`:

```python
@dataclass
class RunResult:
    length: int; allele: str; allele_slug: str; dataset: str
    parquet_path: Path | None; raw_path: Path | None
    peptide_count: int; chunk_count: int
    strong_binders: int; weak_binders: int
    error: str | None
    @property
    def succeeded(self) -> bool: ...

@dataclass
class PredictionResult:
    runs: list[RunResult]; elapsed_seconds: float
    @property
    def failures(self) -> list[RunResult]: ...
    def frame(self) -> pl.DataFrame:  # the summary as data, for a pipeline to log
```

A caller therefore gets machine-readable outcomes rather than a printed table, and decides
itself whether a partial failure is fatal (or sets `strict=True` to have it raised).
`prepare(config)` likewise returns a `PrepareResult` of per-length output paths and counts.

## 9. CLI

`netmhcpan-mp` (also `python -m netmhcpan_multiprocessor`), a click group of two commands
whose options map 1:1 onto config keys.

```
netmhcpan-mp predict
  --input-folder / --output-folder     (required unless given in --config)
  --allele            repeatable; number or slug
  --length            repeatable; default all lengths found
  --mode              kind of peptide set [default: proteome]
  --source            input handler for that mode [default: the mode's own]
  --netmhcpan         path to the netMHCpan executable [env: NETMHCPAN]
  --chunk-size / --workers
  --raw-output-folder
  --config            TOML or JSON file of the same keys
  --no-validate-alleles / --no-overwrite / --strict
  --quiet

netmhcpan-mp prepare
  --mode              [default: proteome]; errors if the mode has no preparer
  --input-path        the FASTA, for mode=proteome
  --output-folder     written as <output>/<N>mers/<dataset>.parquet
  --length            repeatable [default: 8 9 10 11]
  --dataset           [default: input file stem]
  --keep-nonstandard
  --config / --quiet
```

`--source`'s valid values depend on `--mode`, which click cannot express statically, so it
is validated after parsing and the error lists that mode's sources; `--help` shows the
mode-to-source table. Behavioural parity with today otherwise: rich progress bar per run,
the summary table, and a non-zero exit when any run fails.

The `--proteome IDENTIFIER` option disappears. The dataset name now comes from input file
names, which is what lets the tool run on any peptide set rather than only on a named
proteome.

## 10. Tests

`pytest`, added as a dev dependency group. The point of making the binary path a config
option is that the pipeline becomes testable without netMHCpan installed:
`tests/fixtures/fake_netmhcpan.py` is an executable stub that parses `-p`/`-a`, emits a
realistic banner, header, prediction rows and summary, and can be told through environment
variables to fail, to print `Error.` and exit 0, to truncate its output, or to sleep.

Unit:
- `slugify`, number/slug/netMHCpan-name conversions, round trips.
- Allele resolution from both forms; unknown allele raises with suggestions; fallback path
  with validation off; `-listMHC` index built once and cached.
- Registry: `proteome` mode and source resolve; `--source` defaults to the mode's own;
  unknown mode lists modes; unknown source lists *that mode's* sources; a mode with no
  preparer is reported clearly by `prepare`.
- `PredictionConfig` / `PrepareConfig` `from_dict` / `from_file`: coercion, defaults,
  unknown keys, bad values, `[predict]` / `[prepare]` tables, CLI-overrides-file, and the
  netMHCpan resolution order (option > env > PATH) with its error.
- `is_prediction_row` on banner, header, prediction, `pos >= 1000`, summary lines.
- `build_predictions_frame`: schema and dtypes, `bind_level` null / `SB` / `WB`, wrong
  field count raises, row order preserved.
- Proteome source discovery: length folders, explicit vs discovered lengths, missing length
  errors, wrong peptide length, missing column, several datasets in one length folder.
- Chunking: exact multiples, remainder, single chunk, ordering.
- `prepare` on a small FASTA fixture: one parquet per length, columns and coordinates,
  non-standard residues filtered (and kept under `--keep-nonstandard`), sequences shorter
  than the length contributing nothing, dataset naming.

Integration, against the fake binary:
- End-to-end `predict()` over a small two-length two-allele fixture: files land at the
  expected paths, row count matches the input, row **order** matches the input, the
  `allele_slug` column matches its folder.
- `prepare` then `predict` over the same temporary folder: the layout one writes is the
  layout the other discovers. This is the test that keeps the two halves honest.
- Chunk reassembly: `chunk_size` smaller than the input, deliberately out-of-order
  completion, output still in input order.
- Failure isolation: one allele fails, the others still write their parquet, the result
  records the error, the CLI exits non-zero, `strict=True` raises instead.
- `Error.` + exit 0 is caught; a short result is caught.
- Raw text written only when `raw_output_folder` is set.
- `overwrite=False` skips existing outputs.
- CLI via `click.testing.CliRunner`: options and `--config` produce identical runs; mixed
  number/slug alleles; `--mode`/`--source` validation errors.

Real-binary tests are marked `@pytest.mark.netmhcpan` and skipped unless `NETMHCPAN` is
set, so a full check stays possible on this machine without the suite depending on it.

## 11. Packaging

`pyproject.toml`: `src/` layout, `hatchling` backend, `[project.scripts] netmhcpan-mp`,
runtime deps `click`/`polars`/`rich`, an optional `proteome` extra carrying `polars-bio`
and `shared-helpers` for the preparation step, and a `dev` group with `pytest`.
`.gitignore` covers `inputs/`, `outputs/`, `tmp/`, `.venv/`, `__pycache__/`, `*.egg-info`,
`.pytest_cache/` — no data of any kind goes to the repo, only code.

## 12. Order of work

1. Package skeleton, `pyproject.toml`, `.gitignore`.
2. `errors.py`, `alleles.py`, `parsing.py` (lifted from the script), `chunking.py`.
3. `modes/base.py` + registry, `modes/proteome/source.py`, `config.py`.
4. `modes/proteome/prepare.py` (lifted from `prepare_proteome.py`).
5. `progress.py`, `runner.py`, `predict.py`.
6. `cli.py`: the group, both commands, the rich reporter.
7. Fake netMHCpan fixture, then the test suite; run it.
8. One real-binary smoke run against `outputs/peptides_parquet`, confirming the new output
   matches the existing `outputs/netmhcpan_parquet` files row for row.
9. Delete `run_netmhcpan.py` and `prepare_proteome.py`.
10. `README.md`: install, both commands, Python usage, config reference, the mode/source
    model and how to add one, output schema, the netMHCpan prerequisite, testing.

## 13. Risks

- **Slug -> allele number is not universally invertible.** Handled by resolving against
  netMHCpan's own allele list rather than by string rules alone; the rule-based fallback is
  documented as covering HLA class I only.
- **Two axes, one pair.** `(mode, source)` carries exactly one entry today, so the registry
  is deliberately thin — a protocol, a dict, and a resolver, no plugin machinery, no entry
  points. The cost if a second mode never appears is one small module; the benefit if the
  input type does turn out to be pre-processed parquet is that it is a new source object
  and nothing else.
- **Preparation is in the package but not of it.** Fenced behind an optional extra and a
  lazy import, with `preparer = None` a supported state, so its eventual removal is a
  deletion rather than an unpick. Worth re-checking at review time that nothing outside
  `modes/proteome/` has grown an import of it.
- **Removing `--proteome` changes the CLI**, as does moving `prepare_proteome.py` behind a
  subcommand. Both are deliberate and are called out in the README.
- **Step 8 needs the real binary** (~1–2 min per allele-length here). If it is unavailable
  when the work runs, it is reported as unverified rather than skipped silently.

## 14. What changed during the build

Everything above was built as planned. Three things were added or settled after the plan
was approved, and one class of defect was found by review afterwards.

**Added on request, mid-build**

- **`--scores el|ba`**, defaulting to `el`. netMHCpan's `-BA` adds binding affinity to the
  eluted ligand output rather than replacing it, so `ba` is the wider output, not the other
  one. Confirmed against the real program: `-BA` adds three columns — `Score_BA`,
  `%Rank_BA` and `Aff(nM)` — between `%Rank_EL` and `BindLevel`. This matters more than it
  sounds: an affinity non-binder row is 16 fields where an eluted-ligand *binder* row is
  15, so rows are told apart by the `<=` marker rather than by width, and reading one
  score set's output as the other's is refused rather than silently misparsed.
- **`--mode` and `--source` as two axes**, registry keyed by `(mode, source)`.
- **Proteome preparation absorbed** into `modes/proteome/`, behind the `proteome` extra and
  a lazy import, as `netmhcpan-mp prepare`.
- **`--verbose` / `-v`** on both commands, for developing a pipeline: the resolved settings,
  what each allele resolved to and which input files were found, printed *before* the run
  so a wrong folder or allele costs a second rather than an hour. It also surfaces what
  netMHCpan wrote to stderr, which was previously discarded on a successful run — that is
  now kept on `RunResult.warnings`, deduplicated across chunks and capped, so a pipeline
  can log it without the flag.

**Found by review, after the tests passed**

The suite was green and the output matched the old script's byte for byte before any of
these were found, which is the argument for the review pass rather than against the tests.

- **A stderr deadlock.** stdout and stderr were both pipes, and stdout was drained to EOF
  before stderr was read. A netMHCpan that wrote more than a pipe buffer of warnings — easy
  over a 5000 peptide chunk — would block on its stderr write, stop producing stdout, and
  hang the read for ever, taking the whole run with it. Confirmed by reproduction. stderr is
  now drained by its own thread.
- **Killing a shell-wrapped netMHCpan.** netMHCpan is a wrapper script, so signalling the
  process left the real binary holding the pipes. The child now gets its own process group
  and the group is killed. The kill also had to move ahead of the pipe waits, since an
  `except` that runs after a `finally` that waits is no kill at all.
- **Failure isolation had holes.** Only this package's own exceptions were caught per chunk,
  so anything else — a caller's progress reporter raising, a polars cast failure while
  building the frame — escaped and lost *every* run's results, in exactly the case the
  contract promises not to. Progress reporting is now guarded outright, since a display
  should never be able to fail a prediction.
- **A failed run left the previous run's parquet in place**, which a later `--no-overwrite`
  retry would then skip as though it held current results — silently joining predictions for
  a peptide set that may have changed. Failed runs now clear their output, and successful
  ones are written to a temporary file and moved into place.
- **An alias regression I introduced.** Rejecting `allele` beside `alleles` in one mapping
  also broke `predict(config, allele=...)`, because the merge sites build a mapping holding
  every canonical field name. Aliases are now resolved at each mapping's own boundary,
  before anything is merged.

**Added because the review asked the right question**

Nothing verified that netMHCpan answered in the order it was asked. A reordered result would
have been internally consistent — the peptide column is netMHCpan's own output — and would
have broken the documented positional join with nothing to show for it. Every chunk's
returned peptides are now checked against the peptides it was sent.
