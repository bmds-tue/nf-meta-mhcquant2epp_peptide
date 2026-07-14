# nf-meta-mhcquant2epp_peptide

A small Nextflow adapter that turns [nf-core/mhcquant](https://github.com/nf-core/mhcquant)'s
own samplesheet plus its per-`Sample_Condition` peptide identification output
(`<Sample>_<Condition>.tsv`) into the `sample,alleles,mhc_class,filename` samplesheet
[nf-core/epitopeprediction](https://github.com/nf-core/epitopeprediction) requires as `--input`.
Implemented entirely with native Nextflow channel operators (`splitCsv`, `map`, `unique`,
`groupTuple`, `join`, `collectFile`) -- no separate join-key column or manual bookkeeping
needed.

The logic lives in a named workflow, `mhcquant2epp`, which reads its inputs from `params`
and writes `samplesheet.csv` (and `conflicts.tsv`, if any) under `--outdir`.

## The problem this solves

`mhcquant` and `epitopeprediction` have genuinely incompatible samplesheet conventions:
`mhcquant` requires exactly `ID, Sample, Condition, ReplicateFileName` (capitalized),
while `epitopeprediction` requires exactly `sample, alleles, mhc_class, filename`
(lowercase) -- and `mhcquant`'s samplesheet is one row per raw MS replicate file, while
its own output is one merged file per `Sample_Condition`. This adapter derives the
`Sample_Condition` join key internally from `Sample`/`Condition` (the same construction
mhcquant's own output filenames use), collapses replicate-level duplicate rows down to
one row per `Sample_Condition`, and reports (rather than silently drops or duplicates)
any rows that disagree on `alleles`/`mhc_class` under the same key.

## Usage

```bash
nextflow run main.nf \
  --samplesheet mhcquant_samplesheet.csv \
  --mhcquant_outdir results/mhcquant \
  --outdir results/epitopeprediction_input
```

`--samplesheet` can be the exact same file passed to mhcquant's own `--input`, as long
as `alleles` and `mhc_class` columns have been added -- mhcquant's schema has no
`additionalProperties: false`, so the extra columns are safely ignored by mhcquant.

## Params

| Param | Default | Description |
|---|---|---|
| `samplesheet` | (required) | Path to the samplesheet mhcquant's `--input` consumes (or a superset of it). One row per raw MS replicate file. Always parsed as comma-separated. |
| `mhcquant_outdir` | (required) | Directory containing mhcquant's per-`Sample_Condition` peptide TSVs, e.g. `HLE_p1011.tsv`. |
| `outdir` | (required) | Directory `samplesheet.csv` (and `conflicts.tsv`, if any) are written to. |
| `sample_col` | `Sample` | Column holding the cell line / sample name. |
| `condition_col` | `Condition` | Column holding the plasmid/condition name. |
| `alleles_col` | `alleles` | Column holding the semicolon-separated HLA allele list, e.g. `A*02:01;A*31:01;B*18:01;B*40:01`. |
| `mhc_class_col` | `mhc_class` | Column holding the MHC class (`I` or `II`). |

`sample_col`/`condition_col` default to mhcquant's own required casing, and
`alleles_col`/`mhc_class_col` default to lowercase (they aren't part of mhcquant's
schema, so there's no casing constraint forcing a particular choice) -- so the same
samplesheet mhcquant consumes can be handed to this adapter directly.

`"${row[sample_col]}_${row[condition_col]}"` is the derived join key, using the same
string construction mhcquant's own output filenames are parsed with, so both sides of
the join are guaranteed to agree without a separately-authored key column.

## Output

- **`samplesheet.csv`** -- `sample,alleles,mhc_class,filename`, one row per unique
  `Sample_Condition` that has both a samplesheet entry and a matching mhcquant output
  TSV. Rows are sorted by `sample`. A `Sample_Condition` present in only one of the two
  inputs (samplesheet or `mhcquant_outdir`) is silently excluded, matching mhcquant's own
  filtering semantics.
- **`conflicts.tsv`** -- written only if two or more rows share a derived
  `Sample_Condition` key but disagree on `alleles`/`mhc_class` (a real data-entry error,
  not the expected replicate-row duplication). Absent when there are no conflicts. The
  main `samplesheet.csv` still gets exactly one row for a conflicting key (the
  first-occurrence value) -- the conflict is surfaced for follow-up, not silently
  resolved by dropping the row.

## Params validation & help

Parameters are described in `nextflow_schema.json` (types, defaults, required
constraints) via the [nf-schema](https://nextflow-io.github.io/nf-schema/) plugin. Every
run validates params against that schema before anything executes, so a missing required
param or a nonexistent input path fails fast with a clear message instead of an obscure
downstream error.

```bash
nextflow run main.nf --help        # all params, grouped
```

## Test

```bash
nf-test test test/main.nf.test
```
