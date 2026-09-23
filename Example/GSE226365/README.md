# GSE226365 BioLRAF example

This directory provides a reproducible BioLRAF analysis example using
the GSE226365 single-cell dataset.

## Directory structure

- `input/`: input cell-level gficf lineage scores.
- `processed/`: processed BioLRAF scores and lineage-state metrics.
- `results/`: interactive HTML and fixed-view PNG visualizations.

## Input

The input file is:

`input/GSE226365_gficf_scores.csv`

It must contain:

- `cell_id`
- `Os`
- `Ch`
- `Ad`
- `MSC`
- `celltype`
- `dataset`

## Processed output

The BioLRAF workflow generates:

`processed/GSE226365_BioLRAF_scores.csv`

This table contains normalized lineage activities, `MSC_ratio`,
`Diff_score`, and `lineage_dominant`.

## Run the example

From the repository root:

```bash
python -m pip install -r requirements.txt
jupyter notebook
