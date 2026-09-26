# GSE226365 BioLRAF example

This directory provides a reproducible BioLRAF analysis example using the GSE226365 single-cell dataset.

## Directory structure

- `input/`: input files used for analysis.
- `processed/`: processed BioLRAF score tables.
- `results/`: interactive HTML plots and fixed-view PNG visualizations.

## Input files

- `GSE226365_data.rds`: annotated GSE226365 query Seurat object.
- `Reference_data.rds`: BioLRAF reference Seurat object.
- `gene_set.rds`: gene-set list used for lineage-associated scoring.
- `GSE226365_gficf_scores.csv`: cell-level gficf lineage-associated scores used by the visualization notebook. This file can be regenerated from the three RDS files using the analysis script.

## Workflow

### 1. Generate gficf scores in R

Run this command from the repository root:

```bash
Rscript Methods/BioLRAF_analysis.R \
  --reference Example/GSE226365/input/Reference_data.rds \
  --query Example/GSE226365/input/GSE226365_data.rds \
  --gene-list Example/GSE226365/input/gene_set.rds \
  --dataset GSE226365 \
  --output Example/GSE226365/input/GSE226365_gficf_scores.csv
