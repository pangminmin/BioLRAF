# GSE226365 BioLRAF example

This directory contains the input, intermediate-file description, and output
locations for the GSE226365 BioLRAF example.

## Required visualization input

Place the cell-level score table at:

```text
input/GSE226365_BioLRAF_scores.csv
```

The file must contain at least:

```text
cell_id, Os, Ch, Ad, MSC, celltype, dataset
```

## Visualization

Open `../../Methods/BioLRAF_visualization_GSE226365.ipynb` from the `Methods`
directory and run all cells. The notebook saves processed scores, interactive
HTML files, and four fixed PNG views for each coloring scheme to `results/`.

## Interactive results

- [Color by cell type](results/GSE226365_3D_ternary_celltype.html)
- [Color by dominant lineage state](results/GSE226365_3D_ternary_lineage_dominant.html)

