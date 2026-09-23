# Processed BioLRAF scores

`GSE226365_BioLRAF_scores.csv` is generated from
`../input/GSE226365_gficf_scores.csv` using the BioLRAF scoring workflow.

The processed table contains the original lineage-associated scores
(`Os`, `Ch`, `Ad`, and `MSC`), cell metadata, normalized lineage activities
(`Os_n`, `Ch_n`, and `Ad_n`), `MSC_ratio`, `Diff_score`, and
`lineage_dominant`.

This processed file is used as the common input for the 3D ternary
visualizations colored by `celltype` and `lineage_dominant`.
