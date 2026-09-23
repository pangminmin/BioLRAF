# BioLRAF GSE226365 upload steps

## 1. Add the method files

Upload the contents of `Methods/` to the `Methods/` directory of the BioLRAF
GitHub repository.

## 2. Add the example directories

Upload `Example/GSE226365/` without moving or renaming its subdirectories.

## 3. Prepare the gficf matrix locally

From the repository root, run:

```bash
Rscript Methods/BioLRAF_analysis.R \
  --reference local_data/tri_int.rda \
  --query local_data/GSE226365_anno.rda \
  --dataset GSE226365 \
  --output Example/GSE226365/intermediate/GSE226365_gficf_matrix.csv
```

## 4. Complete lineage scoring

Use the BioLRAF NMF and gene-set scoring workflow to convert the gficf matrix
into:

```text
Example/GSE226365/input/GSE226365_BioLRAF_scores.csv
```

The file must contain `cell_id`, `Os`, `Ch`, `Ad`, `MSC`, `celltype`, and
`dataset`.

## 5. Install Python dependencies

```bash
python -m pip install -r requirements.txt
```

## 6. Run the notebook

Start Jupyter from the BioLRAF repository root, open
`Methods/BioLRAF_visualization_GSE226365.ipynb`, and run all cells. The
notebook also supports being started directly from `Methods/`.

## 7. Commit the generated result files

Upload the HTML and PNG files from `Example/GSE226365/results/` to the same
directory in GitHub.

## 8. Link the results from the root README

Use the GitHub Pages URLs:

```text
https://pangminmin.github.io/BioLRAF/Example/GSE226365/results/GSE226365_3D_ternary_celltype.html
https://pangminmin.github.io/BioLRAF/Example/GSE226365/results/GSE226365_3D_ternary_lineage_dominant.html
```
