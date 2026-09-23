# BioLRAF

**BioLRAF (Biomaterial-induced Lineage Remodeling Assessment Framework)** is a reference-guided single-cell framework for visualizing and quantitatively assessing mesenchymal stem cell (MSC) lineage remodeling induced by biomaterials.

Rather than evaluating differentiation using only selected lineage markers or population-averaged measurements, BioLRAF represents each cell within a shared lineage-state space spanning osteogenic, chondrogenic, and adipogenic states. It provides complementary measurements of lineage-associated activity, differentiation status, and dominant lineage-state composition, thereby enabling comparisons of heterogeneous cellular responses across datasets, experimental conditions, and biomaterial systems.

## BioLRAF workflow

The BioLRAF workflow consists of the following main steps:

1. **Preparation of the single-cell expression matrix**

   A cell-by-gene expression matrix is generated from a processed single-cell RNA-sequencing dataset. Cell-level metadata, including dataset identity, experimental condition, and cell type, are retained for downstream analysis.

2. **gficf transformation and lineage-associated activity estimation**

   The expression matrix is transformed using the gene-frequency–inverse-cell-frequency (gficf) approach. Reference-derived gene sets are then used to estimate osteogenic (`Os`), chondrogenic (`Ch`), adipogenic (`Ad`), and undifferentiated MSC (`MSC`) activities for each individual cell.

3. **Normalization of lineage-associated activities**

   The three differentiation-associated activities are normalized as:

   $$
   Os_n = \frac{Os}{Os + Ch + Ad}
   $$

   $$
   Ch_n = \frac{Ch}{Os + Ch + Ad}
   $$

   $$
   Ad_n = \frac{Ad}{Os + Ch + Ad}
   $$

   These normalized values describe the relative position of each cell within the osteogenic–chondrogenic–adipogenic lineage-state space.

4. **Quantification of differentiation status**

   The relative undifferentiated MSC activity is calculated as:

   $$
   MSC\_ratio = \frac{MSC}{Os + Ch + Ad + MSC}
   $$

   The differentiation score is subsequently defined as:

   $$
   Diff\_score = 1 - MSC\_ratio
   $$

   A higher `Diff_score` indicates a stronger departure from the undifferentiated MSC-associated state.

5. **Classification of dominant lineage states**

   Cells located near the center of the ternary lineage-state space are classified as `Non-dominant`. For the remaining cells, the largest normalized lineage-associated activity determines whether a cell is classified as:

   * `Os-dominant`
   * `Ch-dominant`
   * `Ad-dominant`
   * `Non-dominant`

6. **Three-dimensional lineage-state visualization**

   Each cell is projected into a three-dimensional ternary space. Its horizontal position represents the relative osteogenic, chondrogenic, and adipogenic activities, while its vertical position represents `Diff_score`.

   Cells can be colored by metadata such as `celltype`, experimental group, dataset, or `lineage_dominant`. The resulting visualization therefore displays lineage direction, differentiation status, and population heterogeneity within a unified coordinate system.

## BioLRAF outputs

BioLRAF generates a cell-level processed table containing the original lineage-associated activities, normalized lineage coordinates, differentiation metrics, dominant lineage classification, and available cell metadata.

The principal output fields include:

| Output                 | Description                                                  |
| ---------------------- | ------------------------------------------------------------ |
| `Os`, `Ch`, `Ad`       | Original osteogenic, chondrogenic, and adipogenic activities |
| `MSC`                  | Undifferentiated MSC-associated activity                     |
| `Os_n`, `Ch_n`, `Ad_n` | Normalized lineage-associated activities                     |
| `MSC_ratio`            | Relative undifferentiated MSC-associated activity            |
| `Diff_score`           | Differentiation status calculated as `1 − MSC_ratio`         |
| `lineage_dominant`     | Dominant lineage-state classification                        |
| `celltype`             | Cell-type annotation                                         |
| `dataset`              | Dataset or sample identity                                   |

The visualization workflow produces:

* Interactive three-dimensional HTML plots;
* Fixed-view PNG images for figure preparation;
* Cell-type-colored lineage-state maps;
* Dominant-lineage-colored lineage-state maps;
* Cell-level BioLRAF score tables for downstream statistical analysis;
* Summaries of dominant lineage-state composition across experimental groups.

Together, these outputs allow users to examine whether a biomaterial primarily changes lineage direction, differentiation status, dominant lineage composition, or the distribution of cellular states within a population.

## Website modules

The BioLRAF interactive website is organized into four modules:

* **Reference Atlas**
  Displays the shared MSC lineage-state reference space and its osteogenic, chondrogenic, and adipogenic organization.

* **Differentiation Models**
  Presents public single-cell datasets representing experimentally induced MSC differentiation processes.

* **Validated Datasets**
  Shows datasets used to evaluate whether BioLRAF recovers expected lineage-associated responses.

* **Novel Biomaterials**
  Presents applications of BioLRAF to biomaterial-induced MSC responses, including comparisons of material architectures and interface designs.

## Interactive website and reproducible example

Explore the complete interactive BioLRAF resource:

[Open the BioLRAF interactive website](https://pangminmin.github.io/BioLRAF/)

Access the reproducible GSE226365 example:

[View the GSE226365 example directory](Example/GSE226365)

The GSE226365 example includes:

* An input cell-level gficf score table;
* A processed BioLRAF score table;
* Interactive three-dimensional visualization results;
* Fixed-view PNG outputs;
* A Jupyter notebook documenting the visualization workflow.

Open the corresponding example notebook:

[BioLRAF visualization notebook for GSE226365](Methods/BioLRAF_visualization_GSE226365.ipynb)

Additional method scripts are available in:

[BioLRAF Methods](Methods)
