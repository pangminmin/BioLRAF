#!/usr/bin/env Rscript

# BioLRAF analysis example for GSE226365
#
# Input files:
#   Reference_data.rds   Reference Seurat object
#   GSE226365_data.rds   GSE226365 query Seurat object
#   gene_set.rds         Gene-set list
#
# Workflow:
# 1. Read the reference object, query object, and gene-set list from RDS files.
# 2. Merge the reference and query objects.
# 3. Join the Seurat v5 RNA count layers.
# 4. Run gficf, PCA, UMAP, and single-cell GSEA.
# 5. Extract and save the GSE226365 query-cell results as an RDS file.
#
# Example:
# Rscript Methods/BioLRAF_analysis.R \
#   --reference Example/GSE226365/input/Reference_data.rds \
#   --query Example/GSE226365/input/GSE226365_data.rds \
#   --gene-list Example/GSE226365/input/gene_set.rds \
#   --dataset GSE226365 \
#   --output Example/GSE226365/processed/GSE226365_mat.rds

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(gficf)
})

print_usage <- function() {
  cat(
    paste0(
      "Usage:\n",
      "  Rscript Methods/BioLRAF_analysis.R \\\n",
      "    --reference <reference_rds> \\\n",
      "    --query <query_rds> \\\n",
      "    --gene-list <gene_set_rds> \\\n",
      "    --dataset <dataset_id> \\\n",
      "    --output <output_rds>\n\n",
      "Arguments:\n",
      "  --reference  RDS file containing the reference Seurat object.\n",
      "  --query      RDS file containing the annotated query Seurat object.\n",
      "  --gene-list  RDS file containing the gene-set list.\n",
      "  --dataset    Dataset label used in query metadata.\n",
      "  --output     Output RDS file for the query-cell results.\n",
      "  --help       Print this message.\n"
    )
  )
}

parse_arguments <- function(args) {
  if (length(args) == 0L || "--help" %in% args) {
    print_usage()
    quit(status = 0L)
  }

  valid_flags <- c(
    "--reference",
    "--query",
    "--gene-list",
    "--dataset",
    "--output"
  )

  if (length(args) %% 2L != 0L) {
    stop(
      "Each command-line flag must be followed by a value.",
      call. = FALSE
    )
  }

  supplied_flags <- args[seq(1L, length(args), by = 2L)]
  unknown_flags <- setdiff(supplied_flags, valid_flags)

  if (length(unknown_flags) > 0L) {
    stop(
      "Unknown argument(s): ",
      paste(unknown_flags, collapse = ", "),
      call. = FALSE
    )
  }

  values <- args[seq(2L, length(args), by = 2L)]
  names(values) <- sub("^--", "", supplied_flags)

  required_arguments <- sub("^--", "", valid_flags)
  missing_arguments <- setdiff(required_arguments, names(values))

  if (length(missing_arguments) > 0L) {
    stop(
      "Missing required argument(s): --",
      paste(missing_arguments, collapse = ", --"),
      call. = FALSE
    )
  }

  as.list(values)
}

read_rds_checked <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " file does not exist: ", path, call. = FALSE)
  }

  readRDS(path)
}

join_rna_layers <- function(object) {
  if (!"RNA" %in% names(object@assays)) {
    stop(
      "The merged Seurat object does not contain an RNA assay.",
      call. = FALSE
    )
  }

  join_layers <- get0(
    "JoinLayers",
    envir = asNamespace("SeuratObject"),
    mode = "function"
  )

  if (is.null(join_layers)) {
    stop(
      "JoinLayers() is unavailable. Install SeuratObject version 5 or later.",
      call. = FALSE
    )
  }

  join_layers(object, assay = "RNA")
}

args <- parse_arguments(commandArgs(trailingOnly = TRUE))

message("[1/7] Reading input RDS files")

reference_data <- read_rds_checked(
  args$reference,
  "Reference"
)

query_data <- read_rds_checked(
  args$query,
  "Query"
)

gene_set <- read_rds_checked(
  args[["gene-list"]],
  "Gene-set"
)

if (!inherits(reference_data, "Seurat")) {
  stop(
    "Reference_data.rds must contain a Seurat object.",
    call. = FALSE
  )
}

if (!inherits(query_data, "Seurat")) {
  stop(
    "GSE226365_data.rds must contain a Seurat object.",
    call. = FALSE
  )
}

if (!is.list(gene_set) || length(gene_set) == 0L) {
  stop(
    "gene_set.rds must contain a non-empty gene-set list.",
    call. = FALSE
  )
}

if (!"dataset" %in% colnames(query_data@meta.data)) {
  message(
    "      Query metadata has no 'dataset' column; assigning GSE226365."
  )

  query_data$dataset <- "GSE226365"
}

query_dataset_values <- unique(
  as.character(query_data$dataset)
)

if (!"GSE226365" %in% query_dataset_values) {
  stop(
    "The query metadata does not contain dataset label 'GSE226365'. ",
    "Available value(s): ",
    paste(query_dataset_values, collapse = ", "),
    call. = FALSE
  )
}

if (!"celltype" %in% colnames(query_data@meta.data)) {
  warning(
    "Query metadata has no 'celltype' column; exporting NA values."
  )

  query_data$celltype <- NA_character_
}

# Mark data provenance before merging, so only query cells are exported.
reference_data$BioLRAF_source <- "reference"
query_data$BioLRAF_source <- "query"

message("[2/7] Merging reference and GSE226365 query objects")

merged_obj <- merge(
  x = reference_data,
  y = query_data,
  add.cell.ids = c("reference", "GSE226365"),
  project = "BioLRAF_GSE226365"
)

rm(reference_data, query_data)
invisible(gc())

message("[3/7] Joining RNA count layers")

merged_obj <- join_rna_layers(merged_obj)

message("[4/7] Extracting RNA counts")

counts_mat <- LayerData(
  object = merged_obj,
  assay = "RNA",
  layer = "counts"
)

if (is.null(rownames(counts_mat)) || is.null(colnames(counts_mat))) {
  stop(
    "The RNA count matrix must contain gene and cell names.",
    call. = FALSE
  )
}

message(
  "      Joint matrix: ",
  nrow(counts_mat),
  " genes x ",
  ncol(counts_mat),
  " cells"
)

message("[5/7] Running gficf, PCA, and UMAP")

gficf_data <- gficf(
  M = counts_mat
)

gficf_data <- runPCA(
  data = gficf_data,
  dim = 10,
  use.odgenes = TRUE
)

gficf_data <- runReduction(
  data = gficf_data,
  reduction = "umap",
  nt = 2,
  n_neighbors = 150
)

message("[6/7] Running single-cell GSEA with gene_set")

gficf_data <- runScGSEA(
  data = gficf_data,
  geneID = "symbol",
  species = "human",
  pathway.list = gene_set,
  nmf.k = 100,
  rescale = "none"
)

if (is.null(gficf_data$scgsea$x)) {
  stop(
    "runScGSEA() did not return a result in gficf_data$scgsea$x.",
    call. = FALSE
  )
}

# Convert scGSEA results to a matrix and align the cell dimension
# with metadata from the merged Seurat object.
scgsea_result <- as.matrix(gficf_data$scgsea$x)
metadata <- merged_obj@meta.data
metadata_cell_ids <- rownames(metadata)

if (
  !is.null(rownames(scgsea_result)) &&
  all(rownames(scgsea_result) %in% metadata_cell_ids)
) {
  # Cells are represented by rows.
  score_mat <- scgsea_result
  metadata_use <- metadata[rownames(score_mat), , drop = FALSE]

} else if (
  !is.null(colnames(scgsea_result)) &&
  all(colnames(scgsea_result) %in% metadata_cell_ids)
) {
  # Cells are represented by columns; transpose to make cells rows.
  score_mat <- t(scgsea_result)
  metadata_use <- metadata[rownames(score_mat), , drop = FALSE]

} else if (nrow(scgsea_result) == nrow(metadata)) {
  # If cell names are absent, use merged-object order after checking
  # that the number of rows matches.
  score_mat <- scgsea_result
  rownames(score_mat) <- metadata_cell_ids
  metadata_use <- metadata

} else if (ncol(scgsea_result) == nrow(metadata)) {
  # If cell names are absent and cells are columns, transpose the matrix.
  score_mat <- t(scgsea_result)
  rownames(score_mat) <- metadata_cell_ids
  metadata_use <- metadata

} else {
  stop(
    "Could not align scGSEA results with cell metadata. ",
    "Check the dimensions and cell names of gficf_data$scgsea$x.",
    call. = FALSE
  )
}

scgsea_mat <- as.data.frame(
  score_mat,
  check.names = FALSE
)

scgsea_mat$dataset <- as.character(metadata_use$dataset)
scgsea_mat$celltype <- as.character(metadata_use$celltype)
scgsea_mat$BioLRAF_source <- as.character(
  metadata_use$BioLRAF_source
)

query_rows <- (
  scgsea_mat$BioLRAF_source == "query" &
    scgsea_mat$dataset == "GSE226365"
)

query_result <- scgsea_mat[
  query_rows,
  ,
  drop = FALSE
]

if (nrow(query_result) == 0L) {
  stop(
    "No GSE226365 query cells were retained.",
    call. = FALSE
  )
}

# This marker was used to select query cells and is not needed in the output.
query_result$BioLRAF_source <- NULL

message(
  "      Exporting ",
  nrow(query_result),
  " query cells and ",
  ncol(query_result) - 2L,
  " score column(s)"
)

output_path <- args$output
output_dir <- dirname(output_path)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

saveRDS(
  query_result,
  file = output_path
)

message("BioLRAF analysis completed successfully.")
message(
  "Output: ",
  normalizePath(output_path, mustWork = FALSE)
)
