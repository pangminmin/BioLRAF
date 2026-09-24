#!/usr/bin/env Rscript

# BioLRAF analysis example for GSE226365
#
# This script:
# 1. Reads the reference Seurat object, query Seurat object, and gene list
#    from RDS files.
# 2. Merges the reference and query objects.
# 3. Joins the Seurat v5 RNA count layers.
# 4. Runs gficf, PCA, UMAP, and single-cell GSEA.
# 5. Extracts and saves the GSE226365 query-cell results as an RDS file.
#
# Example:
# Rscript Methods/BioLRAF_analysis.R \
#   --reference Example/GSE226365/input/tri_int.rds \
#   --query Example/GSE226365/input/GSE226365_anno_final.rds \
#   --gene-list Example/GSE226365/input/bulk_gene_list.rds \
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
      "    --gene-list <gene_list_rds> \\\n",
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

message("[1/7] Reading reference, query, and gene-list RDS files")

reference_obj <- read_rds_checked(
  args$reference,
  "Reference"
)

query_obj <- read_rds_checked(
  args$query,
  "Query"
)

bulk_gene_list <- read_rds_checked(
  args[["gene-list"]],
  "Gene-list"
)

if (!inherits(reference_obj, "Seurat")) {
  stop(
    "The reference RDS file must contain a Seurat object.",
    call. = FALSE
  )
}

if (!inherits(query_obj, "Seurat")) {
  stop(
    "The query RDS file must contain a Seurat object.",
    call. = FALSE
  )
}

if (!is.list(bulk_gene_list) || length(bulk_gene_list) == 0L) {
  stop(
    "The gene-list RDS file must contain a non-empty list.",
    call. = FALSE
  )
}

if (!"dataset" %in% colnames(query_obj@meta.data)) {
  message(
    "      Query metadata has no 'dataset' column; assigning ",
    args$dataset
  )

  query_obj$dataset <- args$dataset
}

query_dataset_values <- unique(as.character(query_obj$dataset))

if (!args$dataset %in% query_dataset_values) {
  stop(
    "Dataset '",
    args$dataset,
    "' was not found in query metadata. Available value(s): ",
    paste(query_dataset_values, collapse = ", "),
    call. = FALSE
  )
}

if (!"celltype" %in% colnames(query_obj@meta.data)) {
  warning(
    "Query metadata has no 'celltype' column; exporting NA values."
  )

  query_obj$celltype <- NA_character_
}

# Add a source label before merging so that reference cells
# are not included in the final query-only output.
reference_obj$BioLRAF_source <- "reference"
query_obj$BioLRAF_source <- "query"

message("[2/7] Merging reference and query objects")

merged_obj <- merge(
  x = reference_obj,
  y = query_obj,
  add.cell.ids = c("reference", args$dataset),
  project = paste0("BioLRAF_", args$dataset)
)

rm(reference_obj, query_obj)
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

message("[6/7] Running single-cell GSEA")

gficf_data <- runScGSEA(
  data = gficf_data,
  geneID = "symbol",
  species = "human",
  pathway.list = bulk_gene_list,
  nmf.k = 100,
  rescale = "none"
)

if (is.null(gficf_data$scgsea$x)) {
  stop(
    "runScGSEA() did not return a result in gficf_data$scgsea$x.",
    call. = FALSE
  )
}

# Convert the scGSEA result to a matrix and align its cell dimension
# with the metadata from the merged Seurat object.
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
  # If cell names are absent, use the merged-object order
  # after checking that the number of rows matches.
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
    scgsea_mat$dataset == args$dataset
)

query_result <- scgsea_mat[
  query_rows,
  ,
  drop = FALSE
]

if (nrow(query_result) == 0L) {
  stop(
    "No query cells were retained for dataset '",
    args$dataset,
    "'.",
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

output_dir <- dirname(args$output)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

saveRDS(
  query_result,
  file = args$output
)

message("BioLRAF analysis completed successfully.")
message(
  "Output: ",
  normalizePath(args$output, mustWork = FALSE)
)
