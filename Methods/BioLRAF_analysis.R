#!/usr/bin/env Rscript

# BioLRAF example analysis: reference-query integration and gficf transformation
#
# This script merges one query scRNA-seq sample with the BioLRAF reference
# object, calculates a joint gficf representation, and exports the query sample
# as a cell-by-feature CSV matrix for downstream BioLRAF scoring or Jupyter
# visualization.
#
# Example:
# Rscript Methods/BioLRAF_analysis.R \
#   --reference Example/input/tri_int.rda \
#   --query Example/input/GSE198482_anno.rda \
#   --dataset GSE198482 \
#   --output Example/output/GSE198482_gficf_matrix.csv

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(gficf)
})

print_usage <- function() {
  cat(
    paste0(
      "Usage:\n",
      "  Rscript BioLRAF_analysis.R \\\n",
      "    --reference <reference_rda> \\\n",
      "    --query <query_rda> \\\n",
      "    --dataset <dataset_id> \\\n",
      "    --output <output_csv>\n\n",
      "Arguments:\n",
      "  --reference  RDA file containing the BioLRAF reference Seurat object.\n",
      "  --query      RDA file containing one annotated query Seurat object.\n",
      "  --dataset    Value in meta.data$dataset identifying the query sample.\n",
      "  --output     Output cell-by-feature gficf CSV file.\n",
      "  --help       Print this message.\n"
    )
  )
}

parse_arguments <- function(args) {
  if (length(args) == 0L || "--help" %in% args) {
    print_usage()
    quit(status = 0L)
  }

  valid_flags <- c("--reference", "--query", "--dataset", "--output")
  if (length(args) %% 2L != 0L) {
    stop("Each command-line flag must be followed by a value.", call. = FALSE)
  }

  supplied_flags <- args[seq(1L, length(args), by = 2L)]
  unknown_flags <- setdiff(supplied_flags, valid_flags)
  if (length(unknown_flags) > 0L) {
    stop(
      "Unknown argument(s): ", paste(unknown_flags, collapse = ", "),
      call. = FALSE
    )
  }

  values <- args[seq(2L, length(args), by = 2L)]
  names(values) <- sub("^--", "", supplied_flags)

  missing_flags <- setdiff(sub("^--", "", valid_flags), names(values))
  if (length(missing_flags) > 0L) {
    stop(
      "Missing required argument(s): --",
      paste(missing_flags, collapse = ", --"),
      call. = FALSE
    )
  }

  as.list(values)
}

load_seurat_object <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " file does not exist: ", path, call. = FALSE)
  }

  object_env <- new.env(parent = globalenv())
  loaded_names <- load(path, envir = object_env)
  is_seurat <- vapply(
    loaded_names,
    function(name) inherits(get(name, envir = object_env), "Seurat"),
    logical(1)
  )

  seurat_names <- loaded_names[is_seurat]
  if (length(seurat_names) != 1L) {
    stop(
      label, " RDA must contain exactly one Seurat object; found ",
      length(seurat_names), ".",
      call. = FALSE
    )
  }

  get(seurat_names[[1L]], envir = object_env)
}

join_rna_layers <- function(object) {
  if (!"RNA" %in% names(object@assays)) {
    stop("The merged Seurat object does not contain an RNA assay.", call. = FALSE)
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

message("[1/6] Loading reference and query Seurat objects")
reference_obj <- load_seurat_object(args$reference, "Reference")
query_obj <- load_seurat_object(args$query, "Query")

if (!"dataset" %in% colnames(query_obj@meta.data)) {
  message("      Query metadata has no 'dataset' column; assigning ", args$dataset)
  query_obj$dataset <- args$dataset
}

query_dataset_values <- unique(as.character(query_obj$dataset))
if (!args$dataset %in% query_dataset_values) {
  stop(
    "Dataset '", args$dataset, "' was not found in query meta.data$dataset. ",
    "Available value(s): ", paste(query_dataset_values, collapse = ", "),
    call. = FALSE
  )
}

if (!"celltype" %in% colnames(query_obj@meta.data)) {
  warning("Query metadata has no 'celltype' column; exporting NA values.")
  query_obj$celltype <- NA_character_
}

message("[2/6] Merging the BioLRAF reference and query sample")
merged_obj <- merge(
  x = reference_obj,
  y = query_obj,
  add.cell.ids = c("Reference", args$dataset),
  project = "BioLRAF_example"
)
merged_obj <- join_rna_layers(merged_obj)

rm(reference_obj, query_obj)
invisible(gc())

message("[3/6] Extracting the joined RNA count matrix")
counts_mat <- LayerData(
  object = merged_obj,
  assay = "RNA",
  layer = "counts"
)

if (is.null(rownames(counts_mat)) || is.null(colnames(counts_mat))) {
  stop("The RNA count matrix must contain gene and cell names.", call. = FALSE)
}

message(
  "      Joint matrix: ", nrow(counts_mat), " genes x ",
  ncol(counts_mat), " cells"
)

message("[4/6] Calculating the joint gficf representation")
gficf_result <- gficf(M = counts_mat)

if (is.null(gficf_result$gficf)) {
  stop("gficf() did not return a component named 'gficf'.", call. = FALSE)
}

gficf_mat <- gficf_result$gficf
if (!identical(dim(gficf_mat), dim(counts_mat))) {
  stop(
    "The gficf matrix dimensions do not match the input count matrix.",
    call. = FALSE
  )
}

if (is.null(rownames(gficf_mat))) {
  rownames(gficf_mat) <- rownames(counts_mat)
}
if (is.null(colnames(gficf_mat))) {
  colnames(gficf_mat) <- colnames(counts_mat)
}

message("[5/6] Extracting the query sample")
metadata <- merged_obj@meta.data
query_cells <- rownames(metadata)[
  as.character(metadata$dataset) == args$dataset
]
query_cells <- intersect(query_cells, colnames(gficf_mat))

if (length(query_cells) == 0L) {
  stop(
    "No cells were retained for dataset '", args$dataset, "'.",
    call. = FALSE
  )
}

# Subset before converting to a dense data frame to reduce memory use.
query_gficf <- gficf_mat[, query_cells, drop = FALSE]
query_matrix <- t(query_gficf)
estimated_gb <- prod(dim(query_matrix)) * 8 / 1024^3
message(
  "      Query matrix: ", nrow(query_matrix), " cells x ",
  ncol(query_matrix), " features"
)
message(
  "      Approximate dense matrix size: ",
  format(round(estimated_gb, 2), nsmall = 2), " GB"
)

query_df <- as.data.frame(as.matrix(query_matrix))
query_metadata <- metadata[rownames(query_df), , drop = FALSE]

output_df <- data.frame(
  cell_id = rownames(query_df),
  query_df,
  celltype = as.character(query_metadata$celltype),
  dataset = as.character(query_metadata$dataset),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

message("[6/6] Writing the query gficf matrix")
output_dir <- dirname(args$output)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

write.csv(
  output_df,
  file = args$output,
  row.names = FALSE,
  quote = TRUE
)

message("BioLRAF example preprocessing completed successfully.")
message("Output: ", normalizePath(args$output, mustWork = FALSE))

