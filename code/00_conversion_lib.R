# Library functions for BARseq MAT-to-RDS conversion.
# Sourced by the per-subject 01_BarSeq_RDSconvert_brain*.R scripts.

.libPaths("/opt/conda/envs/r4-base/lib/R/library/")
library(hdf5r)
library(Matrix)
library(SingleCellExperiment)

PROCESS_NAME <- "processed_MAT2RDS"

# Read a v7.3 BARseq .mat file and return a SingleCellExperiment.
convert_v7_filtneurons <- function(mat_file, mode = "r") {
  h5_file <- H5File$new(mat_file, mode = mode)

  n_genes <- length(readDataSet(h5_file[["filt_neurons/expmat/jc"]])) - 1
  gene_refs <- readDataSet(h5_file[["filt_neurons/genes"]])
  gene_ids <- gene_refs$dereference()
  genes <- lapply(gene_ids, readDataSet)
  genes <- lapply(genes, intToUtf8)
  genes <- unlist(genes[1:n_genes])

  sample_name <- as.character(readDataSet(h5_file[["filt_neurons/id"]]))
  position <- readDataSet(h5_file[["filt_neurons/pos"]])
  depth <- readDataSet(h5_file[["filt_neurons/depth"]])
  angle <- readDataSet(h5_file[["filt_neurons/angle"]])
  slice <- as.character(readDataSet(h5_file[["filt_neurons/slice"]]))
  fov_position <- readDataSet(h5_file[["filt_neurons/pos40x"]])
  is_barcoded <- readDataSet(h5_file[["filt_neurons/is_barcoded"]])
  batch_num <- as.character(readDataSet(h5_file[["filt_neurons/batch_num"]]))
  CCF <- readDataSet(h5_file[["filt_neurons/CCF"]])
  CCFano <- readDataSet(h5_file[["filt_neurons/CCFano"]])

  stopifnot(length(batch_num) == length(slice), length(slice) == length(sample_name))
  uid <- paste(batch_num, slice, sample_name, sep = "_")

  expr <- Matrix::sparseMatrix(
    i = readDataSet(h5_file[["filt_neurons/expmat/ir"]]) + 1,
    p = readDataSet(h5_file[["filt_neurons/expmat/jc"]]),
    x = as.vector(readDataSet(h5_file[["filt_neurons/expmat/data"]])),
    dims = c(length(sample_name), length(genes)),
    dimnames = list(sample_name, genes)
  )
  expr <- t(expr)

  metadata <- data.frame(
    slice = as.vector(slice),
    pos_x = position[, 1],
    pos_y = position[, 2],
    fov_x = fov_position[, 1],
    fov_y = fov_position[, 2],
    angle = as.vector(angle),
    depth_x = depth[, 1],
    depth_y = depth[, 2],
    barcode = as.vector(is_barcoded),
    batch_num = as.vector(batch_num),
    uid = uid,
    CCF_AP = CCF[, 1],
    CCF_DV = CCF[, 2],
    CCF_ML = CCF[, 3],
    CCFano = as.vector(CCFano),
    stringsAsFactors = FALSE
  )

  SingleCellExperiment(assays = list(counts = expr), colData = metadata)
}

# End-to-end per-subject pipeline:
#   /data/<input_asset>/BARseq/combined_neurons_clust_CCFv2.mat
#     -> /results/<input_asset>_processed_MAT2RDS_<timestamp>/{
#          combined_neurons_clust_CCFv2.rds,
#          combined_neurons_clust_CCFv2_uid.rds,
#          DBHfilteredneurons_clust_CCFv2_uid.rds
#        }
convert_subject <- function(input_asset) {
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  out_dir <- file.path("/results", paste0(input_asset, "_", PROCESS_NAME, "_", timestamp))
  dir.create(out_dir, recursive = TRUE)
  mat_file <- file.path("/data", input_asset, "BARseq", "combined_neurons_clust_CCFv2.mat")

  cat("Input asset:", input_asset, "\n")
  cat("Mat file:   ", mat_file, "\n")
  cat("Output dir: ", out_dir, "\n\n")

  # Step 1: convert .mat to SCE and save
  sce <- convert_v7_filtneurons(mat_file)
  saveRDS(sce, file.path(out_dir, "combined_neurons_clust_CCFv2.rds"))

  # Sanity check the converted file
  sce <- readRDS(file.path(out_dir, "combined_neurons_clust_CCFv2.rds"))

  str(sce)
  print(dim(sce))
  print(colnames(sce))
  print(rownames(sce))
  print(head(assay(sce)))
  print(colData(sce))
  print(rowData(sce))

  if (!"uid" %in% colnames(colData(sce))) {
    stop("uid is not found in colData of the SCE object!")
  }
  if (length(unique(colData(sce)$uid)) != ncol(sce)) {
    stop("uid is not unique or does not match the number of columns in the SCE object!")
  }
  colnames(sce) <- colData(sce)$uid
  print(head(colnames(sce)))

  # Step 2: drop "unused-" placeholder genes
  unused_genes <- grepl("unused-", rownames(sce))
  sce <- sce[!unused_genes, ]
  print(dim(sce))

  # Step 3: drop duplicate HYB-cycle genes. For each gene name that appears
  # more than once in rownames, keep the lower-indexed row (the real
  # measurement) and drop the higher-indexed row (the blank HYB cycle). This
  # is the same logic that lived in the original brain 3 and brain 4 scripts,
  # consolidated here. The `if (length(rows_to_drop) > 0)` wrap is brain 4's
  # variant; for inputs that contain duplicates -- which both brains' actual
  # inputs do -- behavior is identical to either original.
  print(length(rownames(sce)))
  print(length(unique(rownames(sce))))
  duplicated_names <- rownames(sce)[duplicated(rownames(sce)) | duplicated(rownames(sce), fromLast = TRUE)]
  for (name in unique(duplicated_names)) {
    indices <- which(rownames(sce) == name)
    cat("Duplicate Name:", name, "\n")
    cat("Indices:", indices, "\n")
  }
  rows_to_drop <- integer(0)
  for (name in unique(duplicated_names)) {
    indices <- which(rownames(sce) == name)
    cat("Duplicate Name:", name, "Indices:", indices, "Dropping Index:", max(indices), "\n")
    rows_to_drop <- c(rows_to_drop, max(indices))
  }
  if (length(rows_to_drop) > 0) {
    sce <- sce[-rows_to_drop, ]
  }
  print(dim(sce))

  saveRDS(sce, file.path(out_dir, "combined_neurons_clust_CCFv2_uid.rds"))

  # Step 4: barcode + Dbh diagnostics
  print(table(sce[["barcode"]]))

  counts_matrix <- assay(sce, "counts")
  if (!"Dbh" %in% rownames(counts_matrix)) {
    stop("Dbh gene is not found in the dataset.")
  }

  barcode_info <- colData(sce)$barcode
  barcoded_cells <- colnames(sce)[barcode_info == 1]
  non_barcoded_cells <- colnames(sce)[barcode_info == 0]
  dbh_counts_barcoded <- counts_matrix["Dbh", barcoded_cells]
  dbh_counts_non_barcoded <- counts_matrix["Dbh", non_barcoded_cells]

  cat("Total barcode 1 cells with Dbh expression > 1:", sum(dbh_counts_barcoded > 1), "\n")
  cat("Total barcode 1 cells with Dbh expression > 2:", sum(dbh_counts_barcoded > 2), "\n")
  cat("Total non-barcode cells with Dbh expression > 1:", sum(dbh_counts_non_barcoded > 1), "\n")
  cat("Total non-barcode cells with Dbh expression > 2:", sum(dbh_counts_non_barcoded > 2), "\n")

  # Step 5: filter to Dbh > 2 and save
  sce_filtered <- sce[, counts_matrix["Dbh", ] > 2]
  print(ncol(sce_filtered))
  saveRDS(sce_filtered, file.path(out_dir, "DBHfilteredneurons_clust_CCFv2_uid.rds"))

  invisible(sce)
}
