######## TEMPORARY: validate refactored output against the previously-published reference asset ########
# This script is included only to confirm the refactor produces RDS files equivalent to
# the prior published output. It will be removed (along with the reference asset attachment)
# in a follow-up cleanup commit before the final release.

.libPaths("/opt/conda/envs/r4-base/lib/R/library/")
suppressPackageStartupMessages(library(SingleCellExperiment))

REFERENCE_MOUNT <- "/data/BARseq_MATtoRDSfiles_brain3_brain4"

compare_rds <- function(new_path, ref_path) {
  cat("Comparing:\n")
  cat("  new: ", new_path, "\n", sep = "")
  cat("  ref: ", ref_path, "\n", sep = "")

  new_obj <- readRDS(new_path)
  ref_obj <- readRDS(ref_path)

  if (identical(new_obj, ref_obj)) {
    cat("  IDENTICAL\n\n")
    return(TRUE)
  }

  cat("  DIFFER -- drilling down:\n")
  cat("    class:       new=", paste(class(new_obj), collapse = ","),
      "  ref=", paste(class(ref_obj), collapse = ","), "\n", sep = "")
  cat("    dim:         ", identical(dim(new_obj), dim(ref_obj)),
      "  (new=", paste(dim(new_obj), collapse = "x"),
      "  ref=", paste(dim(ref_obj), collapse = "x"), ")\n", sep = "")
  cat("    colnames:    ", identical(colnames(new_obj), colnames(ref_obj)), "\n", sep = "")
  cat("    rownames:    ", identical(rownames(new_obj), rownames(ref_obj)), "\n", sep = "")

  cd_new <- as.data.frame(colData(new_obj))
  cd_ref <- as.data.frame(colData(ref_obj))
  cat("    colData:     ", identical(cd_new, cd_ref), "\n", sep = "")
  if (!identical(cd_new, cd_ref)) {
    common_cols <- intersect(colnames(cd_new), colnames(cd_ref))
    for (col in common_cols) {
      cat("      colData$", col, ": ", identical(cd_new[[col]], cd_ref[[col]]), "\n", sep = "")
    }
    new_only <- setdiff(colnames(cd_new), colnames(cd_ref))
    ref_only <- setdiff(colnames(cd_ref), colnames(cd_new))
    if (length(new_only) > 0) cat("      colData new only: ", paste(new_only, collapse = ","), "\n", sep = "")
    if (length(ref_only) > 0) cat("      colData ref only: ", paste(ref_only, collapse = ","), "\n", sep = "")
  }

  rd_new <- as.data.frame(rowData(new_obj))
  rd_ref <- as.data.frame(rowData(ref_obj))
  cat("    rowData:     ", identical(rd_new, rd_ref), "\n", sep = "")
  if (!identical(rd_new, rd_ref)) {
    cat("      rowData new ncol=", ncol(rd_new), " cols=(", paste(colnames(rd_new), collapse = ","), ")\n", sep = "")
    cat("      rowData ref ncol=", ncol(rd_ref), " cols=(", paste(colnames(rd_ref), collapse = ","), ")\n", sep = "")
  }

  cat("    assayNames:  new=(", paste(assayNames(new_obj), collapse = ","),
      ")  ref=(", paste(assayNames(ref_obj), collapse = ","),
      ")  identical=", identical(assayNames(new_obj), assayNames(ref_obj)), "\n", sep = "")

  meta_eq <- identical(metadata(new_obj), metadata(ref_obj))
  cat("    metadata:    ", meta_eq, "\n", sep = "")
  if (!meta_eq) {
    cat("      metadata new:\n")
    print(metadata(new_obj))
    cat("      metadata ref:\n")
    print(metadata(ref_obj))
  }

  counts_new <- assay(new_obj, "counts")
  counts_ref <- assay(ref_obj, "counts")
  cat("    counts:      ", identical(counts_new, counts_ref), "\n", sep = "")
  if (!identical(counts_new, counts_ref)) {
    cat("      counts class: new=", paste(class(counts_new), collapse = ","),
        "  ref=", paste(class(counts_ref), collapse = ","), "\n", sep = "")
    cat("      counts dim:   ", identical(dim(counts_new), dim(counts_ref)), "\n", sep = "")
    if (identical(dim(counts_new), dim(counts_ref))) {
      n_diff <- sum(as.matrix(counts_new) != as.matrix(counts_ref))
      cat("      n_differing_elements: ", n_diff, "\n", sep = "")
    }
  }

  # Final fallback: all.equal walks the full structure and reports textual diffs
  ae <- tryCatch(all.equal(new_obj, ref_obj),
                 error = function(e) paste("all.equal error:", conditionMessage(e)))
  cat("    all.equal:\n")
  if (isTRUE(ae)) {
    cat("      TRUE (all.equal sees no semantic difference -- divergence is in object identity / attributes only)\n")
  } else {
    for (msg in ae) cat("      ", msg, "\n", sep = "")
  }

  cat("\n")
  return(FALSE)
}

cat("=== Validation against reference asset ===\n")
cat("Reference mount: ", REFERENCE_MOUNT, "\n\n", sep = "")

if (!dir.exists(REFERENCE_MOUNT)) {
  cat("Reference asset not mounted at ", REFERENCE_MOUNT,
      " - skipping validation.\n", sep = "")
} else {
  ref_files <- list.files(REFERENCE_MOUNT, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  cat("Found ", length(ref_files), " reference RDS file(s):\n", sep = "")
  for (p in ref_files) cat("  ", p, "\n", sep = "")
  cat("\n")

  new_files <- list.files("/results", pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  new_files <- new_files[grepl("_processed_MAT2RDS_", new_files)]
  cat("Found ", length(new_files), " new RDS file(s):\n", sep = "")
  for (p in new_files) cat("  ", p, "\n", sep = "")
  cat("\n")

  results <- list()
  for (new_path in new_files) {
    m <- regmatches(new_path, regexpr("\\d{6}", new_path))
    subject <- if (length(m) > 0) m else NA_character_
    if (is.na(subject)) {
      cat("WARN: could not extract 6-digit subject id from ", new_path, " -- skipping\n", sep = "")
      next
    }

    base <- basename(new_path)
    ref_candidates <- ref_files[basename(ref_files) == base & grepl(subject, ref_files)]

    if (length(ref_candidates) == 0) {
      cat("WARN: no reference match for subject ", subject, " basename ", base, " -- skipping\n", sep = "")
      next
    }
    if (length(ref_candidates) > 1) {
      cat("WARN: multiple reference matches for subject ", subject, " basename ", base, ":\n", sep = "")
      for (p in ref_candidates) cat("  ", p, "\n", sep = "")
      cat("Using first.\n")
    }

    ok <- compare_rds(new_path, ref_candidates[1])
    results[[length(results) + 1]] <- list(new_path = new_path, ref_path = ref_candidates[1], ok = ok)
  }

  cat("=== Summary ===\n")
  n_pass <- sum(sapply(results, function(r) r$ok))
  n_fail <- length(results) - n_pass
  cat(n_pass, " IDENTICAL, ", n_fail, " DIFFER\n", sep = "")
  if (n_fail > 0) {
    cat("\nMismatches:\n")
    for (r in results) {
      if (!r$ok) cat("  - ", r$new_path, "\n", sep = "")
    }
  }
}
