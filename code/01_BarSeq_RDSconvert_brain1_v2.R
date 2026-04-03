######## Brain #1 - 669594 ########
####### Export data from MATLAB to rds format to allow for transcriptomic analyses #######
# script order 01
# r4-base environment to run
# load libraries for the run
.libPaths("/opt/conda/envs/r4-base/lib/R/library/")
library(hdf5r)
library(Matrix)
library(SingleCellExperiment)
############################################################################################################################################################################################################
# Function to convert .mat to .rds file and specify input and output file names
convert = function() {
  sce = convert_v7_filtneurons ("/data/BARseq669594/filt_neurons_20240321.mat")
  dir.create("/results/BARseq_669594/", recursive = TRUE)  # Ensure directory exists
  saveRDS (sce, "/results/BARseq_669594/combined_neurons_clust_CCFv2.rds")    
}

#load and activate the conversion function
convert_v7_filtneurons = function(mat_file, mode = "r") {
  h5_file = H5File$new(mat_file, mode = mode)  # Added mode parameter to open read-only
  
  n_genes = length(readDataSet(h5_file[["filt_neurons/expmat/jc"]]))-1
  gene_refs = readDataSet(h5_file[["filt_neurons/genes"]])
  gene_ids = gene_refs$dereference()
  genes = lapply(gene_ids, readDataSet)
  genes = lapply(genes, intToUtf8)
  genes = unlist(genes[1:n_genes])
  
  sample_name = as.character(readDataSet(h5_file[["filt_neurons/id"]]))
  position = readDataSet(h5_file[["filt_neurons/pos"]])
  depth = readDataSet(h5_file[["filt_neurons/depth"]])
  angle = readDataSet(h5_file[["filt_neurons/angle"]])
  slice = as.character(readDataSet(h5_file[["filt_neurons/slice"]]))
  fov_position = readDataSet(h5_file[["filt_neurons/pos40x"]])
  is_barcoded = readDataSet(h5_file[["filt_neurons/is_barcoded"]])
  #batch_num = as.character(readDataSet(h5_file[["filt_neurons/batch_num"]]))
  CCF = readDataSet(h5_file[["filt_neurons/CCF"]])
  CCFano = readDataSet(h5_file[["filt_neurons/CCFano"]])
  
  # Sanity check for length consistency
  stopifnot (length(slice) == length(sample_name)) #(length(batch_num) == length(slice)),
  # Construct uid as "batch_slice_id"
  uid <- paste(slice, sample_name, sep = "_") #batch_num,
  
  
  expr = Matrix::sparseMatrix(
    i=readDataSet(h5_file[["filt_neurons/expmat/ir"]])+1,
    p=readDataSet(h5_file[["filt_neurons/expmat/jc"]]),
    x=as.vector(readDataSet(h5_file[["filt_neurons/expmat/data"]])),
    dims=c(length(sample_name), length(genes)),
    dimnames = list(sample_name, genes)
  )
  expr = t(expr)
  
  metadata = data.frame(
    slice = as.vector(slice),
    pos_x = position[,1],
    pos_y = position[,2],
    fov_x = fov_position[,1],
    fov_y = fov_position[,2],
    angle = as.vector(angle),
    depth_x = depth[,1],
    depth_y = depth[,2],
    barcode = as.vector(is_barcoded),
    #batch_num = as.vector(batch_num),
    uid = uid,
    CCF_AP = CCF[,1],
    CCF_DV = CCF[,2],
    CCF_ML = CCF[,3],
    CCFano = as.vector(CCFano),
    stringsAsFactors = FALSE
  )
  
  sce = SingleCellExperiment(assays = list(counts = expr), colData = metadata)
  return(sce)
}

#execute convert function
convert()

############################################################################################################################################################################################################
# Sanity check the converted file
sce <- readRDS('/results/BARseq_669594/combined_neurons_clust_CCFv2.rds')

#check parameters to ensure .mat to .rds conversion was correct
str(sce) #provides human readable file structure
dim(sce)
colnames(sce) #cells
rownames(sce) #genes
head(assay(sce)) #counts matrix
colData(sce) #cell metadata
rowData(sce) #gene metadata

# Ensure uid exists in colData
if (!"uid" %in% colnames(colData(sce))) {
  stop("uid is not found in colData of the SCE object!")
}
# Ensure uid is unique and matches the number of columns
if (length(unique(colData(sce)$uid)) != ncol(sce)) {
  stop("uid is not unique or does not match the number of columns in the SCE object!")
}
# Rename colnames using uid
colnames(sce) <- colData(sce)$uid
# Verify the changes
print(head(colnames(sce)))

#remove "empty" genes from further processing
gene_names <- rownames(sce)
unused_genes <- grepl("unused-", gene_names)
sce <- sce[!unused_genes, ]
dim(sce)

#check if there are duplicate genes and exclude those from HYB cycles
length(rownames(sce))
length(unique(rownames(sce)))
# Find duplicated names
duplicated_names <- rownames(sce)[duplicated(rownames(sce)) | duplicated(rownames(sce), fromLast = TRUE)]
# Get the indices of all duplicates
duplicate_indices <- which(rownames(sce) %in% duplicated_names)
# Print the duplicate names, their indices, and associated values
for (name in unique(duplicated_names)) {
  indices <- which(rownames(sce) == name)
  values <- sce[name, ]  # Extract values for the duplicate rows
  cat("Duplicate Name:", name, "\n")
  cat("Indices:", indices, "\n")
}
# Identify the rows to drop (higher index for each duplicate comes from HYB cycle which is blank)
rows_to_drop <- c()
for (name in unique(duplicated_names)) {
  indices <- which(rownames(sce) == name)
  cat("Duplicate Name:", name, "Indices:", indices, "Dropping Index:", max(indices), "\n")
  rows_to_drop <- c(rows_to_drop, max(indices))  # Keep the higher index
}

# Drop the rows with higher indices (only if duplicates exist)
if (length(rows_to_drop) > 0) {
  sce <- sce[-rows_to_drop, ]
}
dim(sce)
rm(duplicated_names, duplicate_indices, indices, gene_names, rows_to_drop, unused_genes, name)

saveRDS(sce, "/results/BARseq_669594/combined_neurons_clust_CCFv2_uid.rds")

############################################################################################################################################################################################################
# Identify total number of barcoded and non-barcoded cells 
table(sce[['barcode']]) 

#check what is the total number of cells with Dbh gene reads detected that also have barcodes 
# Extract count matrix
counts_matrix <- assay(sce, "counts")
# Ensure Dbh gene is present in the dataset
if (!"Dbh" %in% rownames(counts_matrix)) {
  stop("Dbh gene is not found in the dataset.")
}
# Extract barcode and duplicate information
barcode_info <- colData(sce)$barcode

# Identify cells where barcode == 1
barcoded_cells <- colnames(sce)[barcode_info == 1]
# Identify cells where barcode == 0
non_barcoded_cells <- colnames(sce)[barcode_info == 0]

# Extract Dbh expression levels for barcode == 1 cells
dbh_counts_barcoded <- counts_matrix["Dbh", barcoded_cells]
# Extract Dbh expression levels for barcode == 0 cells
dbh_counts_non_barcoded <- counts_matrix["Dbh", non_barcoded_cells]

# Count how many barcode 1 cells have Dbh expression > 1 and > 2
num_cells_dbh_above_1 <- sum(dbh_counts_barcoded > 1)
num_cells_dbh_above_2 <- sum(dbh_counts_barcoded > 2)
# Count how many barcode == 0 cells have Dbh expression > 1 and > 2
num_cells_dbh_above_1_non_barcode <- sum(dbh_counts_non_barcoded > 1)
num_cells_dbh_above_2_non_barcode <- sum(dbh_counts_non_barcoded > 2)

# Print results
cat("Total barcode 1 cells with Dbh expression > 1:", num_cells_dbh_above_1, "\n")
cat("Total barcode 1 cells with Dbh expression > 2:", num_cells_dbh_above_2, "\n")
cat("Total non-barcode cells with Dbh expression > 1:", num_cells_dbh_above_1_non_barcode, "\n")
cat("Total non-barcode cells with Dbh expression > 2:", num_cells_dbh_above_2_non_barcode, "\n")

############################################################################################################################################################################################################
# Identify and subset out only the cells with Dbh 2+ expression
non_zero_dbh_cells <- colnames(sce)[counts_matrix["Dbh", ] > 2]
# Subset the sce object to include only cells with Dbh 2+ expression
sce_filtered <- sce[, non_zero_dbh_cells]
ncol(sce_filtered)
# Save the filtered sce object to a file
saveRDS(sce_filtered, file = "/results/BARseq_669594/DBHfilteredneurons_clust_CCFv2_uid.rds")
