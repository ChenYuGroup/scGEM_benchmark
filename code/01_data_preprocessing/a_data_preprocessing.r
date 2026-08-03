# ================================
# Single-cell RNA-seq data preprocessing
# Data loading + filtering + Linnorm normalization + SAVER imputation
# ================================

rm(list = ls())

# Load required packages
library(SingleCellExperiment)
library(Linnorm)
library(SAVER)


# ================================
# 1. Load raw expression matrix
# ================================

# Rows represent genes, columns represent cells
cell_data <- read.csv(
  "D:\\9_22_dropout\\GSE81861_Cell_Line_FPKM.csv"
)

# Remove Ensembl gene version numbers
# Example: ENSG000001234.5 -> ENSG000001234
cell_data$X <- sub(
  ".*(ENSG[0-9]+)\\..*",
  "\\1",
  cell_data$X
)

# Set gene IDs as row names
rownames(cell_data) <- cell_data$X

# Remove gene ID column
cell_data <- cell_data[, -1]


# Convert to expression matrix
raw_expression <- as.matrix(cell_data)


# ================================
# 2. Create SingleCellExperiment object
# ================================

sce <- SingleCellExperiment(
  assays = list(counts = raw_expression)
)


# ================================
# 3. Quality control and gene filtering
# ================================

# Filter genes and cells based on expression level
# keep1: retain genes with average expression > 1 among expressed cells
# keep2: retain genes detected in more than 10 cells
# keep3: retain cells with more than 100 detected genes

gene_filter <- function(sce){

  keep1 <- apply(
    counts(sce),
    1,
    function(x) mean(x[x > 0]) > 1
  )

  keep2 <- rowSums(
    counts(sce) > 0
  ) > 10

  keep3 <- colSums(
    counts(sce) > 0
  ) > 100


  sce <- sce[(keep1 & keep2), keep3]

  return(sce)
}


# Apply filtering
sce <- gene_filter(sce)


# Extract filtered expression matrix
filtered_counts <- counts(sce)

dim(filtered_counts)



# ================================
# 4. Linnorm normalization
# ================================

# Perform Linnorm normalization
# The output is a normalized expression matrix
linnorm_data <- Linnorm(
  filtered_counts
)

linnorm_data <- as.matrix(linnorm_data)



# ================================
# 5. SAVER imputation
# ================================

# Perform gene expression imputation using SAVER
# estimates.only = TRUE returns only the estimated expression matrix

saver_imputed_data <- saver(
  filtered_counts,
  ncores = 8,
  size.factor = 1,
  estimates.only = TRUE
)


# Convert output to matrix format
saver_imputed_data <- as.matrix(
  saver_imputed_data
)



# ================================
# 6. Save processed data
# ================================

save(
  filtered_counts,
  linnorm_data,
  saver_imputed_data,
  file = "D:\\9_22_dropout\\Linnorm_SAVER_processed_data.RData"
)