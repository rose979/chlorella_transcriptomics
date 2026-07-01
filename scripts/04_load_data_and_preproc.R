# 04_load_data_and_preproc.R
# load the data and perform preprocessing steps

# Load count matrix
counts_raw <- read.table(count_file, header = TRUE, row.names = 1, sep = "\t",
                         check.names = FALSE)

# Load metadata 
metadata <- read.table(metadata_file, header = TRUE, row.names = 1, sep = ",")

# Ensure sample names match between count matrix and metadata
common_samples <- intersect(colnames(counts_raw), rownames(metadata))
counts_raw <- counts_raw[, common_samples]
metadata <- metadata[common_samples, ]

# Display data dimensions
cat("Count matrix dimensions:", dim(counts_raw), "\n")

# display metadata dimensions
cat("Metadata dimensions:", dim(metadata), "\n")

# display number of matching samples
cat("Number of matching samples:", length(common_samples), "\n")

# display first few rows of the count matrix
kable(head(counts_raw[, 1:min(6, ncol(counts_raw))]),
      caption = "First few rows of count matrix") %>%
  kable_styling(bootstrap_options = "striped")

# display the first rows of metadata table
kable(head(metadata), caption = "Sample metadata") %>%
  kable_styling(bootstrap_options = "striped")


