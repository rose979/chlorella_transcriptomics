# 01_load_packages.R
# Initial script to load required packages and loads them in R session

package_names <- c("DESeq2", "edgeR", "limma", "dplyr", "tidyr", "tibble", "ggplot2", 
                   "pheatmap", "RColorBrewer", "ggrepel", "plotly", "corrplot",
                   "clusterProfiler", "enrichplot", "topGO", "GenomicFeatures",
                   "rtracklayer", "Biostrings", "DT", "knitr", "kableExtra")


# create a function that installs packages if they are not already installed, uses BioCManager,
# because it handles package reposotories dynamically

install_if_missing <- function(packages) {
  # 1. Ensure BiocManager is installed, as it handles CRAN, Bioc, and GitHub
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("Installing BiocManager first...")
    install.packages("BiocManager")
  }
  
  # 2. Identify which packages are actually missing
  missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  
  # 3. Install missing packages using BiocManager's smart routing
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    
    # update = FALSE prevents it from trying to update your entire R library
    # ask = FALSE prevents interactive prompts blocking the script
    BiocManager::install(missing_packages, update = FALSE, ask = FALSE)
  } else {
    message("All packages are already installed.")
  }
}

# run the function on the package vector
install_if_missing(package_names)

# load the stringtie function, important for parsing information about transcripts into R format
source("stringtie_parser.R")

# load the required packages into the session, use a loop that iterates over the package names
for (package in package_names) {
  require(package, character.only = TRUE)
}

