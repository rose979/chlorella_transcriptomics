# 03_conf_analysis.R
# configuration of anaysis variables

# set Analysis parameters
min_count <- 5         # Minimum total counts per gene
min_samples <- 3        # Minimum samples with counts above threshold
padj_threshold <- 0.05  # Adjusted p-value threshold for significance
lfc_threshold <- 1      # Log2 fold change threshold

print(paste("Using min_count =", min_count, ", min_samples =", min_samples,
            ", padj_threshold =", padj_threshold, ", lfc_threshold =", lfc_threshold))

# set DESeq analysis parameters for differential expression analysis
design_formula <- ~ ammonium # Replicates are irrelevant so not included


########################################################################################
# maybe write a function, that selects the configuration of deseq out of the listed cases?
# Contrasts for differential expression
contrasts_list <- list("treatment_vs_control" = c("condition", "treatment", "control"))
 
contrasts_list <- list("ammonium" = "ammonium")  # Contrast for ammonium concentration if numeric variable
# With DESeq2, a numeric covariate reports the log2FC per one unit of that variable.

# If you code ammonium as 0 and 600, the coefficient is “per 1 unit ammonium,” so the beta will look ~0 (because the real effect is spread over 600 units). That’s why the log2FCs look tiny.
contrasts_list <- list("ammonium" = c("ammonium", "600", "0"))
########################################################################################