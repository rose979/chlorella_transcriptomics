# Load required libraries
if (!require(stringr, quietly = TRUE)) {
  install.packages("stringr")  
  library(stringr)
}

# Function to parse StringTie GTF file (Base R version - no dplyr dependency)
parse_stringtie_gtf <- function(gtf_file) {
  
  # Read the GTF file, skipping comment lines
  cat("Reading GTF file...\n")
  gtf_lines <- readLines(gtf_file)
  gtf_lines <- gtf_lines[!grepl("^#", gtf_lines)]
  
  if (length(gtf_lines) == 0) {
    stop("No valid GTF lines found in file")
  }
  
  # Parse into data frame
  gtf_df <- read.table(text = gtf_lines, sep = "\t", stringsAsFactors = FALSE,
                       col.names = c("seqname", "source", "feature", "start", 
                                     "end", "score", "strand", "frame", "attribute"))
  
  cat(sprintf("Loaded %d GTF entries\n", nrow(gtf_df)))
  
  # Filter for transcript entries
  transcript_rows <- gtf_df[gtf_df$feature == "transcript", ]
  if (nrow(transcript_rows) == 0) {
    stop("No transcript entries found in GTF file")
  }
  
  cat(sprintf("Found %d transcript entries\n", nrow(transcript_rows)))
  
  # Function to extract specific attribute from the attribute string
  extract_attributes <- function(attribute)
  {
    attribute.list <- unlist(strsplit(attribute, "; +"))
    unlist(lapply(attribute.list, function(x) {value=strsplit(x, " ")[[1]][2]; names(value)=strsplit(x, " ")[[1]][1]; value}))
  }
  
  all.attributes <- lapply(transcript_rows$attribute, extract_attributes)
  all.attributes <- do.call(rbind, lapply(all.attributes, function(x) {
    attr_names <- c("gene_id", "transcript_id", "reference_id", "ref_gene_id", "cov", "FPKM", "TPM")
    attr_values <- sapply(attr_names, function(name) if (name %in% names(x)) x[[name]] else NA)
    return(attr_values)
  }))
  
  result <- cbind(dplyr::select(transcript_rows, !attribute), all.attributes)

  return(result)
}

# Function to create mapping tables (Base R version)
create_transcript_mappings <- function(transcript_data) {
  
  # Create basic mapping table (remove duplicates)
  mapping_table <- transcript_data[, c("transcript_id", "gene_id", "reference_id", "ref_gene_id")]
  mapping_table <- unique(mapping_table)
  
  # Create named vectors for quick lookups
  mappings <- list(
    # StringTie transcript ID to StringTie gene ID
    transcript_to_gene = setNames(mapping_table$gene_id, mapping_table$transcript_id),
    
    # StringTie transcript ID to reference transcript ID
    transcript_to_ref = setNames(mapping_table$reference_id, mapping_table$transcript_id),
    
    # StringTie transcript ID to reference gene ID
    transcript_to_ref_gene = setNames(mapping_table$ref_gene_id, mapping_table$transcript_id),
    
    # StringTie gene ID to reference gene ID
    gene_to_ref_gene = setNames(mapping_table$ref_gene_id, mapping_table$gene_id),
    
    # Reference transcript ID to StringTie transcript ID (reverse lookup)
    ref_to_transcript = setNames(mapping_table$transcript_id, mapping_table$reference_id),
    
    # Reference gene ID to StringTie gene ID (reverse lookup)
    ref_gene_to_gene = setNames(mapping_table$gene_id, mapping_table$ref_gene_id)
  )
  
  # Remove NA mappings
  mappings <- lapply(mappings, function(x) x[!is.na(x) & !is.na(names(x))])
  
  return(list(table = mapping_table, vectors = mappings))
}

# Main function to process StringTie GTF file
process_stringtie_gtf <- function(gtf_file) {
  
  # Parse the GTF file
  transcript_data <- parse_stringtie_gtf(gtf_file)
  
  # Create mappings
  mappings <- create_transcript_mappings(transcript_data)
  
  # Print summary
  summarize_mappings(transcript_data, mappings)
  
  # Return everything
  return(list(
    data = transcript_data,
    mapping_table = mappings$table,
    lookup_vectors = mappings$vectors
  ))
}

# Example usage:
# result <- process_stringtie_gtf("path/to/transcripts.gtf")
# 
# # Access the data
# transcript_data <- result$data
# mapping_table <- result$mapping_table
# lookups <- result$lookup_vectors
# 
# # Example lookups:
# lookups$transcript_to_gene["STRG.1.1"]      # Get StringTie gene ID
# lookups$transcript_to_ref["STRG.1.1"]       # Get reference transcript ID
# lookups$transcript_to_ref_gene["STRG.1.1"]  # Get reference gene ID
# 
# # Get all transcripts for a reference gene
# ref_gene <- "83844-DSWC.00g000010-v1.0.a2"
# transcripts_for_gene <- mapping_table[mapping_table$ref_gene_id == ref_gene & !is.na(mapping_table$ref_gene_id), ]$transcript_id

# Alternative: Using rtracklayer (if you prefer - requires Bioconductor)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("rtracklayer")

parse_stringtie_with_rtracklayer <- function(gtf_file) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("rtracklayer package not available. Install with: BiocManager::install('rtracklayer')")
  }
  
  # Load GTF file
  gtf <- rtracklayer::import(gtf_file)
  
  # Convert to data frame
  gtf_df <- as.data.frame(gtf)
  
  # Filter for transcripts and create mapping
  transcript_mapping <- gtf_df[gtf_df$type == "transcript", c("transcript_id", "gene_id", "reference_id", "ref_gene_id")]
  transcript_mapping <- unique(transcript_mapping)
  
  return(transcript_mapping)
}

# Usage example for rtracklayer approach:
# mapping_rt <- parse_stringtie_with_rtracklayer("path/to/transcripts.gtf")