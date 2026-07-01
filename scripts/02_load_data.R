# 02_load_data.R
# load required files for analysis

# File paths, load count matrice and metadata table
count_file <- "C:/Users/fraro/Documents/Bioinformatics/chlorella_transcriptomics/data/processed/counts_genes.tsv" # Count matrix file
metadata_file <- "C:/Users/fraro/Documents/Bioinformatics/chlorella_transcriptomics/data/processed/metadata_C.csv"   # Sample metadata file

# functional annotation using predicted protein mappings
functional_annotation_file <- "C:/Users/fraro/Documents/Bioinformatics/chlorella_transcriptomics/data/processed/DSW-C_consensus_annotation.tsv"
go_annotation_file <- "C:/Users/fraro/Documents/Bioinformatics/chlorella_transcriptomics/data/processed/descriptions-n-goterms.txt" 

#########################################################################################
#kegg_annotation_file <- "path/to/kegg_annotations.txt"  # kegg or cog annotation?
cog_annotation_file <- "C:/Users/fraro/Documents/Bioinformatics/chlorella_transcriptomics/data/processed/cog.txt"
##########################################################################################

# COG annotations from https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/ 03.09.2025 
# COG category definition 
cog_def_file <- "C:/Users/fraro/Documents/Bioinformatics/chlorella_transcriptomics/data/processed/cog-24.fun.tab"

# COG term mapping and description
cog_descr_file <- "C:/Users/fraro/Documents/Bioinformatics/chlorella_transcriptomics/data/processed/cog-24.def.tab" 
