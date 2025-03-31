# Project: ironProjectAnnaWawer
# Purpose: This script will import the RNAseq counts and metadata, plot library size, filter counts, generate DGEList and create MDS plot.
# Author: Melanie Smith
# Date Created: 2025-03-25
# Last Modified: 2025-03-25

# Script Overview ----
# This script does the following:
# 1. Load necessary libraries
# 2. Set up project directories
# 3. Import data
# 4. Perform data analysis/processing
# 5. Generate outputs/visualisations

# Dependency Management ----
# Use renv for package management (recommended)
# Run once at project start: 
# install.packages("renv")
# renv::init()

# Library Loading ----
# Always use library() or pacman for loading
# Suppress startup messages with suppressPackageStartupMessages()
# if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,    # Data manipulation and visualization
  ggrepel,      # not included in "tidyverse
  edgeR
)



# Project Directory Setup ----
project_root <- "D:/VM_Projects/ironReceptor_placenta" # needs the "windows style" directory path (D:)
dir.create(file.path(project_root, "output"), showWarnings = FALSE)

# Path Variables ----
gen_data_dir <- "D:/VM_Projects/dataFolder" # this is where I keep files to be used across multiple projects
data_dir <- file.path(project_root, "rawData")
output_dir <- file.path(project_root, "output")

# File variables
SAGCQA0492_readCounts_file <- file.path(data_dir, "SAGCQA0492_readCounts.txt")
jointMetaData_file <- file.path(data_dir, "20230405_jointMetaData.csv")
gene_info_file <- file.path(gen_data_dir, "gencode_v29_gene_id_symbol_chr_biotype.csv")

output_librarySize_file <- file.path(output_dir, "rawLibrarySize.pdf")

# Logging and Error Handling ----
# Optional but recommended for tracking
# log_file <- file.path(output_dir, paste0("log_", Sys.Date(), ".txt"))
# sink(log_file, append = TRUE, split = TRUE)

# Configuration ----
# Set global options
options(
  scipen = 999,           # Disable scientific notation
  stringsAsFactors = FALSE # Prevent automatic factor conversion
)

# Function Definitions ----
# Define any custom functions here
clean_data <- function(df) {
  # Example data cleaning function
  df %>%
    janitor::clean_names() %>%
    filter(!is.na(some_column))
}

# Negate the magrittr %in% function
`%notin%` <- Negate(`%in%`)

# Main Script ----

## Import and merge metadata
jointMetaData <- read_csv(file = jointMetaData_file)
## add a cohort column to the jointMetaData table
jointMetaData <- jointMetaData %>%
  mutate(cohort = case_when(
    str_starts(samplename, "STP") ~ "STOP",
    str_starts(samplename, "LMH") ~ "Lyell",
    TRUE ~ NA_character_  # For any other prefix
  ))

## Gene ID versus GeneSymbol
gene_info <- read_csv(file = gene_info_file)
## Drop duplicated gene IDs
gene_info <- gene_info[!duplicated(gene_info$ensembl_gene_id), ]


## Import count file
SAGCQA0492_readCounts <- read.delim(file = SAGCQA0492_readCounts_file, sep = "\t")

# Check the basic structure
str(SAGCQA0492_readCounts)
# originalNames <- data.frame(colmnName = names(SAGCQA0492_readCounts))
# 
# subbedNames <- data.frame(subbedNames = gsub("^X\\.scratch\\.user\\.smit1924\\.annaWaverIronProject\\.aligned_data\\.|_S[0-9]+_GRCh38_Aligned\\.sortedByCoord\\.out\\.bam$", 
#                         "", 
#                         names(SAGCQA0492_readCounts)))
# 
# temp <- cbind(originalNames, subbedNames)
# 
# temp2 <- dplyr::left_join(temp, subbedNames, by = "subbedNames")
# 
# temp3 <- dplyr::left_join(temp2, jointMetaData, by = c("subbedNames" = "ULN"))
# 
# names <- data.frame(subbedNames = gsub("^X\\.scratch\\.user\\.smit1924\\.annaWaverIronProject\\.aligned_data\\.|_S[0-9]+_GRCh38_Aligned\\.sortedByCoord\\.out\\.bam$", 
#                                    "", 
#                                    names(SAGCQA0492_readCounts)))
# 
# temp4 <- dplyr::left_join(temp3, names, by = "subbedNames")

# clean the column names
names(SAGCQA0492_readCounts) <- gsub("^X\\.scratch\\.user\\.smit1924\\.annaWaverIronProject\\.aligned_data\\.|_S[0-9]+_GRCh38_Aligned\\.sortedByCoord\\.out\\.bam$", 
                              "", 
                              names(SAGCQA0492_readCounts))


## there are two rows in this count data that we don't need
## -"X.scratch.user.smit1924.annaWaverIronProject.aligned_data.SAGCFN_22_01736_S71_GRCh38_Aligned.sortedByCoord.out.bam" seems to be an additional sample that isn't ours
## - "Undetermined_GRCh38_Aligned.sortedByCoord.out.bam" are the transcripts not assigned to a sample and can be removed

# remove the columns we don't need
SAGCQA0492_readCounts <- dplyr::select(SAGCQA0492_readCounts, -SAGCFN_22_01736, -Undetermined_GRCh38_Aligned.sortedByCoord.out.bam)

# Create a mapping from ULN to samplename
uln_to_samplename <- setNames(jointMetaData$samplename, jointMetaData$ULN)

# Get current column names
current_colnames <- colnames(SAGCQA0492_readCounts)

# Replace column names where matches are found
# For non-matching columns (like possibly the first column with gene names/IDs), keep the original
new_colnames <- sapply(current_colnames, function(col) {
  if (col %in% names(uln_to_samplename)) {
    return(uln_to_samplename[col])
  } else {
    return(col)
  }
})

# Apply the new column names
colnames(SAGCQA0492_readCounts) <- new_colnames

# Plot the raw library sizes
## library size information for the paper
# calculate the median raw library size pre-deduplication
median_lib_size <- SAGCQA0492_readCounts %>%
  tibble::column_to_rownames("Geneid") %>%
  dplyr::select(all_of(jointMetaData$samplename)) %>%
  colSums() %>%
  median()


data.frame(colSums(SAGCQA0492_readCounts[, -1])) %>%
  rename(lib.size = colSums.SAGCQA0492_readCounts....1..) %>%
  tibble::rownames_to_column("samplename") %>%
  mutate(color = case_when(
    lib.size < 5000000 ~ "red",
    lib.size >= 5000000 & lib.size < 10000000 ~ "blue",
    lib.size >= 10000000 ~ "darkgreen"
  )) %>%
  ggplot(aes(x = reorder(samplename, -lib.size), y = lib.size, fill = color)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 5000000, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 10000000, linetype = "dashed", color = "darkgreen") +
  labs(y = "Library size (total number of mapped and quantified reads)",
       x = "Samplename") +
  scale_fill_identity() +
  coord_flip() +
  ggtitle("Library Size")
# ggsave(filename = output_librarySize_file,
#        width = 7,
#        height = 10,
#        units = "in",
#        dpi = 300)

## See if there is a visual correlation between RNA concentration and library size
data.frame(colSums(SAGCQA0492_readCounts[, -1])) %>%
  rename(lib.size = colSums.SAGCQA0492_readCounts....1..) %>%
  tibble::rownames_to_column("samplename") %>%
  dplyr::left_join(., jointMetaData, by = "samplename") %>%
  ggplot(aes(x = lib.size, y = Conc_ng_ul)) +
  geom_point()

## See if there is a visual correlation between RNA concentration and cluster detection
data.frame(colSums(SAGCQA0492_readCounts[, -1])) %>%
  rename(lib.size = colSums.SAGCQA0492_readCounts....1..) %>%
  tibble::rownames_to_column("samplename") %>%
  dplyr::left_join(., jointMetaData, by = "samplename") %>%
  ggplot(aes(x = lib.size, y = total_clusters_passing_filters_million)) +
  geom_point()

# Close Logging ----
# sink()

# Session Information ----
# Always good to capture for reproducibility
writeLines(capture.output(sessionInfo()), 
           file.path(output_dir, "session_info.txt"))