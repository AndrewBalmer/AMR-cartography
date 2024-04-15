
### Read in genetic data for S. pneumoniae
remove(list = ls())  # Remove all objects in the current workspace

# Install and load necessary packages
# install.packages("tidyverse")
# install.packages("curl")

# Load required packages
library(curl)
library(tidyverse)
library(seqinr)
library(ape)
library(ggdendro)

# Set working directory
setwd("/Users/ajb306/AMR-cartography/data/")

# Set DOI
doi <- "10.1186/s12864-017-4017-7"

# Construct the URL using the DOI
url <- paste0("https://doi.org/", doi)

# Print the URL
print(url)

# Use the curl_download function to download the file
file_urls <- c(
  "https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-017-4017-7/MediaObjects/12864_2017_4017_MOESM4_ESM.csv",
  "https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-017-4017-7/MediaObjects/12864_2017_4017_MOESM5_ESM.csv"
)

# Replace the file names as you want them to be saved locally
file_names <- c(
  "PBP_Sequence_dataset1.csv",
  "PBP_Sequence_dataset2.csv"
)

# Create a function to download files
download_files <- function(url, file_name) {
  curl_download(url, destfile = file_name, quiet = FALSE)
}

# Download files
for (i in seq_along(file_urls)) {
  download_files(file_urls[i], file_names[i])
}

# Read datasets into R
PBPseq_d1 <- read.csv("../data/PBP_Sequence_dataset1.csv", header=TRUE, sep=",", skip = 0)
PBPseq_d2 <- read.csv("../data/PBP_Sequence_dataset2.csv", header=TRUE, sep=",", skip = 0)

# Combine datasets
PBPseq <- rbind(PBPseq_d1, PBPseq_d2)
PBPseq_d1 <- NULL
PBPseq_d2 <- NULL

# Replace "TRUE" and "FALSE" with "T" and "F" respectively in the dataset
for (i in 1:ncol(PBPseq)) {
  PBPseq[, i] <- str_replace_all(PBPseq[, i], "TRUE", "T")
  PBPseq[, i] <- str_replace_all(PBPseq[, i], "FALSE", "F")
}

# Remove rows with missing values
PBPseq <- na.omit(PBPseq)
colnames(PBPseq)

# Read metadata table
tablemic_meta <- read.csv("../data/meta_data_Spneumoniae.csv", header=TRUE, sep=",", skip = 0)

# Perform a left join on LABID to merge the datasets
PBPseq <- left_join(tablemic_meta, PBPseq, by = c("LABID" = "LABID"))

# Extract a subset of the data for a test
date_test <- filter(PBPseq, LABID == 20154514 | LABID == 20153451)
date_test <- date_test[, 9:922]
date_test <- as.data.frame(t(date_test))

# Create a new column "Same" based on the equality of two columns
date_test$Same <- ifelse(date_test$V1 == date_test$V2, "Same", "Not-same")

# Create a lookup table for replacements
replacement_table <- data.frame(from = c("1930", "1938", "1940", "1941", "1959", "2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2013", "2016", "2018", "2021", "2027"),
                                to = c("30", "38", "40", "41", "59", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "13", "16", "18", "21", "27"))

# Loop through the replacement table and apply substitutions
for (i in seq(nrow(replacement_table))) {
  PBPseq$PT <- gsub(replacement_table$from[i], replacement_table$to[i], PBPseq$PT)
}

# Replace '/' with '-'
PBPseq$PT <- gsub("/", "-", PBPseq$PT)

# Save the dataset
save(PBPseq, file = "tablemic_pneumo_gen_3628.RData")

# Define samples to exclude - these samples have a clear indel mutation and so were removed from the genetic map
samples_to_exclude <- c("20156696", "20162849", "20151885", "20153985", "20154509", "2013224047", "2013218247", "2014200662", "5869-99", "2513-99")

# Exclude samples
PBPseq <- PBPseq %>%
  filter(!LABID %in% samples_to_exclude)

# Select relevant columns
PBPseq <- as.data.frame(PBPseq[, 9:922])

# Calculate genetic distances
pbp_dist <- as.data.frame(as.matrix(dist.gene(PBPseq, method = "pairwise", pairwise.deletion = F, variance = F)))

# also generate a small version just for illustration purposes
pbp_dist_200 <- as.data.frame(as.matrix(dist.gene(slice_sample(PBPseq, n = 200), method = "pairwise", pairwise.deletion = F, variance = F)))

# Save LABID information
LABID_for3628_isolates <- PBPseq$LABID
save(LABID_for3628_isolates, file = "LABID_for3628_isolates.RData")

# Save datasets
save(PBPseq, file = "tablemic_pneumo_gen_3628.RData")
save(pbp_dist, file = "tablemic_pneumo_3628_meta_gen_distance_matrix.RData")
save(pbp_dist_200, file = "tablemic_pneumo_3628_meta_gen_distance_matrix_200_subsample.RData")
