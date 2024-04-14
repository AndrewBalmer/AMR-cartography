
# Clear the environment
remove(list = ls())

# Install and load necessary packages
#install.packages("tidyverse")
#install.packages("curl")

# Load required packages
library(tidyverse)
library(curl)

# Set working directory
setwd("/Users/ajb306/AMR-cartography/data")

# Replace "YOUR_DOI" with the actual DOI
doi <- "10.1186/s12864-017-4017-7"

# Construct the URL using the DOI (this is an example, actual URLs may vary)
url <- paste0("https://doi.org/", doi)

# Print the URL (optional)
print(url)

# Use the curl_download function to download the file
file_urls <- c(
  "https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-017-4017-7/MediaObjects/12864_2017_4017_MOESM1_ESM.csv",
  "https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-017-4017-7/MediaObjects/12864_2017_4017_MOESM2_ESM.csv"
)

# Replace the file names as you want them to be saved locally
file_names <- c(
  "MIC_pneumo.csv",
  "MIC_pneumo2.csv"
)

# Create a function to download files
download_files <- function(url, file_name) {
  curl_download(url, destfile = file_name, quiet = FALSE)
}

# Download files
for (i in seq_along(file_urls)) {
  download_files(file_urls[i], file_names[i])
}

# Read MIC data from two CSV files
tablemic <- read.csv("../data/MIC_pneumo.csv", header=TRUE, sep=",", skip = 0)
tablemic2 <- read.csv("../data/MIC_pneumo2.csv", header=TRUE, sep=",", skip = 0)

# Subset the columns with MIC values
tablemic <- tablemic[,1:8]
tablemic2 <- tablemic2[,1:8]

# Combine the two datasets
tablemic <- rbind(tablemic, tablemic2)
tablemic2 <- NULL

# Group by 'PT' and add counts
tablemic <- tablemic %>% group_by(PT) %>% add_count()

colnames(tablemic) <- c("LABID","PT","Penicillin","Amoxicillin","Meropenem","Cefotaxime","Ceftriaxone","Cefuroxime")

# Select relevant columns and remove rows with NA values
tablemic <- tablemic[!is.na(tablemic$Penicillin), ]
tablemic <- tablemic[!is.na(tablemic$Amoxicillin), ]
tablemic <- tablemic[!is.na(tablemic$Meropenem), ]
tablemic <- tablemic[!is.na(tablemic$Cefotaxime), ]
tablemic <- tablemic[!is.na(tablemic$Ceftriaxone), ]
tablemic <- tablemic[!is.na(tablemic$Cefuroxime), ]

# Extract metadata and MIC values
tablemic_meta <- as.data.frame(tablemic[,1:8]) 
tablemic <- as.data.frame(tablemic[,3:8])

# Convert MIC values to numeric and perform log2 transformation
for (i in 1:ncol(tablemic)) {
  tablemic[,i] <- log2(as.numeric(tablemic[,i]))
  tablemic[,i] <- round(tablemic[,i])
}

# Some of the PBP-types had been converted into dates by excel and need to be replaced
# Create a lookup table for replacements
replacement_table <- data.frame(from = c("1930", "1938", "1940", "1941", "1959", "2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2013", "2016", "2018", "2021", "2027"),
                                to = c("30", "38", "40", "41", "59", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "13", "16", "18", "21", "27"))

# Loop through the replacement table and apply substitutions
for (i in seq(nrow(replacement_table))) {
  tablemic_meta$PT <- gsub(replacement_table$from[i], replacement_table$to[i], tablemic_meta$PT)
}

# Replace '/' with '-'
tablemic_meta$PT <- gsub("/", "-", tablemic_meta$PT)

# Save the dataframe as a CSV file
write.csv(tablemic, file = "/Users/ajb306/AMR-cartography/data/MIC_table_Spneumoniae.csv", row.names = FALSE)
write.csv(tablemic_meta, file = "/Users/ajb306/AMR-cartography/data/meta_data_Spneumoniae.csv", row.names = FALSE)

# Generate distance matrix
dist_pne <- dist(slice_sample(tablemic)) # Use this analysis if you would like to replicate the full analysis (be aware this may take a long time to run).
#dist_pne <- dist(slice_sample(tablemic, n = 200))

# Save file
save(dist_pne, file="/Users/ajb306/AMR-cartography/data/phenotype_distance_matrix_Spneumoniae.Rdata")

