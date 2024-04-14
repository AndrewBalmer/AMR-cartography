
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

# Extract metadata and MIC values
tablemic_meta <- as.data.frame(tablemic[,1:2])
tablemic <- as.data.frame(tablemic[,3:8])

# Rename columns
colnames(tablemic) <- c("Penicillin","Amoxicillin","Meropenem","Cefotaxime","Ceftriaxone","Cefuroxime")

# Convert MIC values to numeric and perform log2 transformation
for (i in 1:ncol(tablemic)) {
  tablemic[,i] <- log2(as.numeric(tablemic[,i]))
  tablemic[,i] <- round(tablemic[,i])
}

# Read additional MIC metadata
tablemicpne <- read.csv("../data/MIC_pneumo.csv", header=TRUE, sep=",", skip = 0)
tablemicpne_meta_data <- read.csv("../data/MIC_S.Pneumo_metadata.csv", header=TRUE, sep=",", skip = 0)

# Separate and organize metadata columns
tablemicpne_meta_data <- separate(tablemicpne_meta_data, PEN, c("MIC1", "PEN_Sign", "PEN"), " ")
tablemicpne_meta_data <- separate(tablemicpne_meta_data, AMO, c("MIC2", "AMO_Sign", "AMO"), " ")
tablemicpne_meta_data <- separate(tablemicpne_meta_data, MER, c("MIC3", "MER_Sign", "MER"), " ")
tablemicpne_meta_data <- separate(tablemicpne_meta_data, TAX, c("MIC4", "TAX_Sign", "TAX"), " ")
tablemicpne_meta_data <- separate(tablemicpne_meta_data, CFT, c("MIC5", "CFT_Sign", "CFT"), " ")
tablemicpne_meta_data <- separate(tablemicpne_meta_data, CFX, c("MIC6", "CFX_Sign", "CFX"), " ")

# Select relevant columns and create metadata and MIC datasets
tablemicpne_meta_data <- select(tablemicpne_meta_data, -contains("MIC"))
tablemicpne_meta_data <- select(tablemicpne_meta_data, -contains("Sign"))
tablemicpne_meta <- as.data.frame(tablemicpne_meta_data[, 1:7])
tablemicpne <- as.data.frame(tablemicpne_meta_data[, 8:13])


# Rename columns
colnames(tablemicpne) <- c("Penicillin","Amoxicillin","Meropenem","Cefotaxime","Ceftriaxone","Cefuroxime")

# Combine metadata and MIC datasets
tablemicpne <- cbind(tablemicpne_meta, tablemicpne)
tablemicpne <- na.omit(tablemicpne)
tablemicpne_meta <- tablemicpne[, 1:7]
tablemicpne <- tablemicpne[, 8:13]

# Join metadata with the initial dataset
tablemic_metadata <- left_join(tablemic_meta, tablemicpne_meta, by = "LABID") 

# Select relevant columns and remove rows with NA values
tablemic <- cbind(tablemic_metadata, tablemic)
tablemic <- tablemic[!is.na(tablemic$Penicillin), ]
tablemic <- tablemic[!is.na(tablemic$Amoxicillin), ]
tablemic <- tablemic[!is.na(tablemic$Meropenem), ]
tablemic <- tablemic[!is.na(tablemic$Cefotaxime), ]
tablemic <- tablemic[!is.na(tablemic$Ceftriaxone), ]
tablemic <- tablemic[!is.na(tablemic$Cefuroxime), ]
tablemic_meta <- tablemic
tablemic <- tablemic[,9:14]

# Save the dataframe as a CSV file
write.csv(tablemic, file = "/Users/ajb306/AMR-cartography/data/MIC_table_Spneumoniae.csv", row.names = FALSE)
write.csv(tablemic_meta, file = "/Users/ajb306/AMR-cartography/data/meta_data_Spneumoniae.csv", row.names = FALSE)


