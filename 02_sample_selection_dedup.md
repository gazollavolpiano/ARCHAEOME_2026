# Sample selection and metadata addition for the phyloseq objects

We need to select a single representative sample per individual (sample with the highest total read count) and add metadata to the phyloseq objects.

The file `physeq_list_11494_metagenomes.rds` is a list that contains the phyloseq objects with the taxonomy data, with multiple samples per individual for some studies.

```R
# load libraries
library(tidyverse)
library(phyloseq)
library(microbiome)

# read new taxonomy data
physeq_list <- readRDS("physeq_list_11494_metagenomes.rds")

#######################################################################################################################
#      For the studies with multiple samples per individual, select the sample with the highest number of reads
#######################################################################################################################

##################
# ShaoY_2019
##################

# Read metadata for the study ShaoY_2019
metadata <- readRDS("metadata_1644_samples_753_individuals.rds") 

# Read read counts and tidy it
read_numbers <- read.table("fastq_stats_1644_samples.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE) %>% select(V1, V4) 
read_numbers$V1 <- basename(read_numbers$V1) %>% gsub("_..fastq", "", .)
read_numbers <- read_numbers %>% group_by(V1) %>% summarise(V4 = mean(V4)) 

setdiff(metadata$NCBI_accession, read_numbers$V1) # no missing samples

# Merge the read count data into the metadata by matching the sample IDs (NCBI_accession)
metadata <- metadata %>% left_join(read_numbers, by = c("NCBI_accession" = "V1")) 
head(metadata) # Looks good!
metadata$V4    # Looks good!

# Get the list of sample names from the phyloseq object for this study
samples_in_physeq <- sample_names(physeq_list[["ShaoY_2019-1644samples"]]) 

# Filter the metadata to retain only the samples that are present in the phyloseq object
metadata <- metadata %>% filter(NCBI_accession %in% samples_in_physeq) 
dim(metadata) # 1644 samples x 20 variables, as expected!

# For each individual (subject_id), select the sample with the highest number of reads
samples_to_keep <- metadata %>% 
  mutate(number_reads = as.numeric(V4)) %>%
  group_by(subject_id) %>% 
  filter(number_reads == max(number_reads)) %>% 
  pull(NCBI_accession)

length(samples_to_keep) # Should be 753 samples (one per individual)

# Set up the sample IDs in the metadata for consistency with the phyloseq object
metadata$SampleID <- metadata$NCBI_accession
rownames(metadata) <- metadata$SampleID

# Update the sample_data in the phyloseq object with the modified metadata
sample_data(physeq_list[["ShaoY_2019-1644samples"]]) <- sample_data(metadata)

# Create a new phyloseq object with only the selected samples (one per individual)
physeq_list[["ShaoY_2019-753samples"]] <- subset_samples(physeq_list[["ShaoY_2019-1644samples"]], SampleID %in% samples_to_keep)
physeq_list[["ShaoY_2019-753samples"]]
# otu_table()   OTU Table:         [ 9080 taxa and 753 samples ]
# sample_data() Sample Data:       [ 753 samples by 21 sample variables ]
# tax_table()   Taxonomy Table:    [ 9080 taxa by 7 taxonomic ranks ]

##################
# MehtaRS_2018
##################

# Read metadata for the study MehtaRS_2018
metadata <- readRDS("metadata_928_samples_308_individuals.rds") 

# Read read counts and tidy it
read_numbers <- read.table("fastq_stats_928_samples.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE) %>% select(V1, V4) 
read_numbers$V1 <- basename(read_numbers$V1) %>% gsub("_..fastq", "", .)
read_numbers <- read_numbers %>% group_by(V1) %>% summarise(V4 = mean(V4)) 

setdiff(metadata$NCBI_accession, read_numbers$V1) # no missing samples

# Merge the read count data into the metadata by matching the sample IDs (NCBI_accession)
metadata <- metadata %>% left_join(read_numbers, by = c("NCBI_accession" = "V1")) 
head(metadata) # Looks good!
metadata$V4    # Looks good!

# Get the list of sample names from the phyloseq object for this study
samples_in_physeq <- sample_names(physeq_list[["MehtaRS_2018-928samples"]])

# Filter the metadata to retain only the samples that are present in the phyloseq object
metadata <- metadata %>% filter(NCBI_accession %in% samples_in_physeq) 
dim(metadata) # 928 samples x 7 variables, as expected!

# For each individual (subject_id), select the sample with the highest number of reads
samples_to_keep <- metadata %>% 
  mutate(number_reads = as.numeric(V4)) %>%
  group_by(subject_id) %>% 
  filter(number_reads == max(number_reads)) %>% 
  pull(NCBI_accession)

length(samples_to_keep) # Should be 308 samples (one per individual)

# Set up the sample IDs in the metadata for consistency with the phyloseq object
metadata$SampleID <- metadata$NCBI_accession
rownames(metadata) <- metadata$SampleID

# Update the sample_data in the phyloseq object with the modified metadata
sample_data(physeq_list[["MehtaRS_2018-928samples"]]) <- sample_data(metadata)

# Create a new phyloseq object with only the selected samples (one per individual)
physeq_list[["MehtaRS_2018-308samples"]] <- subset_samples(physeq_list[["MehtaRS_2018-928samples"]], SampleID %in% samples_to_keep)
physeq_list[["MehtaRS_2018-308samples"]]
# otu_table()   OTU Table:         [ 8644 taxa and 308 samples ]
# sample_data() Sample Data:       [ 308 samples by 8 sample variables ]
# tax_table()   Taxonomy Table:    [ 8644 taxa by 7 taxonomic ranks ]

##################
# VatanenT_2016
##################
  
# Read metadata for the study VatanenT_2016
metadata <- readRDS("metadata_785_samples_212_individuals.rds") 

# Read read counts and tidy it
read_numbers <- read.table("fastq_stats_785_samples.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE) %>% select(V1, V4) 
read_numbers$V1 <- basename(read_numbers$V1) %>% gsub("_..fastq", "", .)
read_numbers <- read_numbers %>% group_by(V1) %>% summarise(V4 = mean(V4)) 

setdiff(metadata$NCBI_accession, read_numbers$V1) # no missing samples

# Merge the read count data into the metadata by matching the sample IDs (NCBI_accession)
metadata <- metadata %>% left_join(read_numbers, by = c("NCBI_accession" = "V1")) 
head(metadata) # Looks good!
metadata$V4    # Looks good!

# Get the list of sample names from the phyloseq object for this study
samples_in_physeq <- sample_names(physeq_list[["VatanenT_2016-785samples"]])

# Filter the metadata to retain only the samples that are present in the phyloseq object
metadata <- metadata %>% filter(NCBI_accession %in% samples_in_physeq) 
dim(metadata) # 785 samples x 22 variables, as expected!

# For each individual (subject_id), select the sample with the highest number of reads
samples_to_keep <- metadata %>% 
  mutate(number_reads = as.numeric(V4)) %>%
  group_by(subject_id) %>% 
  filter(number_reads == max(number_reads)) %>% 
  pull(NCBI_accession)

length(samples_to_keep) # Should be 212 samples (one per individual)

# Set up the sample IDs in the metadata for consistency with the phyloseq object
metadata$SampleID <- metadata$NCBI_accession
rownames(metadata) <- metadata$SampleID

# Update the sample_data in the phyloseq object with the modified metadata
sample_data(physeq_list[["VatanenT_2016-785samples"]]) <- sample_data(metadata)

# Create a new phyloseq object with only the selected samples (one per individual)
physeq_list[["VatanenT_2016-212samples"]] <- subset_samples(physeq_list[["VatanenT_2016-785samples"]], SampleID %in% samples_to_keep)
physeq_list[["VatanenT_2016-212samples"]]
# otu_table()   OTU Table:         [ 7800 taxa and 212 samples ]
# sample_data() Sample Data:       [ 212 samples by 23 sample variables ]
# tax_table()   Taxonomy Table:    [ 7800 taxa by 7 taxonomic ranks ]

##################
# BackhedF_2015
##################
  
# Read metadata for the study BackhedF_2015
metadata <- readRDS("metadata_400_samples_200_individuals.rds") 

# Read read counts and tidy it
read_numbers <- read.table("fastq_stats_400_samples.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE) %>% select(V1, V4) 
read_numbers$V1 <- basename(read_numbers$V1) %>% gsub("_..fastq", "", .)
read_numbers <- read_numbers %>% group_by(V1) %>% summarise(V4 = mean(V4)) 

setdiff(metadata$NCBI_accession, read_numbers$V1) # no missing samples

# Merge the read count data into the metadata by matching the sample IDs (NCBI_accession)
metadata <- metadata %>% left_join(read_numbers, by = c("NCBI_accession" = "V1")) 
head(metadata) # Looks good!
metadata$V4    # Looks good!

# Get the list of sample names from the phyloseq object for this study
samples_in_physeq <- sample_names(physeq_list[["BackhedF_2015-400samples"]])

# Filter the metadata to retain only the samples that are present in the phyloseq object
metadata <- metadata %>% filter(NCBI_accession %in% samples_in_physeq) 
dim(metadata) # 400 samples x 18 variables, as expected!

# For each individual (subject_id), select the sample with the highest number of reads
samples_to_keep <- metadata %>% 
  mutate(number_reads = as.numeric(V4)) %>%
  group_by(subject_id) %>% 
  filter(number_reads == max(number_reads)) %>% 
  pull(NCBI_accession)

length(samples_to_keep) # Should be 200 samples (one per individual)

# Set up the sample IDs in the metadata for consistency with the phyloseq object
metadata$SampleID <- metadata$NCBI_accession
rownames(metadata) <- metadata$SampleID

# Update the sample_data in the phyloseq object with the modified metadata
sample_data(physeq_list[["BackhedF_2015-400samples"]]) <- sample_data(metadata)

# Create a new phyloseq object with only the selected samples (one per individual)
physeq_list[["BackhedF_2015-200samples"]] <- subset_samples(physeq_list[["BackhedF_2015-400samples"]], SampleID %in% samples_to_keep)
physeq_list[["BackhedF_2015-200samples"]]
# otu_table()   OTU Table:         [ 7853 taxa and 200 samples ]
# sample_data() Sample Data:       [ 200 samples by 19 sample variables ]
# tax_table()   Taxonomy Table:    [ 7853 taxa by 7 taxonomic ranks ]

##################
# HMP_2019_ibdmdb
##################
  
# Read metadata for the study HMP_2019_ibdmdb
metadata <- readRDS("metadata_1585_samples_130_individuals.rds") 

# Read read counts and tidy it
read_numbers <- read.table("fastq_stats_1584_samples.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE) %>% select(V1, V4) 
read_numbers$V1 <- basename(read_numbers$V1) %>% gsub("_..fastq", "", .)
read_numbers <- read_numbers %>% group_by(V1) %>% summarise(V4 = mean(V4)) 

setdiff(metadata$sample_id, read_numbers$V1) # no missing samples, this time we use sample_id instead of NCBI_accession
# "HSM5MD8B" is missing... that's ok!

# Merge the read count data into the metadata by matching the sample IDs (sample_id)
metadata <- metadata %>% left_join(read_numbers, by = c("sample_id" = "V1")) 
head(metadata) # Looks good!
metadata$V4    # Looks good!

# Get the list of sample names from the phyloseq object for this study
samples_in_physeq <- sample_names(physeq_list[["HMP_2019_ibdmdb-1572samples"]])

# Filter the metadata to retain only the samples that are present in the phyloseq object
metadata <- metadata %>% filter(sample_id %in% samples_in_physeq) 
dim(metadata) # 1572 samples x 14 variables, as expected!

# For each individual (subject_id), select the sample with the highest number of reads
samples_to_keep <- metadata %>% 
  mutate(number_reads = as.numeric(V4)) %>%
  group_by(subject_id) %>% 
  filter(number_reads == max(number_reads)) %>% 
  pull(sample_id)

length(samples_to_keep) # Should be 130 samples (one per individual)

# Set up the sample IDs in the metadata for consistency with the phyloseq object
metadata$SampleID <- metadata$sample_id
rownames(metadata) <- metadata$SampleID

# Update the sample_data in the phyloseq object with the modified metadata
sample_data(physeq_list[["HMP_2019_ibdmdb-1572samples"]]) <- sample_data(metadata)

# Create a new phyloseq object with only the selected samples (one per individual)
physeq_list[["HMP_2019_ibdmdb-130samples"]] <- subset_samples(physeq_list[["HMP_2019_ibdmdb-1572samples"]], SampleID %in% samples_to_keep)
physeq_list[["HMP_2019_ibdmdb-130samples"]]
# otu_table()   OTU Table:         [ 8356 taxa and 130 samples ]
# sample_data() Sample Data:       [ 130 samples by 15 sample variables ]
# tax_table()   Taxonomy Table:    [ 8356 taxa by 7 taxonomic ranks ]

##################
# BrooksB_2017
##################
  
# Read metadata for the study BrooksB_2017
metadata <- readRDS("metadata_408_samples_30_individuals.rds") 

# Read read counts and tidy it
read_numbers <- read.table("fastq_stats_408_samples.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE) %>% select(V1, V4) 
read_numbers$V1 <- basename(read_numbers$V1) %>% gsub("_..fastq", "", .)
read_numbers <- read_numbers %>% group_by(V1) %>% summarise(V4 = mean(V4)) 

setdiff(metadata$NCBI_accession, read_numbers$V1) # no missing samples

# Merge the read count data into the metadata by matching the sample IDs (NCBI_accession)
metadata <- metadata %>% left_join(read_numbers, by = c("NCBI_accession" = "V1")) 
head(metadata) # Looks good!
metadata$V4    # Looks good!

# Get the list of sample names from the phyloseq object for this study
samples_in_physeq <- sample_names(physeq_list[["BrooksB_2017-408samples"]])

# Filter the metadata to retain only the samples that are present in the phyloseq object
metadata <- metadata %>% filter(NCBI_accession %in% samples_in_physeq) 
dim(metadata) # 408 samples x 18 variables, as expected!

# For each individual (subject_id), select the sample with the highest number of reads
samples_to_keep <- metadata %>% 
  mutate(number_reads = as.numeric(V4)) %>%
  group_by(subject_id) %>% 
  filter(number_reads == max(number_reads)) %>% 
  pull(NCBI_accession)

length(samples_to_keep) # Should be 30 samples (one per individual)

# Set up the sample IDs in the metadata for consistency with the phyloseq object
metadata$SampleID <- metadata$NCBI_accession
rownames(metadata) <- metadata$SampleID

# Update the sample_data in the phyloseq object with the modified metadata
sample_data(physeq_list[["BrooksB_2017-408samples"]]) <- sample_data(metadata)

# Create a new phyloseq object with only the selected samples (one per individual)
physeq_list[["BrooksB_2017-30samples"]] <- subset_samples(physeq_list[["BrooksB_2017-408samples"]], SampleID %in% samples_to_keep)
physeq_list[["BrooksB_2017-30samples"]]
# otu_table()   OTU Table:         [ 2671 taxa and 30 samples ]
# sample_data() Sample Data:       [ 30 samples by 19 sample variables ]
# tax_table()   Taxonomy Table:    [ 2671 taxa by 7 taxonomic ranks ]

#######################################################################################################################
#      For the studies with a single sample per individual, add the metadata to the phyloseq object
#######################################################################################################################

for(study_name in c("AsnicarF_2021-1088samples", "SchirmerM_2016-471samples", "JieZ_2017-385samples", 
                    "QinN_2014-237samples", "LiJ_2017-196samples", "KarlssonFH_2013-145samples", 
                    "MetaCardis_2020_a-1726samples", "YachidaS_2019-609samples")){

  # Remove the number of samples from the study name
  study_name_simple <- strsplit(study_name, "-")[[1]][1]
  # Read metadata 
  metadata_path <- list.files(paste0("/labs/mericlab/curatedMetagenomicData/",study_name_simple), pattern = "rds", full.names = TRUE)
  metadata <- readRDS(metadata_path)
  # Filter to retain only the samples that are present in the phyloseq object
  samples_in_physeq <- sample_names(physeq_list[[study_name]])

  if(study_name_simple == "QinN_2014"){
    metadata$sample_id <- gsub("-", ".", metadata$sample_id) # harmonize sample_id format
  }

  if(study_name_simple %in% c("YachidaS_2019", "SchirmerM_2016", "QinN_2014", "MetaCardis_2020_a")){
    metadata <- metadata %>% filter(sample_id %in% samples_in_physeq)
    metadata$SampleID <- metadata$sample_id
    rownames(metadata) <- metadata$SampleID   
  } else {
    metadata <- metadata %>% filter(NCBI_accession %in% samples_in_physeq)
    metadata$SampleID <- metadata$NCBI_accession
    rownames(metadata) <- metadata$SampleID
  }
  cat(study_name, " - size of metadata is:", dim(metadata), "\n")
  # Add metadata to the phyloseq object
  sample_data(physeq_list[[study_name]]) <- sample_data(metadata)
}

# Save the updated phyloseq objects
saveRDS(physeq_list, "physeq_list_metagenomes_with_metadata_and_single_sample_per_individual.rds")

```
