# Check abundance from Halobacteriota

```R
# load libraries
library(tidyverse)
library(phyloseq)
library(microbiome)

# read new taxonomy data
physeq_list <- readRDS("physeq_list_metagenomes_with_metadata_and_single_sample_per_individual.rds")

names(physeq_list)

# List the phyloseq objects to be plotted
studies <- c("AsnicarF_2021-1088samples", "JieZ_2017-385samples", "KarlssonFH_2013-145samples", "LiJ_2017-196samples", 
             "MetaCardis_2020_a-1726samples", "QinN_2014-237samples", "SchirmerM_2016-471samples",
             "YachidaS_2019-609samples", "ShaoY_2019-753samples", "MehtaRS_2018-308samples",
             "VatanenT_2016-212samples", "BackhedF_2015-200samples", "HMP_2019_ibdmdb-130samples", "BrooksB_2017-30samples")

##############################################################################
# 1. Filter for samples with Halophilic archaea, keep bacteria as is, calculate relative abundance 
################################################################################

# table to store the abundance tables for each study
abundance_per_study <- data.frame()

for(study_name in studies) {

  # Print the study name for CMD and create the taxa_counts
  cat("Processing study:", study_name, "\n")

  # Filter for the study being analyzed
  phyloseq_obj <- physeq_list[[study_name]]

  # Melt the phyloseq object and remove rows with zero abundance
  melt_phyloseq <- psmelt(phyloseq_obj) %>% filter(Abundance > 0) 

  # List the samples with Halophilic Archaea
  archaea_samples <- melt_phyloseq %>% filter(Phylum == "Halobacteriota") %>% pull(SampleID) %>% unique()

  # Cat the number of samples with Halophilic Archaea
  cat("Number of samples with Halophilic Archaea:", length(archaea_samples), "\n")

  # Keep individuals with Halophilic Archaea detected from the phyloseq object and calculate relative abundance with microbiome package
  phyloseq_obj <- prune_samples(archaea_samples, phyloseq_obj)
  phyloseq_obj <- microbiome::transform(phyloseq_obj, "compositional")
  melt_phyloseq <- psmelt(phyloseq_obj) %>% mutate(Abundance = Abundance * 100) # Convert to percentage

  # Extract abundance of Halophilic Archaea for each sample
  archaea_abundance <- melt_phyloseq %>% filter(Phylum == "Halobacteriota") %>% filter(Abundance > 0) 

  # Add to the abundance table for the study
  abundance_per_study <- rbind(abundance_per_study, data.frame(Study = study_name, SampleID = archaea_abundance$SampleID, Abundance = archaea_abundance$Abundance, Species = archaea_abundance$Species))
}

# Calculate summary statistics for the abundance of Halophilic Archaea across all studies
summary_stats <- abundance_per_study %>% group_by(Study) %>%
  summarise(Mean_Abundance = mean(Abundance), Median_Abundance = median(Abundance), Min_Abundance = min(Abundance), Max_Abundance = max(Abundance), SD_Abundance = sd(Abundance), Count_Samples = n())

# Save the abundance table and summary statistics to CSV files
write.csv(summary_stats, "halophilic_archaea_abundance_summary_stats.csv", row.names = FALSE)

# Now, do the summary statistics at the species level
species_summary_stats <- abundance_per_study %>% group_by(Study, Species) %>%
  summarise(Mean_Abundance = mean(Abundance), Median_Abundance = median(Abundance), Min_Abundance = min(Abundance), Max_Abundance = max(Abundance), SD_Abundance = sd(Abundance), Count_Samples = n())    

# Save the species-level summary statistics to a CSV file
write.csv(species_summary_stats, "halophilic_archaea_abundance_species_summary_stats.csv", row.names = FALSE)
```