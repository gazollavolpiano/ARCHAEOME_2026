# Archaea tree highlighting the major species found in the human gut

We will use the GTDB r214 archaeal tree and metadata to create a tree highlighting the major archaeal species found in the human gut.

```bash
# download the tree from GTDB r214
wget https://data.gtdb.ecogenomic.org/releases/release214/214.0/ar53_r214.tree

# download the GTDB r214 metadata for archaea and decompress it:
wget https://data.gtdb.ecogenomic.org/releases/release214/214.0/ar53_metadata_r214.tar.gz 
tar -xvzf ar53_metadata_r214.tar.gz
rm ar53_metadata_r214.tar.gz
```

Use R to plot the tree with the presence of the species in human gut samples. 

```R
# library
library(tidyverse)

# read tree and inspect the tip labels
tree <- ape::read.tree("ar53_r214.tree")
tree$tip.label [1:2] # "GB_GCA_002502135.1" "GB_GCA_002495465.1"
length(tree$tip.label) # 4416 tips

# read the metadata file
metadata <- read.csv("ar53_metadata_r214.tsv", sep = "\t", stringsAsFactors = FALSE) %>%
  select(accession, gtdb_taxonomy) %>%
  filter(accession %in% tree$tip.label)

# add species column from gtdb_taxonomy 
metadata$species <- gsub(".*;s__", "", metadata$gtdb_taxonomy)

# list of major archaeal species found in the human gut
major_taxa <- c("Methanobrevibacter_A smithii", "Methanobrevibacter_A smithii_A", "Methanobrevibacter_A sp900766745", "Methanobrevibacter_A woesei", 
                "Methanosphaera stadtmanae", "Methanosphaera sp900322125", "Methanomassiliicoccus_A intestinalis", "Methanocorpusculum sp001940805", 
                "UBA71 sp006954425", "UBA71 sp905187815", "UBA71 sp944320485")

# add column to the metadata indicating if the species is present in the human gut
metadata$in_gut <- ifelse(metadata$species %in% major_taxa, "Major species found in human gut", "Other species")

# add the phylum column to the metadata
metadata$phylum <- gsub(".*;p__", "", metadata$gtdb_taxonomy)
metadata$phylum <- gsub(";.*", "", metadata$phylum)

# save metadata to plot the tree with microreact
write.csv(metadata, "ar53_r214_metadata_for_Microreact.csv", row.names = FALSE)
```

Most of the work is done in Microreact, where we can upload the tree and the metadata to create an interactive tree highlighting the major archaeal species found in the human gut.