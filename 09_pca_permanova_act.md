# PCA and PERMANOVA to test bacterial community structure vs archaeal community types

The goal is to test if different bacterial communities can be linked to the archaeal community types.

```R
# load libraries
library(tidyverse)
library(phyloseq)
library(ggplot2)

# read taxonomy data
physeq <- readRDS("phyloseq.rds")

# read health data with read depth
load("health_data.RData")
read_depth <- health$covariates %>% mutate(ReadDepth = ReadDepthMillions*10^6) %>% select(Barcode, ReadDepth)

# melt the phyloseq object
melt_phyloseq <- psmelt(physeq) %>% filter(Abundance > 0) 

# list the people with Archaea 
archaea_samples <- melt_phyloseq %>% filter(Domain == "Archaea") %>% filter(Abundance > 0) %>% pull(Sample) %>% unique()

# filter the phyloseq object 
physeq <- prune_samples(archaea_samples, physeq)

# dominant archeaea for each individual
dominant <- subset_taxa(physeq, Domain == "Archaea") %>%
  psmelt() %>%
  filter(Abundance > 0) %>%
  group_by(Sample) %>%
  summarize(ACT = as.character(Species[which.max(Abundance)]))

# create the ACT groups
dominant$ACT <- ifelse(dominant$ACT %in% c("Methanobrevibacter_A smithii",  "Methanobrevibacter_A smithii_A",  "Methanobrevibacter_A sp900766745"), dominant$ACT, "Other")

# filter for bacteria and agglomerate to genus level
physeq_genus <- tax_glom(subset_taxa(physeq, Domain == "Bacteria"), taxrank = "Genus") 

# calculate prevalence
prevalence <- microbiome::prevalence(physeq_genus)

# filter for taxa with prevalence threshold as 5% of total samples
keepTaxa <- prevalence >= 0.05
physeq_genus <- prune_taxa(keepTaxa, physeq_genus) 

# calculate % of variance explained by each PC
x <- microbiome::transform(physeq_genus, "clr")
d <- t(microbiome::abundances(x)) 
pca <- prcomp(d) # we want samples as rows and taxa as columns for this function
eig <- (pca$sdev)^2
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig <- data.frame(eigenvalue = eig, variance = variance, cumvariance = cumvar)
pc1<-paste0(round(eig[1,2], 1), "%")
pc2<-paste0(round(eig[2,2], 1), "%")

# convert PCA scores to a data frame and add SampleID as a column
pca_scores <- as.data.frame(pca$x)
pca_scores$SampleID <- rownames(pca_scores)

# merge the PCA scores with the sample metadata and the ACT groups
pca_data <- left_join(pca_scores, dominant, by = c("SampleID" = "Sample"))

# create a plot 
p <- ggplot(pca_data, aes(x = PC1, y = PC2, fill = ACT, color = ACT)) +  
        geom_point(size=2, shape=21, color = "black") +
        scale_fill_manual(values = c("#7E6C9E", "#181A34", "#F4ED76", "#E5E5DF"), name = "Dominant Archaea") +
        scale_color_manual(values = c("#7E6C9E", "#181A34", "#F4ED76", "#E5E5DF"), name = "Dominant Archaea") +
        labs(x = pc1, y = pc2) + 
        theme_classic() +
        theme(legend.position = "none")

svg("PCA.svg", width = 3, height = 3)
print(p)
dev.off()

# PERMANOVA

# add dominant info to the distance data 
d <- d %>% as.data.frame()
d$SampleID <- rownames(d)
d <- d %>% left_join(pca_data%>%select(SampleID,ACT))

# add read depth data as ReadDepth_log10
read_depth <- read_depth %>%  mutate(ReadDepth_log10 = log10(ReadDepth)) %>% select(Barcode, ReadDepth_log10)
d <- d %>% left_join(read_depth, by = c("SampleID" = "Barcode"))

# Separete the metadata and d
metadata <- d %>% select(SampleID, ACT, ReadDepth_log10) %>% column_to_rownames("SampleID") %>% mutate(ACT = as.factor(ACT))
d <- d %>% select(-ACT, -ReadDepth_log10) %>% column_to_rownames("SampleID") %>% as.matrix()

# overall PERMANOVA

# Run adonis2 controlling for ReadDepth_log10
overall_perm <- vegan::adonis2(d ~ ACT + ReadDepth_log10, data = metadata, permutations=1000, method = "euclidean")  %>% as.data.frame()

# Extract R2 and P-value for ACT
ACT_stats <- data.frame(Factor = "ACT", R2 = overall_perm["ACT", "R2"], P = overall_perm["ACT", "Pr(>F)"], Comparision = "Overall")

# pairwise PERMANOVA

# calculate for each pairwise comparison
pairwise_perm <- data.frame()
for (comp in list(c("Methanobrevibacter_A smithii", "Methanobrevibacter_A smithii_A"), c("Methanobrevibacter_A smithii", "Methanobrevibacter_A sp900766745"), c("Methanobrevibacter_A smithii_A", "Methanobrevibacter_A sp900766745"))) {
# Print the comparison being made
cat("\nComparing:", comp[1], "vs", comp[2], "...\n")

# Subset metadata to only these two categories
keep_rows <- which(metadata$ACT %in% comp)
subset_meta <- metadata[keep_rows, ]

# Subset distance matrix accordingly
subset_d <- d[keep_rows, , drop = FALSE]

# Run adonis2 with same formula you used before
tmp <- vegan::adonis2(subset_d ~ ACT + ReadDepth_log10, data = subset_meta, permutations=1000, method = "euclidean")  %>% as.data.frame()

# Extract R2 and P-value for ACT
tmp <- data.frame(Factor = "ACT", R2 = tmp["ACT", "R2"], P = tmp["ACT", "Pr(>F)"], Comparision = paste(comp[1], "vs", comp[2]))

# Add to the pairwise_perm data frame
pairwise_perm <- rbind(pairwise_perm, tmp)
}

# Combine the overall and pairwise results
combined_results <- rbind(ACT_stats, pairwise_perm)

# Save the results so we can add the other cohorts and correct for multiple testing 
write.csv(combined_results, "PERMANOVA.csv", row.names = FALSE)
```

We afterwards combined the results from the other cohorts and corrected for multiple testing using FDR.