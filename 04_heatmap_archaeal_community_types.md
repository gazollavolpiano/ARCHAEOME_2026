# Heatmap of archaeal community types 

This is an example of the heatmap created for FINRISK 2002 for participants with detectable archaea (n = 3,156). 

Samples are grouped by their dominant archaeal species, highlighting distinct archaeal community types. The same processing steps were applied as for the public cohorts.

We will also export a file with the abundance table so that it can be shared in the manuscript

```R
# load libraries
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ComplexHeatmap)
library(circlize)  
library(grid) 

############################################################
# 1. filter for samples with archaeal detections
#
# 2. organize taxa as follows:  
#    a. retain species that are part of a predefined list of major taxa (prevalence ≥0.1% in at least 10 of the 15 metagenomic studies)
#    b. for other archaeal species:
#       - aggregate the counts at the genus level
#       - for Methanobrevibacter_A, Methanosphaera, Methanomassiliicoccus_A, Methanocorpusculum, and UBA71, keep the aggregated genus-level data
#       - for all other genera, combine them into a single group 
#################################################################################

# read taxonomy data
physeq <- readRDS("phyloseq.rds")

# list major (prevalent) taxa to be kept at species level
major_taxa <- c("Methanobrevibacter_A smithii", "Methanobrevibacter_A smithii_A", "Methanobrevibacter_A sp900766745", "Methanobrevibacter_A woesei", 
                "Methanosphaera stadtmanae", "Methanosphaera sp900322125", "Methanomassiliicoccus_A intestinalis", "Methanocorpusculum sp001940805", 
                "UBA71 sp006954425", "UBA71 sp905187815", "UBA71 sp944320485")

# genera to be grouped
major_genera <- c("Methanobrevibacter_A", "Methanosphaera", "Methanocorpusculum", "Methanomassiliicoccus_A", "UBA71")

# melt the phyloseq object and remove rows with zero abundance
melt_phyloseq <- psmelt(physeq) %>% filter(Abundance > 0) 

# list the samples with Archaea
archaea_samples <- melt_phyloseq %>% filter(Domain == "Archaea") %>% pull(Sample) %>% unique()

# keep individuals with Archaea detected
melt_phyloseq <- melt_phyloseq %>% filter(Sample %in% archaea_samples)

# for Archaea:
# 1. keep rows with major_taxa
archaea_major <- melt_phyloseq %>% 
                  filter(Domain == "Archaea", Species %in% major_taxa) %>%
                  mutate(Domain = "Archaea", Taxa = Species) %>%
                  select(Sample, Domain, Taxa, Abundance) 

# 2. for the other rows, aggregate at the genus level (per sample)
archaea_other <- melt_phyloseq %>% 
                filter(Domain == "Archaea", !(Species %in% major_taxa)) %>%
                group_by(Sample, Genus) %>%
                summarise(Abundance = sum(Abundance), .groups = "drop") %>%
                mutate(Domain = "Archaea", Taxa = Genus) %>%
                select(Sample, Domain, Taxa, Abundance)

# 3. among aggregated rows, keep those genera that are in major_genera
archaea_other_major <- archaea_other %>% filter(Taxa %in% major_genera)

# 4. for genera not in major_genera, aggregate them into one group called "Rare Archaea"
archaea_other_rare <- archaea_other %>%
  filter(!(Taxa %in% major_genera)) %>%
  group_by(Sample) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  mutate(Taxa = "Rare Archaea", Domain = "Archaea")

# combine all the Archaea parts
archaea_final <- bind_rows(
  archaea_major,       # species-level info for major taxa (untouched)
  archaea_other_major, # aggregated at genus level for major genera
  archaea_other_rare   # aggregated as "Rare" for all other genera
)

# save the archaea_final abundance table, so it can be shared
# anonymize Samples by replacing them with generic IDs (e.g. Sample1, Sample2, etc.)
archaea_final %>%
  mutate(Sample = paste0("Sample", as.integer(factor(Sample)))) %>%
  arrange(desc(Abundance)) %>%
  write_csv("archaea_abundance_table.csv")

# add bacteria component
bacteria_component <- melt_phyloseq %>% 
                  filter(Domain == "Bacteria") %>%
                  mutate(Domain = "Bacteria", Taxa = Species) %>%
                  select(Sample, Domain, Taxa, Abundance)

# combine archaea and bacteria
abundance_table_study <- bind_rows(archaea_final, bacteria_component) %>%
                mutate(Taxa = paste(Domain, Taxa, sep = "; ")) %>%
              select(Sample, Taxa, Abundance)

# get a wide version of the abundance_table_study
abundance_table_study <- abundance_table_study %>%
                pivot_wider(names_from = Sample, values_from = Abundance, values_fill = 0) %>%
                as.data.frame() 

###############################################################################################
# calculate the CLR transformed abundances and filter for Archaea
# convert cells that were originally zero back to NA (plotting reasons only)
###########################################################################################

# list the archaea taxa present in the abundance table
archaea_taxa <- grep("Archaea", abundance_table_study[, 1], value = TRUE)

## step 0: a "zero mask" for all rows/columns *before* transformations
##         so we know exactly which entries were originally zero
##         - We'll use the same shape as the main count matrix (rows = taxa, columns = samples)
##         - Then later, we can revert these positions to NA in the final matrices
otu_counts <- as.matrix(abundance_table_study[, -1])
rownames(otu_counts) <- abundance_table_study[, 1]

original_zero_mask <- (otu_counts == 0)

# prepare the CLR data for Archaea
# step 1: compositional transformation: for each sample, divide by total counts
#         a small constant 1e-32 is used to avoid division by zero
otu_compositional <- apply(otu_counts, 2, function(x) {
  x / max(sum(x), 1e-32)
})

# step 2: adjust for zeros: add a small constant ONLY to the cells that are exactly zero,
#         so the pseudocount does not shift the per-sample geometric mean (and therefore
#         does not change the CLR values of the taxa that were actually present)
zero_cells <- (otu_compositional == 0)
if (any(zero_cells)) {
  min_nonzero <- min(otu_compositional[otu_compositional > 0])
  otu_compositional[zero_cells] <- min_nonzero / 2
}

# Step 3: Compute the CLR transformation
#         apply(..., 2, ...) operates column-by-column (samples) and returns the
#         results as columns, so the output is already taxa (rows) x samples (cols).
#         The two t() calls below cancel out and leave the matrix in that orientation.
clr_matrix <- t(apply(otu_compositional, 2, function(x) {
  log(x) - mean(log(x))
}))
clr_matrix <- t(clr_matrix)

# attach original row names
rownames(clr_matrix) <- rownames(otu_counts)

# subset to the Archaea rows
archaea_clr_matrix <- clr_matrix[archaea_taxa, , drop = FALSE]

# step 4: revert the originally-zero cells to NA (plotting only)
#         archaea_zero_mask is already in the same row/column order as
#         archaea_clr_matrix, so it can be used to index directly.
archaea_zero_mask <- original_zero_mask[archaea_taxa, , drop = FALSE]
archaea_clr_matrix[archaea_zero_mask] <- NA

############################################################
#  plot heatmap 
###############################################################

# extract the matrix for the study
taxonomic_table <- archaea_clr_matrix

# add missing taxa ("Archaea; Methanocorpusculum"), so it can match to the public cohorts
extra <- matrix(NA, nrow = 1, ncol = ncol(taxonomic_table))
rownames(extra) <- "Archaea; Methanocorpusculum"
colnames(extra) <- colnames(taxonomic_table)
taxonomic_table <- rbind(taxonomic_table, extra)

# modify the taxonomic table row names
rownames(taxonomic_table) <- gsub("Archaea; ", "", rownames(taxonomic_table))
rownames(taxonomic_table) <- gsub("Rare Archaea", "Other Archaea", rownames(taxonomic_table))
rownames(taxonomic_table) <- gsub("^Methanobrevibacter_A$", "Methanobrevibacter_A sp.", rownames(taxonomic_table))
rownames(taxonomic_table) <- gsub("^Methanosphaera$", "Methanosphaera sp.", rownames(taxonomic_table))
rownames(taxonomic_table) <- gsub("^Methanocorpusculum$", "Methanocorpusculum sp.", rownames(taxonomic_table))
rownames(taxonomic_table) <- gsub("^Methanomassiliicoccus_A$", "Methanomassiliicoccus_A sp.", rownames(taxonomic_table))
rownames(taxonomic_table) <- gsub("^UBA71$", "UBA71 sp.", rownames(taxonomic_table))

# list the column order
order <-c("Methanobrevibacter_A smithii", "Methanobrevibacter_A smithii_A", "Methanobrevibacter_A sp900766745", "Methanobrevibacter_A woesei", "Methanobrevibacter_A sp.",
          "Methanosphaera stadtmanae", "Methanosphaera sp900322125", "Methanosphaera sp.", "Methanomassiliicoccus_A intestinalis", "Methanomassiliicoccus_A sp.", 
          "Methanocorpusculum sp001940805", "Methanocorpusculum sp.", "UBA71 sp006954425", "UBA71 sp905187815", "UBA71 sp.", "UBA71 sp944320485", "Other Archaea")

order <- order[order %in% rownames(taxonomic_table)]

# reorder the rows
taxonomic_table <- taxonomic_table[order,]

# rename archaea 
rownames(taxonomic_table) <- gsub("(.*)", "*\\1*", rownames(taxonomic_table))
rownames(taxonomic_table) <- gsub(" sp.\\*", "* sp.", rownames(taxonomic_table))
rownames(taxonomic_table) <- gsub("\\*Other Archaea\\*", "Other Archaea", rownames(taxonomic_table))
rownames(taxonomic_table) <- gsub(" (sp\\d+)\\*", "* \\1", rownames(taxonomic_table))
rownames(taxonomic_table) <- gsub("\\*UBA71\\*", "UBA71", rownames(taxonomic_table))
rownames(taxonomic_table) <- gsub("\\*Methanobrevibacter_A smithii_A\\*", "**Ca. *M. intestini* (*Methanobrevibacter_A smithii_A*)**", rownames(taxonomic_table))
rownames(taxonomic_table) <- gsub("\\*Methanobrevibacter_A smithii\\*", "***Methanobrevibacter_A smithii***", rownames(taxonomic_table))
rownames(taxonomic_table) <- gsub("\\*Methanobrevibacter_A\\* sp900766745", "***Methanobrevibacter_A* sp900766745**", rownames(taxonomic_table))

# define a color function for the heatmap
col_fun <- colorRamp2(c(-5, 15), c("white", "red"))

# get the title for this study
title <- "FINRISK 2002 \n (n=3,156)"

# create an annotation for the dominant species to add at the top of the heatmap
dominant_species <- apply(taxonomic_table, 2, function(x) rownames(taxonomic_table)[which.max(x)])

# split columns by which dominant species they have
column_groups <- split(seq_len(ncol(taxonomic_table)), dominant_species)

# build the new column order:
new_col_order <- c()

for(sp in rownames(taxonomic_table)) {
  # the columns (samples) for which sp is the dominant taxon
  cols_in_group <- column_groups[[sp]]
  
  # If no columns in this group, skip
  if (length(cols_in_group) == 0) {
    next
  }
  
# subset the matrix to these columns (all archaea rows, so the full row set),
# because we want to cluster these columns among themselves using all row data
submat <- taxonomic_table[, cols_in_group, drop = FALSE]
  
  # If there's only 1 column, there's nothing to cluster, so just keep it
  if (length(cols_in_group) == 1) {
    new_col_order <- c(new_col_order, cols_in_group)
  } else {
    # 1. Compute distance among columns: 
    #    we need the transpose to treat columns as “items to cluster.”
    d <- dist(t(submat))  
    
    # 2. Hierarchical clustering (method of your choice: e.g. "ward.D2" or "complete" etc.)
    hc <- hclust(d, method = "ward.D2")
    
    # 3. Extract the order of columns from the dendrogram
    hc_order <- hc$order
    
    # 4. The actual column indices in the new order
    group_col_order <- cols_in_group[hc_order]
    
    # 5. Append to the final global vector
    new_col_order <- c(new_col_order, group_col_order)
  }
}

# reorder the matrix columns by the newly built order
taxonomic_table <- taxonomic_table[, new_col_order, drop = FALSE]

# also reorder dominant_species so it matches
dominant_species <- dominant_species[new_col_order]

# again, I will tidy up the dominant species names...
dominant_species <- gsub("\\*\\*\\*Methanobrevibacter_A smithii\\*\\*\\*", "*M. smithii*-dominated", dominant_species)
dominant_species <- gsub("\\*\\*Ca. \\*M. intestini\\* .*", "Ca. *M. intestini*-dominated", dominant_species)
dominant_species <- gsub("\\*\\*\\*Methanobrevibacter_A\\* sp900766745\\*\\*", "*M.* sp900766745-dominated", dominant_species)
dominant_species <- ifelse(dominant_species %in% c("*M. smithii*-dominated", "Ca. *M. intestini*-dominated", "*M.* sp900766745-dominated"), dominant_species, "Diverse")

# the factor levels are the rownames in the order they appear in `taxonomic_table`
dominant_species <- factor(dominant_species, levels = c ("*M. smithii*-dominated", "Ca. *M. intestini*-dominated", "*M.* sp900766745-dominated", "Diverse"))

# generate color-blind friendly colors 
species_colors <- c("*M. smithii*-dominated" = "#A472F7", "Ca. *M. intestini*-dominated" = "#17214F", "*M.* sp900766745-dominated"  = "#FFFF54", "Diverse" = "#E6E6E6")

# build annotation for the top of the heatmap
ha_top <- HeatmapAnnotation(
"Archaeal community type" = dominant_species, 
show_annotation_name = FALSE,
col = list("Archaeal community type" = species_colors),
annotation_legend_param = list(
  "Archaeal community type" = list(
    # Manually specify the labels as plotmath expressions:
    labels = c(
      expression(italic("M. smithii") * "-dominated"),
      expression("Ca. " * italic("M. intestini") * "-dominated"),
      expression(italic("M.") * " sp900766745-dominated"),
      "Diverse"
    ),
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 8)
  )
)
)

# build Heatmap object 
ht <- Heatmap(taxonomic_table,
              name = "CLR transformed abundance",
              column_title = title,
              cluster_columns = FALSE,            
              row_labels = gt_render(rownames(taxonomic_table)),
              heatmap_legend_param = list(direction = "horizontal"),
              col = col_fun,
              na_col = "white",  # color for NA cells, here using the same color that col_fun would use for 0
              column_title_gp = gpar(fontsize = 8),
              column_title_rot = 45,
              show_row_names = TRUE,
              row_names_max_width = unit(20, "cm"),
              show_column_names = FALSE,
              row_names_side = "left",
              cluster_rows = FALSE,
              border = "black",
              top_annotation = ha_top)

# save the heatmap to files
svg("heatmap_archcommunity_types.svg", width = 16, height = 8)
draw(ht,
     heatmap_legend_side = "bottom",
     padding = unit(c(2, 40, 2, 2), "mm"))  # bottom, left, top, right
dev.off()

saveRDS(ht, file="heatmap_archcommunity_types.rds")
```
