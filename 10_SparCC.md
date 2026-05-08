# SparCC / FastSpar correlation analysis pipeline

Add the following script to the folder. This script will run fastspar on the data.

```bash
#!/bin/bash
# fastspar_pipeline_dec2025.sh
# A master script to run the fastspar analysis pipeline using Slurm.
#
# Usage:
#   ./fastspar_pipeline_dec2025.sh your_otu_table.tsv
#
# The script performs:
#  1. Submits the main fastspar job.
#  2. Submits the bootstrap job for p-value estimation.
#  3. Submits the job to infer correlations for each bootstrap replicate.
#  4. Submits the p-value calculation job.
#  5. Moves all output files and directories into a consolidated results directory.
#
# Each step is submitted with a dependency so that it starts only after the previous job finishes.

# --- Input Validation ---
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <otu_table_file.tsv>"
    exit 1
fi

OTU_FILE="$1"
if [ ! -f "$OTU_FILE" ]; then
    echo "Error: File '$OTU_FILE' not found!"
    exit 1
fi

# Derive a base name from the OTU file.
# For example, if the file is 'otu_table_6_archaea_1125_bacteria.tsv', then BASE will be 'otu_table_6_archaea_1125_bacteria'
BASE=$(basename "$OTU_FILE" .tsv)

# --- Job Submissions with Dependencies ---

# 1. Run fastspar on the OTU table with 200 iterations.
echo "Submitting fastspar main job..."
jobid1=$(sbatch \
  --partition=standard \
  --time=00-00:20:00 \
  --mem=50GB \
  -c 10 \
  --job-name=fastspar_main \
  --wrap="fastspar --yes --iterations 200 --otu_table $OTU_FILE --correlation ${BASE}_median_correlation.tsv --covariance ${BASE}_median_covariance.tsv --threads 10" | awk '{print $4}')
echo "Submitted fastspar main job with JobID: ${jobid1}"

# 2. Bootstrap step: calculate exact p-values using 10,000 bootstrap replicates.
BOOTSTRAP_DIR="bootstrap_counts_10k_${BASE}"
mkdir -p "$BOOTSTRAP_DIR"
echo "Submitting bootstrap job (10,000 replicates)..."
jobid2=$(sbatch \
  --dependency=afterok:${jobid1} \
  --mem=100G \
  --time=00-02:00:00 \
  -c 30 \
  --job-name=fastspar_bootstrap \
  --wrap="fastspar_bootstrap --otu_table $OTU_FILE --number 10000 --prefix ${BOOTSTRAP_DIR}/${BASE} --threads 30" | awk '{print $4}')
echo "Submitted bootstrap job with JobID: ${jobid2}"

# Brief pause to ensure directories exist.
sleep 5

# 3. Infer correlations for each bootstrap replicate.
BOOTSTRAP_CORR_DIR="bootstrap_correlation_${BASE}"
mkdir -p "$BOOTSTRAP_CORR_DIR"
echo "Submitting bootstrap correlation inference job..."
jobid3=$(sbatch \
  --dependency=afterok:${jobid2} \
  --partition=long \
  --mem=100G \
  --time=02-00:00:00 \
  -c 40 \
  --job-name=fastspar_correlations \
  --wrap="find ${BOOTSTRAP_DIR} -type f | parallel -j 8 fastspar --yes --threads 5 --otu_table {} --correlation ${BOOTSTRAP_CORR_DIR}/cor_{/} --covariance ${BOOTSTRAP_CORR_DIR}/cov_{/} -i 5" | awk '{print $4}')
echo "Submitted bootstrap correlation job with JobID: ${jobid3}"

# 4. Calculate exact p-values from the bootstrap correlations.
echo "Submitting p-values calculation job..."
jobid4=$(sbatch \
  --dependency=afterok:${jobid3} \
  --mem=100G \
  --partition=standard \
  --time=00-02:00:00 \
  -c 20 \
  --job-name=fastspar_pvalues \
  --wrap="fastspar_pvalues --threads 20 --otu_table $OTU_FILE --correlation ${BASE}_median_correlation.tsv --prefix ${BOOTSTRAP_CORR_DIR}/cor_${BASE}_ --permutations 10000 --outfile ${BASE}_pvalues.tsv" | awk '{print $4}')
echo "Submitted p-values calculation job with JobID: ${jobid4}"

# 5. Final consolidation: move all output files and directories to a single results directory.
FINAL_DIR="${BASE}_results"
echo "Submitting consolidation job to move all output files to ${FINAL_DIR}..."
jobid5=$(sbatch \
  --dependency=afterok:${jobid4} \
  --job-name=fastspar_consolidate \
  --wrap="mkdir -p ${FINAL_DIR}; \
          mv ${BASE}_median_correlation.tsv ${BASE}_median_covariance.tsv ${BASE}_pvalues.tsv ${FINAL_DIR}/; \
          mv ${BOOTSTRAP_DIR} ${FINAL_DIR}/; \
          mv ${BOOTSTRAP_CORR_DIR} ${FINAL_DIR}/" | awk '{print $4}')

echo "Submitted consolidation job with JobID: ${jobid5}"
# Optionally, display file counts (these commands execute immediately and may show incomplete results if jobs are still running).
echo "Pipeline submission complete."
echo "Final results will be consolidated in the directory: ${FINAL_DIR}"
```

`your_otu_table.tsv` should be replaced with the OTU table file

To plot use R.

```R
# load libraries
library(tidyverse)

# create a function to read and tidy the results, provide:
fastspar_archaea_extract <- function(path, analysis, archaea_taxa) {
  # read the results
  correlation  <- read.delim(paste0(path, "/", analysis, "_results/", analysis, "_median_correlation.tsv"), header = FALSE, comment = "#")
  rownames(correlation) <- correlation[,1]
  correlation[,1] <- NULL
  colnames(correlation) <- rownames(correlation)
  
  pvalues <- read.delim(paste0(path, "/", analysis, "_results/", analysis, "_pvalues.tsv"), header = FALSE, comment = "#")
  rownames(pvalues) <- pvalues[,1]
  pvalues[,1] <- NULL
  colnames(pvalues) <- rownames(pvalues)
  
  # mask upper triangle with NAs then melt, excluding upper triangle values
  correlation[upper.tri(correlation, diag=TRUE)] <- NA
  pvalues[upper.tri(pvalues, diag=TRUE)] <- NA
  correlation <- reshape2::melt(as.matrix(correlation), na.rm=TRUE, varnames=c('taxa_1', 'taxa_2'), value.name='correlation')
  pvalues <- reshape2::melt(as.matrix(pvalues), na.rm=TRUE, varnames=c('taxa_1', 'taxa_2'), value.name='pvalue')

  # merge correlations and pvalues
  check1 <- table(correlation$taxa_1 == pvalues$taxa_1)
  check2 <- table(correlation$taxa_2 == pvalues$taxa_2)
  if(check1 == FALSE | check2 == FALSE){
    stop("Error: taxa names do not match between correlation and pvalues files")
  }
  correlation_results <- cbind(correlation, pvalues[,3])
  colnames(correlation_results)[4] <- "pvalue"

  # keep only correlations of taxa_2 with Archaea 
  correlation_results <- correlation_results %>% filter(taxa_1 %in% archaea_taxa | taxa_2 %in% archaea_taxa) 
  correlation_results <- correlation_results %>% mutate(is_archaea_vs_archaea = taxa_1 %in% archaea_taxa & taxa_2 %in% archaea_taxa) %>% filter(is_archaea_vs_archaea == FALSE) %>% select(-is_archaea_vs_archaea)

  # organize archaea in a single column
  correlation_results <- correlation_results %>% 
  mutate(taxa_1=as.character(taxa_1), taxa_2=as.character(taxa_2)) %>%
  mutate(archaea = ifelse(taxa_1 %in% archaea_taxa, taxa_1, taxa_2)) %>%
  mutate(taxa_2 = ifelse(taxa_1 %in% archaea_taxa, taxa_2, taxa_1)) %>%
  select(archaea, taxa_2, correlation, pvalue)
      
  return(correlation_results)
}

# Loop through the different studies:
analysis <- c("fastspar_adult_AsnicarF_2021_n845", "fastspar_adult_JieZ_2017_n188", "fastspar_adult_QinN_2014_n121", "fastspar_adult_YachidaS_2019_n178",
                "fastspar_senior_KarlssonFH_2013_n122", "fastspar_adult_BackhedF_2015_n73", "fastspar_adult_MehtaRS_2018_n225", "fastspar_adult_SchirmerM_2016_n300",
                "fastspar_senior_FR02_n405", "fastspar_senior_MetaCardis_2020_a_n236", "fastspar_adult_FR02_n2731", "fastspar_adult_MetaCardis_2020_a_n1091",
                "fastspar_adult_ShaoY_2019_n124", "fastspar_senior_JieZ_2017_n67", "fastspar_senior_YachidaS_2019_n138")

corrMet <- data.frame()
for (s in c("Archaea; Methanobrevibacter_A smithii_A", "Archaea; Methanobrevibacter_A smithii", "Archaea; Methanobrevibacter_A sp900766745")) {
  cat("Analysing Species", s, "\n")
  for(i in analysis){
    cat("Processing", i, "\n")

    # use function to extract the results for the desired Archaea taxa
    tmp <- fastspar_archaea_extract("/labs/workspace/users/camilagv/Archaeal_Colonization/SparCC/01_new_run_dec2025/", i, s) 
    if(nrow(tmp) == 0){
      cat("No correlations found for", s, "in", i, ", check if this species is present in the data\n")
      next
    }

    # split tmp into two dataframes: one with archaea and one with bacteria
    bacteria <- tmp %>% filter(!grepl("Archaea", taxa_2))
    archaea <- tmp %>% filter(grepl("Archaea", taxa_2))

    # adjust p-value for multiple testing
    bacteria$FDR <- p.adjust(bacteria$pvalue, method = "fdr")
    archaea$FDR <- p.adjust(archaea$pvalue, method = "fdr")

    # merge the two dataframes
    tmp <- rbind(bacteria, archaea)

    # add Study
    tmp$Study <- i
    
    # bind results
    corrMet <- rbind(corrMet, tmp)
  }
}

# Create a new column with the combination of the two taxa being compared
corrMet$Comparision <- paste(corrMet$archaea, corrMet$taxa_2, sep = "|")
corrMet <- corrMet %>% select(Comparision, Study, correlation, pvalue, FDR)

# Select comparisions that have significant correlations according to FDR values and also |correlation| > 0.4
comparision_to_filter <- corrMet %>% filter(FDR < 0.05 & abs(correlation) > 0.4) %>% arrange(desc(correlation)) %>% pull(Comparision) %>% unique()
length(comparision_to_filter) 

wide_corr <- corrMet %>%
  mutate(correlation = round(correlation, 2), pvalue = round(pvalue, 4), FDR = round(FDR, 4)) %>%
  filter(Comparision %in% comparision_to_filter) %>%
  mutate(Comparision = factor(Comparision, levels = comparision_to_filter)) %>%
  arrange(Comparision, Study) %>%
  pivot_wider(
    names_from = Study,
    values_from = c(correlation, pvalue, FDR),
    names_glue = "{.value} {Study}",
    values_fn = ~ .[1]
  )

# Extract the bacteria names with significant correlations
significant_bacteria <- comparision_to_filter %>% 
  str_split("\\|") %>% 
  unlist() %>% 
  unique() %>% 
  .[!grepl("Archaea", .)] %>% 
  .[!grepl("Methanobrevibacter", .)]

# Tidy and filter for plotting
corrMet_filt <- corrMet  %>%
            mutate(Bacteria = str_extract(Comparision, "Bacteria; (.*)")) %>%
            filter(Bacteria %in% significant_bacteria)

corrMet_filt <- corrMet_filt %>% 
                mutate(Archaea = gsub("Archaea; (.*)\\|.*", "\\1", Comparision)) %>% 
                mutate(Bacteria = gsub("Bacteria; (.*)", "\\1", Bacteria)) %>%
                select(Archaea, Bacteria, Study, correlation, pvalue, FDR) 

# Tidy study names
corrMet_filt$age_group <- ifelse(grepl("senior", corrMet_filt$Study), "Seniors", "Adults")
corrMet_filt$Study <- gsub("fastspar_(adult|senior)_(.*)_n(.*)", "\\2 (n=\\3)", corrMet_filt$Study)
corrMet_filt$Study <- gsub("n=2731", "n=2,731", corrMet_filt$Study)
corrMet_filt$Study <- gsub("n=1091", "n=1,091", corrMet_filt$Study)
length(unique(corrMet_filt$Study)) # 15

# Tidy archaea names
corrMet_filt$Archaea <- gsub("Methanobrevibacter_A smithii_A", "M. intestini", corrMet_filt$Archaea) 
corrMet_filt$Archaea <- gsub("Methanobrevibacter_A smithii", "M. smithii", corrMet_filt$Archaea) 
corrMet_filt$Archaea <- gsub("Methanobrevibacter_A sp900766745", "M. sp900766745", corrMet_filt$Archaea) 

# ensure factors are in order 
corrMet_filt %>% filter(Archaea=="M. smithii") %>% arrange(correlation) %>% pull(Bacteria) %>% unique() -> significant_bacteria_levels
significant_bacteria_levels <- c(significant_bacteria_levels[significant_bacteria_levels != "Christensenella"], "Christensenella")
corrMet_filt$Bacteria <- factor(corrMet_filt$Bacteria, levels=significant_bacteria_levels)

adult_studies <- corrMet_filt %>% filter(age_group == "Adults") %>% pull(Study) %>% unique()
senior_studies <- corrMet_filt %>% filter(age_group == "Seniors") %>% pull(Study) %>% unique()

# make sure Study is ordered sensibly
corrMet_filt_adults <- corrMet_filt %>% filter(Study %in% adult_studies)

corrMet_filt_adults$Study <- factor(
  corrMet_filt_adults$Study,
  levels = unique(corrMet_filt_adults$Study[corrMet_filt_adults$Study %in% adult_studies])
)

svg("fastspar_metagenomics_correlations_adult_by_archaea.svg", width = 8, height = 5)
ggplot(
  corrMet_filt_adults,
  aes(x = Study, y = Bacteria, fill = correlation)
) +
  geom_tile(color = "grey80") +
  facet_grid(. ~ Archaea, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low      = "#0000F5",
    mid      = "white",
    high     = "#C00000",
    midpoint = 0,
    name     = "ρ"
  ) +
  geom_text(
    data = corrMet_filt_adults %>% 
      filter(FDR < 0.05),
    aes(label = "*"),
    color = "white",
    size  = 4
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(face = "italic"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank(),
    strip.text  = element_text(face = "italic", size = 8)
  )
dev.off()


# make sure Study is ordered sensibly
corrMet_filt_seniors <- corrMet_filt %>% filter(Study %in% senior_studies)

corrMet_filt_seniors$Study <- factor(
  corrMet_filt_seniors$Study,
  levels = unique(corrMet_filt_seniors$Study[corrMet_filt_seniors$Study %in% senior_studies])
)

svg("fastspar_metagenomics_correlations_senior_by_archaea.svg", width = 5, height = 5)
ggplot(
  corrMet_filt_seniors,
  aes(x = Study, y = Bacteria, fill = correlation)
) +
  geom_tile(color = "grey80") +
  facet_grid(. ~ Archaea, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low      = "#0000F5",
    mid      = "white",
    high     = "#C00000",
    midpoint = 0,
    name     = "ρ"
  ) +
  geom_text(
    data = corrMet_filt_seniors %>% 
      filter(FDR < 0.05),
    aes(label = "*"),
    color = "white",
    size  = 4
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(face = "italic"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank(),
    strip.text  = element_text(face = "italic", size = 8)
  )
dev.off()
```