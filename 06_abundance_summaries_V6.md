# Generate abundance summaries for CMD cohorts

We will compute the values among SAMPLES WITH DETECTION for Methanobrevibacter smithii, Methanobrevibacter intestini or Methanobrevibacter sp900766745. This means that we will subset the data to include only samples where the taxon is present and then calculate the median, IQR, min, and max relative abundance values for those samples. 

We are getting two different types of relative abundance values:
1) Methanobrevibacter values are fractions of the whole bacterial + archaeal community: dividing each taxon by the total community per sample (x / sum(x) over everything)
2) Methanobrevibacter values are fractions within the archaea, we divide by the total archaeal abundance per sample instead
3) We also calculate the proportion of samples where each taxon is the sole archaeon detected, among samples with archaeal detection, and the median relative abundance of that taxon among samples where it is co-detected with other archaeal taxa

```R
# libraries
library(tidyverse)
library(phyloseq)

# read CMD phyloseq object
physeq_list <- readRDS("physeq_list_metagenomes_with_metadata_and_single_sample_per_individual.rds")

# filter for cohorts of interest 
selected <- c("KarlssonFH_2013-145samples", "AsnicarF_2021-1088samples", "MetaCardis_2020_a-1726samples",
              "MehtaRS_2018-308samples", "SchirmerM_2016-471samples",
              "JieZ_2017-385samples", "QinN_2014-237samples", "YachidaS_2019-609samples",
              "BackhedF_2015-200samples", "HMP_2019_ibdmdb-130samples", "VatanenT_2016-212samples",
              "LiJ_2017-196samples", "ShaoY_2019-753samples", "BrooksB_2017-30samples")

physeq_list <- physeq_list[selected]

# define target species 
target_species <- c("Methanobrevibacter_A smithii", "Methanobrevibacter_A smithii_A", "Methanobrevibacter_A sp900766745")

# helper: summarise distribution among detection-positive values only
summ_detected <- function(x) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) == 0) {
    return(tibble(
      n_detected = 0L,
      median = NA_real_, min = NA_real_, max = NA_real_,
      IQR_Lower = NA_real_, IQR_Upper = NA_real_
    ))
  }
  tibble(
    n_detected = length(x),
    median = median(x),
    min = min(x),
    max = max(x),
    IQR_Lower = as.numeric(quantile(x, 0.25, names = FALSE)),
    IQR_Upper = as.numeric(quantile(x, 0.75, names = FALSE))
  )
}

relative_archaea_abundance <- map_dfr(names(physeq_list), function(study) {
  cat("Processing study:", study, "\n")

  phy <- physeq_list[[study]]

  # Relative abundance per sample (%). Guard against empty samples.
  phy <- transform_sample_counts(phy, function(x) {
    s <- sum(x)
    if (s == 0) return(x)
    (x / s) * 100
  })

  m <- psmelt(phy)

  # 1) Species-specific summaries among samples where THAT species is detected (>0)
  species_sum <- m %>%
    filter(Species %in% target_species) %>%
    group_by(Sample, Species) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(Species) %>%
    summarise(summ_detected(Abundance), .groups = "drop") %>%
    tidyr::complete(
      Species = target_species,
      fill = list(
        n_detected = 0L,
        median = NA_real_, min = NA_real_, max = NA_real_,
        IQR_Lower = NA_real_, IQR_Upper = NA_real_
      )
    ) %>%
    mutate(Study = study, Taxon = Species) %>%
    select(Study, Taxon, median, min, max, IQR_Lower, IQR_Upper, n_detected)

  # 2) Total Archaea per sample, summarised among samples with archaeal detection (>0)
  total_arch_per_sample <- m %>%
    filter(Domain == "Archaea") %>%
    group_by(Sample) %>%
    summarise(Total_Archaea = sum(Abundance, na.rm = TRUE), .groups = "drop")

  total_sum <- total_arch_per_sample %>%
    summarise(summ_detected(Total_Archaea)) %>%
    mutate(
      Study = study,
      Taxon = "Total Archaea (summed across taxa)"
    ) %>%
    select(Study, Taxon, median, min, max, IQR_Lower, IQR_Upper, n_detected)

  bind_rows(species_sum, total_sum)
})

write.csv(relative_archaea_abundance, "relative_archaea_abundance_summary_by_study_detected_only.csv")

# ---- Within-archaea relative abundance (each sample's archaeome = 100%) ----
relative_within_archaea <- map_dfr(names(physeq_list), function(study) {
  cat("Processing study (within-archaea):", study, "\n")
  phy <- physeq_list[[study]]

  # keep only Archaea, THEN renormalise -> fractions of the archaeome
  phy_arch <- subset_taxa(phy, Domain == "Archaea")
  phy_arch <- transform_sample_counts(phy_arch, function(x) {
    s <- sum(x)
    if (s == 0) return(x)        # samples with no archaea stay 0 -> excluded below
    (x / s) * 100
  })

  m <- psmelt(phy_arch)

  m %>%
    filter(Species %in% target_species) %>%
    group_by(Sample, Species) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(Species) %>%
    summarise(summ_detected(Abundance), .groups = "drop") %>%
    tidyr::complete(
      Species = target_species,
      fill = list(n_detected = 0L, median = NA_real_, min = NA_real_,
                  max = NA_real_, IQR_Lower = NA_real_, IQR_Upper = NA_real_)
    ) %>%
    mutate(Study = study, Taxon = Species) %>%
    select(Study, Taxon, median, min, max, IQR_Lower, IQR_Upper, n_detected)
})

write.csv(relative_within_archaea, "relative_within_archaea_abundance_summary_by_study_detected_only.csv")

# ---- Table S3C: composition among co-detected samples ----

# among detection-positive samples, how often is this taxon the sole archaeon?
per_sample_within_archaea <- map_dfr(names(physeq_list), function(study) {
  cat("Per-sample (within-archaea):", study, "\n")
  phy <- physeq_list[[study]]

  phy_arch <- subset_taxa(phy, Domain == "Archaea")
  phy_arch <- transform_sample_counts(phy_arch, function(x) {
    s <- sum(x)
    if (s == 0) return(x)
    (x / s) * 100
  })

  psmelt(phy_arch) %>%
    group_by(Sample, Species) %>%
    summarise(WithinArchaea = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(Study = study)
})

sole_detection_summary <- per_sample_within_archaea %>%
  filter(WithinArchaea > 0) %>%
  group_by(Study, Sample) %>%
  mutate(n_arch_taxa = n()) %>%
  ungroup() %>%
  filter(Species %in% target_species) %>%
  group_by(Study, Species) %>%
  summarise(
    n_detected = n(),
    prop_sole  = mean(n_arch_taxa == 1),
    median_co  = median(WithinArchaea[n_arch_taxa > 1]),
    n_co       = sum(n_arch_taxa > 1),
    .groups = "drop"
  )

table_S3C <- sole_detection_summary %>%
  tidyr::complete(
    Study   = Study,
    Species = target_species,
    fill = list(n_detected = 0L, n_co = 0L,
                prop_sole = NA_real_, median_co = NA_real_)
  ) %>%
  mutate(
    Study      =  Study,
    Species    = factor(Species, levels = target_species),
    pct_sole   = round(prop_sole * 100, 1),
    median_co  = round(median_co, 2),
    median_co  = if_else(n_co == 0, NA_real_, median_co),
    flag_low_n = if_else(n_co > 0 & n_co < 5, "*", "")
  ) %>%
  arrange(Study, Species) %>%
  select(Study, Taxon = Species,
         n_detected, pct_sole, n_co, median_co, flag_low_n)

write.csv(table_S3C, "table_S3C_within_archaea_codetected.csv", row.names = FALSE)
```
