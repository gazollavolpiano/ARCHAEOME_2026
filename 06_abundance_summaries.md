# Generate abundance summaries for CMD cohorts

We will compute the values among SAMPLES WITH DETECTION for Methanobrevibacter smithii, Methanobrevibacter intestini or Methanobrevibacter sp900766745. This means that we will subset the data to include only samples where the taxon is present and then calculate the median, IQR, min, and max relative abundance values for those samples. 

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
```
