# Estimate archaeal genome coverage 

We want to compute among archaea-positive samples only: per-species estimated coverage, genus-level upper bound and summaries per cohort and pooled

coverage = (reads_species x read_length x k) / genome_size_bp 

--> k = 2 for a read pair
--> read length varies across CMD (some are 100 bp, some 150 bp)

FINRISK 2002 will not be included in this analysis, it is extremely low coverage

```R
# libraries
library(tidyverse)
library(phyloseq)

# ---- inputs ----

physeq_list <- readRDS("physeq_list_metagenomes_with_metadata_and_single_sample_per_individual.rds")

selected <- c("KarlssonFH_2013-145samples", "AsnicarF_2021-1088samples", "MetaCardis_2020_a-1726samples",
              "MehtaRS_2018-308samples", "SchirmerM_2016-471samples",
              "JieZ_2017-385samples", "QinN_2014-237samples", "YachidaS_2019-609samples",
              "BackhedF_2015-200samples", "HMP_2019_ibdmdb-130samples", "VatanenT_2016-212samples",
              "LiJ_2017-196samples", "ShaoY_2019-753samples", "BrooksB_2017-30samples")

physeq_list <- physeq_list[selected]

target_species <- c("Methanobrevibacter_A smithii",
                    "Methanobrevibacter_A smithii_A",
                    "Methanobrevibacter_A sp900766745")

# Median genome sizes (bp) from Table S7
genome_sizes <- tribble(
  ~Species,                            ~genome_bp,
  "Methanobrevibacter_A smithii",         1.725e6,
  "Methanobrevibacter_A smithii_A",       1.905e6,
  "Methanobrevibacter_A sp900766745",     1.573e6
)

# Genome size used for the genus-level upper bound: SMALLEST of the three
GENUS_GENOME_BP <- min(genome_sizes$genome_bp)

# Per-cohort read length (bp)
read_lengths <- tribble(
  ~Study,                            ~read_len,
  "KarlssonFH_2013-145samples",         100,
  "AsnicarF_2021-1088samples",          150,
  "MetaCardis_2020_a-1726samples",      150,
  "MehtaRS_2018-308samples",            100,
  "SchirmerM_2016-471samples",          100,
  "JieZ_2017-385samples",               100,
  "QinN_2014-237samples",               100,
  "YachidaS_2019-609samples",           150,
  "BackhedF_2015-200samples",           100,
  "HMP_2019_ibdmdb-130samples",         100,
  "VatanenT_2016-212samples",           100,
  "LiJ_2017-196samples",                100,
  "ShaoY_2019-753samples",              100,
  "BrooksB_2017-30samples",             150
)

K_PAIRED <- 2   # 2 if counts are read PAIRS, 1 if single reads. VERIFY.

# ---- extract raw archaeal read counts (no relative-abundance transform) ----
archaeal_counts <- map_dfr(names(physeq_list), function(study) {
  cat("Counting reads:", study, "\n")
  phy <- physeq_list[[study]]

  # keep counts as-is; psmelt on the untransformed object gives raw read counts
  m <- psmelt(phy)

  # per-sample, per-species Methanobrevibacter counts
  spp <- m %>%
    filter(Species %in% target_species) %>%
    group_by(Sample, Species) %>%
    summarise(reads = sum(Abundance, na.rm = TRUE), .groups = "drop")

  # per-sample total archaeal reads, used to define "archaea-positive"
  arch_total <- m %>%
    filter(Domain == "Archaea") %>%
    group_by(Sample) %>%
    summarise(archaeal_reads = sum(Abundance, na.rm = TRUE), .groups = "drop")

  spp %>%
    left_join(arch_total, by = "Sample") %>%
    mutate(Study = study)
})

# ---- 1) per-species coverage ----

cov_species <- archaeal_counts %>%
  filter(archaeal_reads > 0) %>%              # archaea-positive samples only
  left_join(genome_sizes, by = "Species") %>%
  left_join(read_lengths, by = "Study") %>%
  mutate(coverage = reads * read_len * K_PAIRED / genome_bp)

# ---- 2) genus-level upper bound ----

cov_genus <- archaeal_counts %>%
  filter(archaeal_reads > 0) %>%
  group_by(Study, Sample) %>%
  summarise(genus_reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
  filter(genus_reads > 0) %>%                 # Methanobrevibacter-positive
  left_join(read_lengths, by = "Study") %>%
  mutate(coverage_ub = genus_reads * read_len * K_PAIRED / GENUS_GENOME_BP)

# ---- 3) summaries ----

summ_cov <- function(x) {
  tibble(
    n         = length(x),
    median    = median(x),
    IQR_Lower = as.numeric(quantile(x, 0.25, names = FALSE)),
    IQR_Upper = as.numeric(quantile(x, 0.75, names = FALSE)),
    max       = max(x),
    pct_ge_1  = 100 * mean(x >= 1),
    pct_ge_10 = 100 * mean(x >= 10),
    pct_ge_20 = 100 * mean(x >= 20)
  )
}

# by cohort x species
cov_by_study_species <- cov_species %>%
  group_by(Study, Species) %>%
  summarise(summ_cov(coverage), .groups = "drop")

# by cohort, genus-level upper bound
cov_by_study_genus <- cov_genus %>%
  group_by(Study) %>%
  summarise(summ_cov(coverage_ub), .groups = "drop")

# pooled across all cohorts (the headline numbers)
cov_pooled_genus   <- summ_cov(cov_genus$coverage_ub)
cov_pooled_species <- cov_species %>%
  group_by(Species) %>%
  summarise(summ_cov(coverage), .groups = "drop")

write.csv(cov_by_study_species, "coverage_by_study_species.csv", row.names = FALSE)
write.csv(cov_by_study_genus, "coverage_by_study_genus_upperbound.csv", row.names = FALSE)

cov_by_study_species %>% 
filter(Study == "AsnicarF_2021-1088samples") %>%
as.data.frame()
#                       Study                          Species   n      median
# 1 AsnicarF_2021-1088samples     Methanobrevibacter_A smithii 846 0.159913043
# 2 AsnicarF_2021-1088samples   Methanobrevibacter_A smithii_A 846 0.000000000
# 3 AsnicarF_2021-1088samples Methanobrevibacter_A sp900766745 846 0.009154482
#   IQR_Lower IQR_Upper       max  pct_ge_1  pct_ge_10 pct_ge_20
# 1         0 7.5953913 203.14539 35.342790 23.6406619 17.021277
# 2         0 0.2777953 100.17024 18.557920 12.6477541  8.628842
# 3         0 0.1192467  14.30102  5.082742  0.3546099  0.000000

# Recompute on detection-positive samples:
cov_species |>
  filter(reads > 0) |>
  group_by(Study, Species) |>
  summarise(n = n(), median = median(coverage),
            pct_ge_10 = 100*mean(coverage >= 10),
            pct_ge_20 = 100*mean(coverage >= 20),
            max = max(coverage), .groups = "drop") %>% 
filter(Study == "AsnicarF_2021-1088samples") %>%
as.data.frame()
#                       Study                          Species   n    median
# 1 AsnicarF_2021-1088samples     Methanobrevibacter_A smithii 534 1.9915652
# 2 AsnicarF_2021-1088samples   Methanobrevibacter_A smithii_A 414 0.3059055
# 3 AsnicarF_2021-1088samples Methanobrevibacter_A sp900766745 424 0.1191036
#    pct_ge_10 pct_ge_20       max
# 1 37.4531835  26.96629 203.14539
# 2 25.8454106  17.63285 100.17024
# 3  0.7075472   0.00000  14.30102
```