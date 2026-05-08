# Create abundance tables for the MICOM model

We want a table with columns: kingdom	xxx genus	species	sample_id	abundance study. Abundance is expressed as count of reads assigned to that taxon in that sample.

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

# table to store the abundance tables for each study
abundance_per_study <- data.frame()

for(study_name in studies) {

  # Print the study name for CMD and create the taxa_counts
  cat("Processing study:", study_name, "\n")

  # Filter for the study being analyzed
  phyloseq_obj <- physeq_list[[study_name]]

  # Melt the phyloseq object and remove rows with zero abundance
  melt_phyloseq <- psmelt(phyloseq_obj) %>% filter(Abundance > 0) %>% select(Domain, Phylum, Class, Order, Family, Genus, Species, SampleID, Abundance)
  melt_phyloseq$Study <- study_name

  # Add to the abundance table for the study
  abundance_per_study <- rbind(abundance_per_study, melt_phyloseq)
}

abundance_per_study %>% 
filter(Domain == "Archaea") %>%
select(SampleID, Study) %>%
distinct() %>%
group_by(Study) %>%
summarise(n_samples = n()) %>%
arrange(desc(n_samples)) 

# Save the abundance table for all studies
write.table(abundance_per_study, "abundance_per_study_for_Alex_30_Mar_2026.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Save as a R object for later use
saveRDS(abundance_per_study, "abundance_per_study_for_Alex_30_Mar_2026.rds")

# Save what is the dominant archaeal species per sample in each study?
dominant_archaea_per_sample <- abundance_per_study %>%
  filter(Domain == "Archaea") %>%
  group_by(Study, SampleID) %>%
  mutate(total_archaeal_reads = sum(Abundance)) %>%
  mutate(Percent_of_Archaeal_Reads = (Abundance / total_archaeal_reads)*100) %>%
  slice_max(order_by = Percent_of_Archaeal_Reads, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Study, SampleID, Dominant_Species = Species, Percent_of_Archaeal_Reads) 

write.table(dominant_archaea_per_sample, "dominant_archaea_per_sample_for_Alex_30_Mar_2026.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
```

# Creating a shorter list for MICOM modeling (PREDICT 1)

- We need to do a filtering down to a manageable size for MICOM modeling 
- Target sample size: 10, 50 and 100 samples per species (sp900766745 and M. smithii) to get initial results quickly?
- Metadata should ideally include sequencing depth, sex, BMI, age, and geographic region — to allow homogenization of the subset

We will use only AsnicarF_2021-1088samples (PREDICT 1).

```R
# load libraries
library(tidyverse)
library(MatchIt) # For matching samples based on covariates

# Read dominant_archaea_per_sample so we can filter the data previously created
dominant_archaea_per_sample <- read.table("dominant_archaea_per_sample_for_Alex_30_Mar_2026.tsv", sep = "\t", header = TRUE)

#################################################################
# AsnicarF_2021-1088samples
#################################################################

# Select samples from the dominant_archaea_per_sample table for the current study
samples_in_study <- dominant_archaea_per_sample %>% filter(Study == "AsnicarF_2021-1088samples") %>% pull(SampleID)
length(samples_in_study) # 846 samples in this study, correct!

# Filter for the study being analyzed and only keep samples that are in the dominant_archaea_per_sample table
metadata <- readRDS("AsnicarF_2021/metadata_1088_samples_1088_individuals.rds") 

# Check data
head(metadata, 2)
#      sample_id        subject_id age age_category gender country number_reads
# 1 SAMEA7041133 predict1_MTG_0001  36        adult female     USA     27999010
# 2 SAMEA7041134 predict1_MTG_0002  24        adult   male     GBR     28573197
#   number_bases minimum_read_length NCBI_accession      BMI family
# 1   4227850510                 151     ERR4330026 37.26849   <NA>
# 2   4314552747                 151     ERR4330027 23.26859   <NA>

# Filter for samples in samples_in_study
metadata <- metadata %>% filter(NCBI_accession %in% samples_in_study)

# Filter for adults only
metadata <- metadata %>% filter(age_category == "adult")
dim(metadata) # 845  12

# Filter for GBR country only (since we have very few samples from USA)
metadata <- metadata %>% filter(country == "GBR")

# Add dominant archaeal species information for this study
metadata <- metadata %>% left_join(dominant_archaea_per_sample, by = c("NCBI_accession" = "SampleID"))
metadata <- metadata %>% rename(SampleID = NCBI_accession)

# Filter for samples with dominant archaeal species being either M. smithii or sp900766745
metadata <- metadata %>% filter(Dominant_Species %in% c("Methanobrevibacter_A smithii", "Methanobrevibacter_A sp900766745"))

# Check how many samples we have for each dominant archaeal species
metadata$Dominant_Species %>% table()
# Methanobrevibacter_A smithii Methanobrevibacter_A sp900766745 
#                          291                              208 

# BMI group breakdown (WHO thresholds)
metadata <- metadata %>%
  mutate(bmi_group = case_when(
    BMI < 18.5             ~ "underweight",
    BMI >= 18.5 & BMI < 25 ~ "normal",
    BMI >= 25  & BMI < 30  ~ "overweight",
    BMI >= 30              ~ "obese",
    TRUE                   ~ NA_character_
  ))

# Remove missing values on BMI, gender, number of reads
metadata <- metadata %>%  filter(!is.na(BMI), !is.na(gender), !is.na(number_reads), !is.na(age))
dim(metadata) # 488 x 16

# Remove underweight (too few samples)
metadata$bmi_group %>% table()
#  normal       obese  overweight underweight 
#     271          71         146          11 

metadata <- metadata %>% filter(bmi_group != "underweight")

# Create a binary indicator (1 for the target species, 0 for the match)
metadata$is_sp900766745 <- ifelse(metadata$Dominant_Species == "Methanobrevibacter_A sp900766745", 1, 0)

# Perform the matching inside each BMI group (normal, overweight and obese) based on the other covariates 
# Because we think that the biological function of these archaea shifts in each BMI condition we will not treat BMI merely as a covariate to be "balanced" globally, so we will create separate matched datasets for each BMI group. 

# COVARIATES TO MATCH ON: number of reads, age, BMI (numerical), gender

df_matched <- data.frame() # Create an empty dataframe to store the matched samples for all BMI groups
for(bmi_group_selected in c("normal", "overweight", "obese")) {
    
    # Filter the data for the current BMI group
    df <- metadata %>% filter(bmi_group == bmi_group_selected) 
    
    set.seed(123) # For reproducibility

    # matching using MatchIt
    m.out <- matchit(is_sp900766745 ~ number_reads + age + BMI, 
                  data = df, 
                  method = "nearest",      # Nearest neighbor matching
                  distance = "glm",        # Uses logistic regression to compute propensity scores
                  exact = ~  gender,       # Forces exact matches 
                  ratio = 1,               # 1-to-1 matching
                  caliper = 0.035)          # Caliper to ensure good matches on the propensity score

    # Extract the new matched dataframe
    df <- match.data(m.out) %>% as.data.frame()

    # Add dataset_stata information back to the matched dataframe
    df$dataset_stata <- paste0(bmi_group_selected, "_AsnicarF_2021-1088samples")
    # Append the matched samples for the current BMI group to the overall matched dataframe
    df_matched <- rbind(df_matched, df)
}

# Check how many samples we have for each dominant archaeal species in the matched dataset
df_matched %>% 
  group_by(bmi_group, Dominant_Species) %>%
  summarise(n_samples = n()) %>% 
  arrange(bmi_group)%>%
  as.data.frame()

#    bmi_group                 Dominant_Species n_samples
# 1     normal     Methanobrevibacter_A smithii        68
# 2     normal Methanobrevibacter_A sp900766745        68
# 3      obese     Methanobrevibacter_A smithii        10
# 4      obese Methanobrevibacter_A sp900766745        10
# 5 overweight     Methanobrevibacter_A smithii        30
# 6 overweight Methanobrevibacter_A sp900766745        30

# Rename dataset_stata to be strata_1, strata_2 and strata_3 for the three BMI groups
df_matched$dataset_stata <- case_when(
  df_matched$dataset_stata == "normal_AsnicarF_2021-1088samples" ~ "strata_1",
  df_matched$dataset_stata == "overweight_AsnicarF_2021-1088samples" ~ "strata_2",
  df_matched$dataset_stata == "obese_AsnicarF_2021-1088samples" ~ "strata_3",
  TRUE ~ df_matched$dataset_stata
)

# Organize columns
df_matched <- df_matched %>% 
  select(dataset_stata, subclass, distance, bmi_group, SampleID, Study, Dominant_Species, number_reads, age, BMI, gender) %>% 
  arrange(dataset_stata, subclass)

# Save the matched dataset for AsnicarF_2021-1088samples
write.table(df_matched, "matched_samples_AsnicarF_2021_for_Alex_07_Apr_2026.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# #################################################################
# # MetaCardis_2020_a-1726samples
# #################################################################

# # Select samples from the dominant_archaea_per_sample table for the current study
# samples_in_study <- dominant_archaea_per_sample %>% filter(Study == "MetaCardis_2020_a-1726samples") %>% pull(SampleID)
# length(samples_in_study) # 1327 samples in this study, correct!

# # Filter for the study being analyzed and only keep samples that are in the dominant_archaea_per_sample table
# metadata <- readRDS("/labs/mericlab/curatedMetagenomicData/MetaCardis_2020_a/metadata_1727_samples_1727_individuals.rds") 

# # Need to prepare read count data for this study to be able to do the matching based on sequencing depth (number of reads)
# fastq <- read.table("/labs/mericlab/curatedMetagenomicData/MetaCardis_2020_a/fastq_stats_1726_samples.txt", sep = "\t", header = FALSE) 
# fastq$SampleID <- gsub(".*fastq\\/(.*)\\_.\\.fastq", "\\1", fastq$V1)
# setdiff(samples_in_study, fastq$SampleID) # All samples in the study are in the fastq file, good!
# fastq$number_reads <- as.numeric(fastq$V4) # Convert number of reads to numeric
# metadata <- metadata %>% left_join(fastq %>% select(SampleID, number_reads), by = c("subject_id" = "SampleID"))

# # Check data
# head(metadata, 2)
# #      sample_id   subject_id antibiotics_current_use study_condition disease
# # 1 M0x10MCx1134 M0x10MCx1134                      no             IGT  IGT;MS
# # 2 M0x10MCx1135 M0x10MCx1135                      no             T2D     T2D
# #   age_category gender country number_bases
# # 1        adult female     FRA   2995319758
# # 2        adult   male     FRA   3003448107
# #                                                                 NCBI_accession
# # 1 ERR4086424;ERR4086425;ERR4086426;ERR4086427;ERR4086428;ERR4086429;ERR4086430
# # 2                       ERR4563248;ERR4563249;ERR4563250;ERR4563251;ERR4563252
# #        BMI                           treatment location disease_subtype
# # 1 44.36589       antihta;thiazidique;at2_inhib    Paris            <NA>
# # 2 27.71967 antidiab;su;metformin;dppiv;insulin    Paris            <NA>
# #   triglycerides hba1c smoke bristol_score hsCRP      LDL number_reads
# # 1      97.94071     6    no             1   1.3 140.6737     21164177
# # 2            NA    NA    no            NA    NA       NA     21367137

# # save metadata as rds for later use
# saveRDS(metadata, "metadata_1727_samples_1727_individuals.rds")

# # Filter for samples in samples_in_study
# metadata <- metadata %>% filter(subject_id %in% samples_in_study)

# # Filter for adults only
# metadata <- metadata %>% filter(age_category == "adult")
# dim(metadata) # 1091   21

# # Filter DNK out because it has no BMI
# metadata$country %>% table()
# # DEU DNK FRA 
# # 291 422 378 

# metadata <- metadata %>% filter(country != "DNK")

# # Add dominant archaeal species information for this study
# metadata <- metadata %>% left_join(dominant_archaea_per_sample, by = c("subject_id" = "SampleID"))
# metadata <- metadata %>% rename(SampleID = subject_id)

# # Filter for samples with dominant archaeal species being either M. smithii or sp900766745
# metadata <- metadata %>% filter(Dominant_Species %in% c("Methanobrevibacter_A smithii", "Methanobrevibacter_A sp900766745"))

# # Check how many samples we have for each dominant archaeal species
# metadata$Dominant_Species %>% table()
# # Methanobrevibacter_A smithii Methanobrevibacter_A sp900766745 
# #                          232                              262 

# # BMI group breakdown (WHO thresholds)
# metadata <- metadata %>%
#   mutate(bmi_group = case_when(
#     BMI < 18.5             ~ "underweight",
#     BMI >= 18.5 & BMI < 25 ~ "normal",
#     BMI >= 25  & BMI < 30  ~ "overweight",
#     BMI >= 30              ~ "obese",
#     TRUE                   ~ NA_character_
#   ))

# # Remove missing values on BMI, gender, number of reads
# metadata <- metadata %>%  filter(!is.na(BMI), !is.na(gender), !is.na(number_reads))
# dim(metadata) # 492  25

# # Check group breakdowns
# metadata$bmi_group %>% table()
#     # normal      obese overweight 
#     #     83        331         78 

# # metadata <- metadata %>% filter(bmi_group != "underweight") # no need

# # Create a binary indicator (1 for the target species, 0 for the match)
# metadata$is_sp900766745 <- ifelse(metadata$Dominant_Species == "Methanobrevibacter_A sp900766745", 1, 0)

# # Perform the matching inside each BMI group (normal, overweight and obese) based on the other covariates 
# # Because we think that the biological function of these archaea shifts in each BMI condition we will not treat BMI merely as a covariate to be "balanced" globally, so we will create separate matched datasets for each BMI group. 

# # COVARIATES TO MATCH ON: number of reads, BMI (numerical), gender, country, antibiotics_current_use, study_condition, country

# df_matched <- data.frame() # Create an empty dataframe to store the matched samples for all BMI groups
# for(bmi_group_selected in c("normal", "overweight", "obese")) {
    
#     # Filter the data for the current BMI group
#     df <- metadata %>% filter(bmi_group == bmi_group_selected) 
    
#     set.seed(123) # For reproducibility

#     # matching using MatchIt
#     m.out <- matchit(is_sp900766745 ~  BMI + number_reads, 
#                   data = df, 
#                   method = "nearest",      # Nearest neighbor matching
#                   distance = "glm",        # Uses logistic regression to compute propensity scores
#                   exact = ~  gender + country + antibiotics_current_use + study_condition,       # Forces exact matches 
#                   ratio = 1,               # 1-to-1 matching
#                   caliper = 0.25)          # Caliper to ensure good matches on the propensity score

#     # Extract the new matched dataframe
#     df <- match.data(m.out) %>% as.data.frame()

#     # Add dataset_stata information back to the matched dataframe
#     df$dataset_stata <- paste0(bmi_group_selected, "_MetaCardis_2020_a-1726samples")
#     # Append the matched samples for the current BMI group to the overall matched dataframe
#     df_matched <- rbind(df_matched, df)
# }

# # Check how many samples we have for each dominant archaeal species in the matched dataset
# df_matched %>% 
#   group_by(bmi_group, Dominant_Species) %>%
#   summarise(n_samples = n()) %>% 
#   arrange(bmi_group)%>%
#   as.data.frame()
# #    bmi_group                 Dominant_Species n_samples
# # 1     normal     Methanobrevibacter_A smithii        13
# # 2     normal Methanobrevibacter_A sp900766745        13
# # 3      obese     Methanobrevibacter_A smithii        84
# # 4      obese Methanobrevibacter_A sp900766745        84
# # 5 overweight     Methanobrevibacter_A smithii        12
# # 6 overweight Methanobrevibacter_A sp900766745        12

# # Rename dataset_stata to be strata_1, strata_2 and strata_3 for the three BMI groups
# df_matched$dataset_stata <- case_when(
#   df_matched$dataset_stata == "normal_MetaCardis_2020_a-1726samples" ~ "strata_1",
#   df_matched$dataset_stata == "overweight_MetaCardis_2020_a-1726samples" ~ "strata_2",
#   df_matched$dataset_stata == "obese_MetaCardis_2020_a-1726samples" ~ "strata_3",
#   TRUE ~ df_matched$dataset_stata
# )

# # Organize columns
# df_matched <- df_matched %>% 
#   select(dataset_stata, subclass, distance, bmi_group, SampleID, Study, Dominant_Species, number_reads, BMI, gender, antibiotics_current_use, study_condition, country) %>% 
#   arrange(dataset_stata, subclass)

# # Save the matched dataset for MetaCardis_2020-1726samples
# write.table(df_matched, "matched_samples_MetaCardis_2020_1726samples_for_Alex_07_Apr_2026.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
```
