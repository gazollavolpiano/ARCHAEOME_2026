#  Additional analysis to address collaborators comments

The following comments were made by collaborators on the manuscript draft:
MI: "Not essential but did you look at whether archeal detection showed signs of plateau’ing with sequencing depth? e.g. it might be useful to give a sense of how many of archea you can detect every additional million reads sequenced etc."
KP: "You could define once that there is a detection limit for archaea and you only detect archaea in those samples where the abundance is higher than the LOD."

## Archaeal detection/abundance vs. sequencing depth in the cohorts

Use the file `archaea_detection_per_sample_with_read_depth_and_age_category.csv` created with `03_age_depth_vs_archaea_detection_V6.md`

```R
# load packages
library(tidyverse)
library(lme4)
library(splines)
library(ggrepel)

# read the data
df <- read.csv("archaea_detection_per_sample_with_read_depth_and_age_category.csv")
dim(df)  # 13666 6

df %>% filter(Study %in% c("FR02")) %>% pull(ReadDepth) %>% summary()
df %>% filter(Study %in% c("AsnicarF_2021-1088samples")) %>% pull(ReadDepth) %>% summary()


df <- df %>% mutate(age_category = factor(age_category, levels = c("newborn","child","schoolage","adult","senior"), labels = c("Infants","Children","Adolescents","Adults","Seniors")))

# pre-scale depth so predict() is trivial and OR rescaling is clean
mu <- mean(df$log10_depth); s <- sd(df$log10_depth)
df <- df %>% mutate(z = (log10_depth - mu)/s)

# ============================================================================
# 1. Model (all data): depth + age + cohort random intercept
# ============================================================================
m_lin2 <- glmer(Detection ~ z + age_category + (1|Study), df, binomial, control = glmerControl(optimizer = "bobyqa"))

# depth OR per 10x reads (z is standardised log10 depth -> divide by s = per 1 log10 unit = per 10x)
b_z    <- fixef(m_lin2)[["z"]]
or_10x <- exp(b_z / s)
ci_z   <- confint(m_lin2, parm = "z", method = "profile")   # use method="profile" for the final version
ci_10x <- exp(ci_z / s)

cat(sprintf("Depth OR per 10x = %.2f (95%% CI %.2f-%.2f)\n", or_10x, ci_10x[1], ci_10x[2]))
# Depth OR per 10x = 2.75 (95% CI 2.32-3.27) ---> same data is shown in Table S2!

# ============================================================================
# 2. Figure data: one point per cohort x age stratum (>=20 samples)
# ============================================================================

strata <- df %>%
  group_by(Study, age_category) %>%
  summarise(n = n(), pos = sum(Detection),
            depthM = median(ReadDepth)/1e6, prev = 100*pos/n,
            lo = 100*binom.test(pos, n)$conf.int[1],
            hi = 100*binom.test(pos, n)$conf.int[2], .groups = "drop") %>%
  filter(n >= 20) %>%
  mutate(Study = str_remove(Study, "-.*"))

# per-age model curve, over each panel's own depth range (free_x)
grid <- strata %>%
  group_by(age_category) %>%
  summarise(depthM = list(10^seq(log10(min(depthM)*0.7),
                                 log10(max(depthM)*1.4), length.out = 200)),
            .groups = "drop") %>%
  unnest(depthM) %>%
  mutate(z = (log10(depthM*1e6) - mu)/s)
grid$fit <- 100*predict(m_lin2, newdata = grid, re.form = NA, type = "response")

pdf("depth_vs_archaea_detection_curves.pdf", width = 7, height = 7)
ggplot() +
  geom_line(data = grid, aes(depthM, fit), colour = "#2b6cb0", linewidth = 1) +
  geom_errorbar(data = strata, aes(depthM, ymin = lo, ymax = hi),
                width = .03, colour = "grey55") +
  geom_point(data = strata, aes(depthM, prev, size = n),
             shape = 21, fill = "#5555f6", colour = "grey20") +
  geom_text_repel(data = strata, aes(depthM, prev, label = Study),
                  size = 2.3, max.overlaps = Inf, seed = 1) +
  facet_wrap(~ age_category, scales = "free_x") +
  scale_x_log10() +
  labs(x = "Cohort median sequencing depth (million reads)", y = "Archaeal detection prevalence (%)", size = "Samples") +
  theme_bw(base_size = 11)
dev.off()
```

## Rarefaction of archaeal detection at different sequencing depths 

Again use the file `archaea_detection_per_sample_with_read_depth_and_age_category.csv` created with `03_age_depth_vs_archaea_detection_V6.md`

```R
# load packages
library(tidyverse)
library(lme4)
library(splines)
library(ggrepel)

# read the data
df <- read.csv("archaea_detection_per_sample_with_read_depth_and_age_category.csv")
dim(df)  # 13666 6

# ============================================================================
# 1. Decide on cohorts and design for rarefaction analysis
# ============================================================================

# what are the cohorts with the deepest sequencing?
df %>%
group_by(Study) %>%
  summarise(depthM = median(ReadDepth)/1e6, n = n()) %>%
  arrange(desc(depthM)) %>%
  print(n=Inf)  
#    Study                         depthM     n
#  1 JieZ_2017-385samples          27.5     385
#  2 AsnicarF_2021-1088samples     26.8    1088
#  3 BrooksB_2017-30samples        24.9      30
#  4 YachidaS_2019-609samples      22.6     609
#  5 QinN_2014-237samples          21.8     237
#  6 BackhedF_2015-200samples      21.6     200
#  7 LiJ_2017-196samples           21.4     196
#  8 MetaCardis_2020_a-1726samples 21.0    1726
#  9 HMP_2019_ibdmdb-130samples    17.3     130
# 10 VatanenT_2016-212samples      14.9     212
# 11 SchirmerM_2016-471samples     14.4     471
# 12 KarlssonFH_2013-145samples    13.9     145
# 13 MehtaRS_2018-308samples       13.2     308
# 14 ShaoY_2019-753samples         10.1     753
# 15 FR02                           0.765  7176 ---> missing 50 because of age missing data

# what are the deepest sequenced samples with archaea detected?
df %>%
  filter(Detection == 1) %>%
  arrange(desc(ReadDepth)) %>%
  mutate(depthM = ReadDepth/1e6) %>%
  select(Study, SampleID, depthM, age_category) %>%
  head(10)
#                         Study   SampleID    depthM age_category
# 1        QinN_2014-237samples      LD.11 109.95238        adult
# 2        QinN_2014-237samples      LD.97  93.91126        adult
# 3   AsnicarF_2021-1088samples ERR4341843  91.16163        adult
# 4        QinN_2014-237samples      LD.75  81.19306        adult
# 5        QinN_2014-237samples      LD.45  78.23605        adult
# 6  KarlssonFH_2013-145samples  ERR260134  73.99702       senior
# 7   AsnicarF_2021-1088samples ERR4334183  72.60216        adult
# 8        QinN_2014-237samples      HD.81  72.28762        adult
# 9   AsnicarF_2021-1088samples ERR4335171  71.46161        adult
# 10 KarlssonFH_2013-145samples  ERR260135  69.21443       senior

# The design will be: 18 samples = (QinN + AsnicarF + KarlssonFH) * (3 low + 3 high archaeal carriage defined by quartiles)
# Compute: 18 * 7 depths * 3 seeds ≈ 380 Kraken2/Bracken runs

# ============================================================================
# 2. Decide on samples from QinN (adults only) + AsnicarF (adults only) + KarlssonFH (seniors only)
# ============================================================================

library(phyloseq)

# read new taxonomy data
physeq_list <- readRDS("physeq_list_metagenomes_with_metadata_and_single_sample_per_individual.rds")

# List the phyloseq objects to be tested 
studies <- c("AsnicarF_2021-1088samples", "KarlssonFH_2013-145samples", "QinN_2014-237samples")

Archaea_prevalence <- data.frame()
for (study_name in studies) {
  cat("Processing study:", study_name, "\n")
  
  phyloseq_obj <- physeq_list[[study_name]]

# Relative abundance per sample (%). Guard against empty samples.
  phyloseq_obj <- transform_sample_counts(phyloseq_obj, function(x) {
    s <- sum(x)
    if (s == 0) return(x)
    (x / s) * 100
  })

  psmelt <- phyloseq_obj %>% psmelt()

  # --- age: adults for QinN & AsnicarF, seniors for KarlssonFH ---
  if (study_name == "KarlssonFH_2013-145samples") psmelt$age_category <- "senior"

  # --- summed archaeal abundance per sample (the new bit) ---
  arch_abund <- psmelt %>%
    filter(Domain == "Archaea") %>%
    group_by(SampleID) %>%
    summarise(archaeal_abundance = sum(Abundance, na.rm = TRUE),
              archaeal_richness  = sum(Abundance > 0),   # n archaeal taxa (full depth)
              .groups = "drop")

  # presence/absence
  archaea_samples <- arch_abund %>% filter(archaeal_abundance > 0) %>% pull(SampleID)

  # --- one row per sample ---
  meta <- psmelt %>%
    select(SampleID, age_category) %>%
    distinct() %>%
    mutate(Study = study_name) %>%
    left_join(arch_abund, by = "SampleID") %>%
    mutate(archaeal_abundance = tidyr::replace_na(archaeal_abundance, 0),
           archaeal_richness  = tidyr::replace_na(archaeal_richness, 0),
           Detection = ifelse(SampleID %in% archaea_samples, 1, 0))

  # --- read depth ---
  file_path <- list.files(paste0("/labs/mericlab/curatedMetagenomicData/", strsplit(study_name, "-")[[1]][1]),
                          pattern = "fastq_stats_", full.names = TRUE)
  fastq_stats <- read.table(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE) %>%
    rename(SampleID = V1, ReadDepth = V4) %>%
    mutate(SampleID = gsub(".*\\/(.*)_1.fastq", "\\1", SampleID),
           SampleID = gsub("-", ".", SampleID)) %>%
    filter(SampleID != "file") %>%
    mutate(ReadDepth = as.numeric(ReadDepth)) %>%
    select(SampleID, ReadDepth)

  meta <- meta %>% left_join(fastq_stats, by = "SampleID")
  Archaea_prevalence <- rbind(Archaea_prevalence, meta)
}

candidates <- Archaea_prevalence %>%
  filter(Detection == 1) %>%
  mutate(depthM = ReadDepth / 1e6) %>%
    # keep only the intended age per cohort
  filter((Study == "KarlssonFH_2013-145samples" & age_category == "senior") |
         (Study %in% c("QinN_2014-237samples","AsnicarF_2021-1088samples") & age_category == "adult")) %>%
  filter(depthM >= 20) %>% # depth floor for rarefaction headroom
  group_by(Study) %>%
  mutate(q = ntile(archaeal_abundance, 4)) %>% # quartiles within cohort
  filter(q %in% c(1, 4)) %>% # Q1 = low carriage, Q4 = high
  mutate(carriage_group = ifelse(q == 1, "low", "high")) %>%
  group_by(Study, carriage_group) %>%
  slice_max(depthM, n = 3) %>% # 3 deepest per group
  ungroup() %>%
  select(Study, SampleID, age_category, depthM, archaeal_abundance, archaeal_richness, carriage_group)

as.data.frame(candidates)   #  18 samples (3 cohorts x 2 groups x 3)

# ============================================================================
# 3. Build a manifest (file + a per-sample depth ladder)
# ============================================================================
BASE <- "/labs/workspace/users/camilagv/Archaeal_Colonization/rarefaction"
RAW  <- file.path(BASE, "raw")                 # where you'll copy the fastqs
dir.create(RAW, recursive = TRUE, showWarnings = FALSE)

ladder_base <- c(0.5, 1, 2, 5, 10, 20, 30, 50)  # million read-pairs (80 dropped: ~full file)

manifest <- candidates %>%
  mutate(
    cohort   = sub("-.*", "", Study),
    base     = gsub("\\.", "-", SampleID),      # HD.52 -> HD-52 ; ERR ids unchanged
    r1       = file.path(RAW, paste0(base, "_1.fastq")),
    r2       = file.path(RAW, paste0(base, "_2.fastq")),
    # rungs strictly below this sample's max; full-depth file handled separately
    rungs    = map_chr(depthM, ~ paste(ladder_base[ladder_base < .x], collapse = ","))) %>%
  transmute(sample = SampleID, cohort, carriage_group, rel_abundance_pct = archaeal_abundance, max_depthM = round(depthM, 1), r1, r2, rungs)

write_tsv(manifest, file.path(BASE, "manifest.tsv"))
```

Now, let's copy the files to the folder and subsample.

```bash
# go to the folder 
cd rarefaction

# copy files
BASE=/labs/workspace/users/camilagv/Archaeal_Colonization/rarefaction
SRC=/labs/mericlab/curatedMetagenomicData
RAW=$BASE/raw

tail -n +2 "$BASE/manifest.tsv" | while IFS=$'\t' read -r sample cohort grp ab maxd r1 r2 rungs; do
  for dest in "$r1" "$r2"; do
    b=$(basename "$dest" .fastq)
    echo "Processing $b"
    src=$SRC/$cohort/fastq/${b}.fastq
    [ -s "$dest" ] && { echo "skip $b"; continue; }
    if   [ -r "$src" ];    then cp "$src" "$dest";              echo "cp   $b"
    elif [ -r "$src.xz" ]; then xz -dc "$src.xz" > "$dest";     echo "unxz $b"
    elif [ -r "$src.gz" ]; then gzip -dc "$src.gz" > "$dest";   echo "gunz $b"
    else echo "MISSING $b" >&2
    fi
  done
done

ls "$RAW" | wc -l  # 36 (18 samples x 2 reads)
```

Add the following script `rarefy_rasusa.sh` to the folder `rarefaction`:

```bash
#!/bin/bash
set -euo pipefail

BASE=/labs/workspace/users/camilagv/Archaeal_Colonization/rarefaction
MANIFEST=$BASE/manifest.tsv
mkdir -p "$BASE/subsampled" "$BASE/logs"

# this array task's manifest row (line 1 = header)
LINE=$(( SLURM_ARRAY_TASK_ID + 1 ))
IFS=$'\t' read -r sample cohort carriage ab maxd r1 r2 rungs < <(sed -n "${LINE}p" "$MANIFEST")

echo "=== task ${SLURM_ARRAY_TASK_ID} | $sample ($cohort, $carriage, ${ab}%) | max ${maxd}M ==="
[ -s "$r1" ] && [ -s "$r2" ] || { echo "FATAL: missing fastq for $sample"; exit 1; }

# 3 seeds for low carriage (near the detection limit, noisy); 2 for high (stable)
if [ "$carriage" = "low" ]; then SEEDS="100 200 300"; else SEEDS="100 200"; fi

IFS=',' read -ra RUNGS <<< "$rungs"
for d in "${RUNGS[@]}"; do
  N=$(python3 -c "print(int(round($d*1e6)))")
  for s in $SEEDS; do
    o1=$BASE/subsampled/${sample}__d${d}M__s${s}_1.fastq.gz
    o2=$BASE/subsampled/${sample}__d${d}M__s${s}_2.fastq.gz
    { [ -s "$o1" ] && [ -s "$o2" ]; } && { echo "  ${d}M s${s} -> exists, skip"; continue; }

    t0=$SECONDS
    rasusa reads --num "$N" --seed "$s" -o "$o1" -o "$o2" "$r1" "$r2"
    g1=$(zcat "$o1" | wc -l | awk '{print $1/4}')
    g2=$(zcat "$o2" | wc -l | awk '{print $1/4}')
    [ "$g1" = "$g2" ] || { echo "FATAL: pair mismatch $g1 vs $g2"; exit 1; }
    echo "  ${d}M s${s} -> $g1 pairs ($(ls -lh "$o1" | awk '{print $5}')) time=$((SECONDS-t0))s"
  done
done
echo "DONE $sample"
```

```bash
# install rasusa in a micromamba environment
micromamba create -n rasusa -c bioconda rasusa
micromamba activate rasusa
rasusa --version # rasusa 4.1.0

# one task first
sbatch --job-name=rarefy --array=1 \
  --output=$PWD/logs/rarefy_%A_%a.log --error=$PWD/logs/rarefy_%A_%a.err \
  --cpus-per-task=2 --mem=8G --time=12:00:00 --partition=standard \
  rarefy_rasusa.sh 
  
# Submitted batch job 2336034
sacct -j 2336034 --format=JobID,MaxRSS,Elapsed,AllocCPUS # trim cores, raise memory?

# remaining tasks now
wc -l manifest.tsv # 19 (header + 18 samples)

sbatch --job-name=rarefy --array=1-18%18 \
  --output=$PWD/logs/rarefy_%A_%a.log --error=$PWD/logs/rarefy_%A_%a.err \
  --cpus-per-task=2 --mem=4G --time=12:00:00 --partition=standard \
  rarefy_rasusa.sh

# Submitted batch job 2336064
tail logs/*err -n 1

# JOBS CANCELLED DUE TO TIME LIMIT: 4, 5, 6, 10, 16, 17, 18

# resume-safe, requeue the 6 that died with more time 
sbatch --job-name=rarefy --array=4,5,6,10,16,17,18%6 \
  --output=$PWD/logs/rarefy_%A_%a.log --error=$PWD/logs/rarefy_%A_%a.err \
  --cpus-per-task=2 --mem=4G --time=24:00:00 --partition=standard \
  rarefy_rasusa.sh

# Submitted batch job 2336238
tail logs/*2336238*err -n 1 # eveything is fine now
ls subsampled/ | wc -l # 664

BASE=/labs/workspace/users/camilagv/Archaeal_Colonization/rarefaction
tail -n +2 $BASE/manifest.tsv | while IFS=$'\t' read -r s c grp ab maxd r1 r2 rungs; do
  nr=$(( $(echo "$rungs" | tr ',' '\n' | wc -l) ))
  ns=$([ "$grp" = low ] && echo 3 || echo 2)
  exp=$(( nr * ns ))
  got=$(ls $BASE/subsampled/${s}__*_1.fastq.gz 2>/dev/null | wc -l)
  printf "%-12s %-5s rungs=%s seeds=%s  expect=%s got=%s %s\n" \
    "$s" "$grp" "$nr" "$ns" "$exp" "$got" "$([ "$exp" = "$got" ] && echo OK || echo SHORT)"
done
# everything is ok
```

Now, we use the Kraken2/Bracken pipeline. Add `step1_slurm_array.sh` to the folder:

```bash
#!/bin/bash
#step1_slurm_array.sh 
#SBATCH -p long
#SBATCH --time=00-04:00:00                       
#SBATCH --mem=50G # Total memory of 50G per task
#SBATCH --cpus-per-task=20
#SBATCH --array=1-332%20

BASE=/labs/workspace/users/camilagv/Archaeal_Colonization/rarefaction
DB=/labs/mericlab/databases/k2_HPRC_20230810/db
SUB=$BASE/subsampled
OUT=$BASE/temporary_files; mkdir -p $OUT

# Pull one single line/entry from the list
tok=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $BASE/sample_names.txt)
r1=$SUB/${tok}_1.fastq.gz
r2=$SUB/${tok}_2.fastq.gz
res=$OUT/${tok}.HPRC.kraken
[ -s "$res" ] && { echo "skip $tok"; exit 0; }

# Initialize shell and activate the micromamba environment
eval "$(micromamba shell hook --shell=bash)"
micromamba activate metagenome_taxonomy 

# Run the appropriate command (samples are paired)
kraken2 --db $DB --threads 20 --paired "$r1" "$r2" --output $res --gzip-compressed

echo "DONE $tok"
```

```bash
# create sample_names.txt 
BASE=/labs/workspace/users/camilagv/Archaeal_Colonization/rarefaction
ls subsampled/*_1.fastq.gz | sed 's#.*/##; s/_1\.fastq\.gz$//' > sample_names.txt
wc -l sample_names.txt # 332 (664/2)

head sample_names.txt -n 3
# ERR260134__d0.5M__s100
# ERR260134__d0.5M__s200
# ERR260134__d0.5M__s300

# send slurm array script to cluster, test one sample
sbatch step1_slurm_array.sh # Submitted batch job 2338362

# count number of files processed
ls temporary_files/*.HPRC.kraken | wc -l # 332

# move files to misc_slurm
mkdir misc_slurm
mv slurm-2338362_*.out misc_slurm
```

Now, use the following script to extract the non-human reads from the samples:

```bash
#!/bin/bash
#step2_slurm_array.sh 
#SBATCH -p long
#SBATCH --mem=15G                  # Use 15G 
#SBATCH --cpus-per-task=1          # Use 1 CPU
#SBATCH --array=1-332%20

BASE=/labs/workspace/users/camilagv/Archaeal_Colonization/rarefaction
SUB=$BASE/subsampled
OUT=$BASE/temporary_files

# Pull one single line/entry from the list
tok=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $BASE/sample_names.txt)
res=$OUT/${tok}.HPRC.kraken
[ -s "$res" ] || { echo "FATAL: missing $res"; exit 1; }      # need the input
[ -s "$r1_noh" ] && { echo "skip $tok"; exit 0; }             # output done -> resume-safe

# initialize your shell and activate the micromamba environment
eval "$(micromamba shell hook --shell=bash)"
micromamba activate metagenome_taxonomy 

# Define input and output file paths
r1=$SUB/${tok}_1.fastq.gz;  r2=$SUB/${tok}_2.fastq.gz
r1_noh=$OUT/nonhuman_${tok}.1.fastq;  r2_noh=$OUT/nonhuman_${tok}.2.fastq

# Check if the sample is paired or single and run the appropriate command
/home/cgazollavolpiano/micromamba/envs/metagenome_taxonomy/bin/python /labs/mericlab/databases/kraken_utils/extract_kraken_reads.py \
    -k "$res" -s1 "$r1" -s2 "$r2" -t 0 -o "$r1_noh" -o2 "$r2_noh" --fastq-output
``` 

Now run:

```bash
# send slurm array script to cluster
sbatch --time=00-10:00:00 step2_slurm_array.sh # batch job 2338695

# count number of files processed
ls temporary_files/*.1.fastq | wc -l # 332
ls temporary_files/*.2.fastq | wc -l # 332

# move files to misc_slurm
mv slurm-2338695_*.out misc_slurm
```

Now add and run the script for classification with gtdb r214 indexes:

```bash
#!/bin/bash
# step3_slurm_array.sh 
#SBATCH -p standard
#SBATCH --time=01-00:00:00
#SBATCH --mem=132G          # Total memory of 132G per task
#SBATCH --cpus-per-task=4

BASE=/labs/workspace/users/camilagv/Archaeal_Colonization/rarefaction
DB=/labs/mericlab/databases/kraken2_indices/gtdb_r214/128gb
TEMP=$BASE/temporary_files
OUT=$BASE/gtdb_results; mkdir -p $OUT

# Pull one single line/entry from the list
tok=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $BASE/sample_names.txt)

# Initialize your shell and activate the micromamba environment
eval "$(micromamba shell hook --shell=bash)"
micromamba activate metagenome_taxonomy 

# Define input and output file paths for Kraken2 results
r1_noh=$TEMP/nonhuman_${tok}.1.fastq;  r2_noh=$TEMP/nonhuman_${tok}.2.fastq

kraken_output=$OUT/${tok}.128gb_kraken_gtdb_output.txt
kraken_report=$OUT/${tok}.128gb_kraken_gtdb_report.txt

[ -s "$kraken_report" ] && { echo "skip $tok"; exit 0; }
[ -s "$r1_noh" ] || { echo "FATAL: missing $r1_noh"; exit 1; }

# Run Kraken2 based on sample type
kraken2 --db "$DB" --threads 4 --output "$kraken_output" --report "$kraken_report" --paired "$r1_noh" "$r2_noh" --confidence 0.1 --use-names
``` 

And then run with:

```bash
# using slurm array to run the script
sbatch --array=1-332%10 step3_slurm_array.sh # Submitted batch job 2341132

# count how many files are done
ls gtdb_results/*.128gb_kraken_gtdb_report.txt | wc -l # 332
mv slurm-2341132_*.out misc_slurm

# Bracken to estimate the abundance for each taxon
micromamba activate metagenome_taxonomy 

AsnicarF_2021 # bracken run with -r 150
KarlssonFH_2013 # bracken run with -r 100
QinN_2014 # bracken run with -r 100

tail -n +2 manifest.tsv | while IFS=$'\t' read -r sample cohort rest; do 
    # 1. Determine read length based on cohort
    read_length=150 
    if [[ "$cohort" == "KarlssonFH_2013" ]] || [[ "$cohort" == "QinN_2014" ]]; then
        read_length=100
    fi
    echo "Checking files for sample: $sample ($cohort) using read length: $read_length"
    # 2. Loop through all matching Kraken output files for this sample
    for infile in gtdb_results/${sample}*.128gb_kraken_gtdb_report.txt; do
        # If no files match the pattern, skip to the next sample
        [[ -e "$infile" ]] || { echo "  -> No matching files found for $sample"; continue; }
        # 3. Dynamically generate the output filename
        outfile="${infile%.128gb_kraken_gtdb_report.txt}.128gb_bracken_25_S_gtdb_report.txt"
        echo "  -> Running Bracken on: $(basename "$infile")"
        bracken -d /labs/mericlab/databases/kraken2_indices/gtdb_r214/128gb -i "$infile" -o "$outfile" -r "$read_length" -l S -t 25
    done
done

# count how many files are done
ls gtdb_results/*.128gb_bracken_25_S_gtdb_report.txt | wc -l # 332

# clean unnecessary files
rm temporary_files gtdb_results -r

# combine outputs 
python /labs/mericlab/databases/kraken_utils/combine_bracken_outputs_mod.py --files gtdb_results/*128gb_bracken_25_S_gtdb_report.txt -o bracken_128gb_25_S_gtdb_report_merged.txt

# start R
micromamba deactivate
micromamba activate R && R
```

```R
# library
library(phyloseq)
library(tidyverse)

# increase the stdout limit
options(width = 200)

# read merged bracken output result and modify it
merged_bracken_output <- read.table("bracken_128gb_25_S_gtdb_report_merged.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
                          select(-matches("frac|taxonomy_id|taxonomy_lvl")) %>% 
                          rename(Species = name) 
merged_bracken_output <- merged_bracken_output %>% rename_all(~ str_remove(., ".128gb_bracken_25_S_gtdb_report.txt_num"))
head(merged_bracken_output)
dim(merged_bracken_output) # 5221x333 samples (1 col is "Species")

# read the taxonomy auxiliary file
taxonomy <- read.table("/labs/mericlab/databases/kraken_utils/gtdb_214_taxnames/gtdb_v214_taxonomy.txt", 
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# add the full taxonomy to the merged bracken output 
merged_bracken_output <- merge(taxonomy, merged_bracken_output, all.y = TRUE)

# create otu table and tax table
otu <- merged_bracken_output %>% 
select(-Domain, -Phylum, -Class, -Order, -Family, -Genus, Species) %>%
column_to_rownames("Species")

tax <- merged_bracken_output %>% 
select(Domain, Phylum, Class, Order, Family, Genus, Species) 
rownames(tax) <- tax$Species

# read file with information about the samples from the studies
sample_info <- read.table("manifest.tsv", header = TRUE) %>% 
  select(sample, cohort, carriage_group, rel_abundance_pct, max_depthM, rungs)

# ensure sample names are the row names of the metadata
sample_raref <- data.frame(sampleID = colnames(otu), sample = gsub("__.*","",colnames(otu)))
sample_info <- left_join(sample_raref, sample_info)
rownames(sample_info) <- sample_info$sampleID
sam <- sample_data(sample_info)

# create complete phyloseq object 
physeq <- phyloseq(
  otu_table(as.matrix(otu), taxa_are_rows = TRUE), 
  tax_table(as.matrix(tax)),
  sam
)

# Double check that everything merged correctly
physeq
# otu_table()   OTU Table:         [ 5221 taxa and 332 samples ]
# sample_data() Sample Data:       [ 332 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 5221 taxa by 7 taxonomic ranks ]

# save 
saveRDS(physeq, "physeq_raref_metagenomes.rds")
```

Now we will use `physeq` object with the rarefaction results and the old `physeq_list_metagenomes_with_metadata_and_single_sample_per_individual.rds`, not rarefied.

```R
# request libs
library(phyloseq); library(tidyverse); library(lme4); library(splines)

# read files
physeq      <- readRDS("physeq_raref_metagenomes.rds")
physeq_list <- readRDS("/labs/workspace/users/camilagv/Archaeal_Colonization/physeq_list_metagenomes_with_metadata_and_single_sample_per_individual.rds")
studies     <- c("AsnicarF_2021-1088samples","KarlssonFH_2013-145samples","QinN_2014-237samples")

# helper: archaeal taxa (robust to "Archaea" vs "d__Archaea")
arch_rows <- function(ps) rownames(tax_table(ps))[grepl("Archaea", as(tax_table(ps),"matrix")[,"Domain"])]

# ============================================================================
# Detection at each rarefied depth (sample x depth x seed)
# ============================================================================
otu  <- as(otu_table(physeq), "matrix"); if (!taxa_are_rows(physeq)) otu <- t(otu)
arch <- arch_rows(physeq)

raref <- tibble(
  sampleID   = colnames(otu),
  arch_reads = colSums(otu[arch, , drop = FALSE]),
  tot_reads  = colSums(otu)) %>%
  mutate(
    sample = sub("__.*", "", sampleID),
    depth  = as.numeric(sub(".*__d([0-9.]+)M__.*", "\\1", sampleID)),   # nominal M read-pairs
    seed   = sub(".*__s", "", sampleID),
    detected = as.integer(arch_reads > 0)) %>%
  left_join(data.frame(sample_data(physeq)) %>%
              distinct(sample, cohort, carriage_group, rel_abundance_pct, max_depthM),
            by = "sample")

tail(as.data.frame(raref),2)
#             sampleID arch_reads tot_reads sample depth seed detected    cohort
# 331 LD.83__d5M__s100        794   4210485  LD.83     5  100        1 QinN_2014
# 332 LD.83__d5M__s200        816   4209964  LD.83     5  200        1 QinN_2014
#     carriage_group rel_abundance_pct max_depthM
# 331           high        0.01906778       34.4
# 332           high        0.01906778       34.4

otu %>% 
as.data.frame() %>%
filter(LD.83__d5M__s100 !=0) %>%
select(LD.83__d5M__s100) %>%
rownames_to_column()  %>%
filter(rowname %in% arch)
#             rowname LD.83__d5M__s100
# 1 UBA71 sp905187815              794

otu %>% 
as.data.frame() %>%
summarise(tot_reads=sum(LD.83__d5M__s100))
#   tot_reads
# 1   4210485

# ============================================================================
# Full-depth anchor point per sample (caps each curve at its true max)
# ============================================================================
anchor <- map_dfr(studies, function(st){
  ps <- physeq_list[[st]]; o <- as(otu_table(ps),"matrix"); if(!taxa_are_rows(ps)) o <- t(o)
  a  <- arch_rows(ps)
  tibble(sample = colnames(o), arch_reads = colSums(o[a,,drop=FALSE]), tot_reads = colSums(o))
}) %>%
  inner_join(raref %>% distinct(sample, cohort, carriage_group, rel_abundance_pct, max_depthM),
             by = "sample") %>%
  transmute(sampleID = paste0(sample,"__full"), arch_reads, tot_reads,
            sample, depth = max_depthM, seed = "full",
            detected = as.integer(arch_reads > 0),
            cohort, carriage_group, rel_abundance_pct, max_depthM)

head(as.data.frame(anchor),1)
#           sampleID arch_reads tot_reads     sample depth seed detected
# 1 ERR4334123__full        630  60789472 ERR4334123  62.9 full        1
#          cohort carriage_group rel_abundance_pct max_depthM
# 1 AsnicarF_2021            low       0.001036364       62.9

det <- bind_rows(raref, anchor) 

det %>% filter(sample == "ERR4334123") %>% arrange(depth) %>% as.data.frame() # ok

# Save 
write.csv(det %>% arrange(cohort, sample, depth), file="raref_result_final.csv")

# ============================================================================
# PLOT A, aggregate detection vs depth, by carriage group (rarefied only, so points share a common x-grid)
# ============================================================================
agg <- det %>% filter(seed != "full") %>%
  group_by(carriage_group, depth) %>%
  summarise(n = n(), k = sum(detected), p = 100*k/n,
            lo = 100*binom.test(k,n)$conf.int[1],
            hi = 100*binom.test(k,n)$conf.int[2], .groups = "drop")


pA <- ggplot(agg, aes(depth, p, colour = carriage_group, fill = carriage_group)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .15, colour = NA) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_x_log10(breaks = c(0.5,1,2,5,10,20,50)) +
  scale_colour_manual(values = c(low = "#dd6b20", high = "#2b6cb0"), aesthetics = c("colour","fill")) +
  labs(x = "Sequencing depth (million read pairs)",
       y = "Subsampled libraries with archaea detected (%)", colour = "Carriage", fill = "Carriage") +
  theme_bw(base_size = 11)

ggsave("rarefaction_detection_aggregate.pdf", pA, width = 6, height = 4.2)

# ============================================================================
# EMPIRICAL LOD, depth needed for reliable detection vs archaeal abundance
# ============================================================================
lod <- det %>% filter(seed != "full") %>%
  group_by(sample, cohort, carriage_group, rel_abundance_pct, depth) %>%
  summarise(all_seeds = all(detected == 1), .groups = "drop") %>%
  group_by(sample, cohort, carriage_group, rel_abundance_pct) %>%
  summarise(lod_depth = { d <- depth[all_seeds]; if (length(d)) min(d) else NA_real_ },
            .groups = "drop")

# ---- fit the LOD-abundance relationship (log-log) ----
lod_fit <- lod %>% filter(is.finite(lod_depth), rel_abundance_pct > 0)
fit <- lm(log10(lod_depth) ~ log10(rel_abundance_pct), data = lod_fit)
b     <- coef(fit)[["log10(rel_abundance_pct)"]]      # slope (per decade of abundance)
ci    <- confint(fit)["log10(rel_abundance_pct)", ]
r2    <- summary(fit)$r.squared

newx  <- tibble(rel_abundance_pct = 10^seq(log10(min(lod_fit$rel_abundance_pct)), log10(max(lod_fit$rel_abundance_pct)), length.out = 100))
pred  <- bind_cols(newx, as_tibble(predict(fit, newx, interval = "confidence"))) %>% mutate(across(c(fit, lwr, upr), ~ 10^.x)) # back-transform to raw depth
pv <- summary(fit)$coefficients["log10(rel_abundance_pct)", "Pr(>|t|)"]

lab <- if (pv < 0.001) sprintf("atop(R^2 == %.2f, italic(P) < 0.001)", r2) else sprintf("atop(R^2 == %.2f, italic(P) == %.3f)", r2, pv)

pLOD <- ggplot(lod_fit, aes(rel_abundance_pct, lod_depth)) +
  geom_ribbon(data = pred, aes(y = fit, ymin = lwr, ymax = upr),
              fill = "grey70", alpha = .3) +
  geom_line(data = pred, aes(y = fit), colour = "grey25", linewidth = 1) +
  geom_point(aes(colour = carriage_group), size = 2.4) +
  scale_x_log10() + scale_y_log10() +
scale_colour_manual(
  values = c(low = "#dd6b20", high = "#2b6cb0"),
  breaks = c("high", "low"),
  labels = c("High (median 0.62%)", "Low (median 0.0005%)"),
  name   = "Archaeal carriage")+
    annotate("text", x = min(lod_fit$rel_abundance_pct), y = max(lod_fit$lod_depth),
           hjust = -3, vjust = 1, size = 3.3, parse = TRUE, label = lab) +
  labs(x = "Archaeal relative abundance at full depth (%)",
       y = "Minimum depth for detection in all seeds (M read pairs)", colour = "Carriage") +
  theme_bw(base_size = 11)

ggsave("rarefaction_LOD_vs_abundance.pdf", pLOD, width = 5.5, height = 5.5)

confint(fit)
#                               2.5 %     97.5 %
# (Intercept)              -0.5717989 -0.1598929
# log10(rel_abundance_pct) -0.5864233 -0.4131116

# ============================================================================
# Some numbers for text....
# ============================================================================

lod_fit %>%
  group_by(carriage_group) %>%
  summarise(n       = n(),
            ab_med  = median(rel_abundance_pct),
            ab_min  = min(rel_abundance_pct),
            ab_max  = max(rel_abundance_pct),
            lod_med = median(lod_depth),
            lod_min = min(lod_depth),
            lod_max = max(lod_depth))

#   carriage_group     n   ab_med    ab_min  ab_max lod_med lod_min lod_max
# 1 high               9 0.615    0.0191    2.70        0.5     0.5       1
# 2 low                8 0.000532 0.0000858 0.00204    20      10        50
```