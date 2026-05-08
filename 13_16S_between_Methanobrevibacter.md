# 16S rRNA difference between M. smithii, M. intestini, and M. sp900766745

We will use barrnap to find rRNA genes and write their coordinates in GFF3. It detects 5S, 16S, and 23S rRNA features from a genome FASTA. 

bedtools getfasta can then extract the sequence for those coordinates directly from the original FASTA, including strand-aware extraction with -s.

```bash
# download the GTDB r226 metadata for archaea and decompress it:
wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/ar53_metadata_r226.tsv.gz
gunzip ar53_metadata_r226.tsv.gz
```

Use R to select the genomes for the major archaeal species found in the human gut:

```R
# library
library(tidyverse)

# read the metadata file
metadata <- read.csv("ar53_metadata_r226.tsv", sep = "\t", stringsAsFactors = FALSE) %>% select(accession, gtdb_taxonomy) 

# add species column from gtdb_taxonomy 
metadata$species <- gsub(".*;s__", "", metadata$gtdb_taxonomy)

# list of major archaeal species found in the human gut
major_taxa <- c("Methanocatella smithii", "Methanocatella smithii_A", "Methanocatella sp900766745")
metadata <- metadata %>% filter(species %in% major_taxa)

#check the number of genomes for each species
table(metadata$species)
#     Methanocatella smithii   Methanocatella smithii_A 
#                        107                         29 
# Methanocatella sp900766745 
#                          6 

# save metadata for the selected genomes
write.csv(metadata, "ar53_r226_metadata_for_142_selected_genomes.csv", row.names = FALSE)

# save simplied acession numbers for the selected genomes 
write.table(gsub(".*GC.(_...........)", "GCA\\1", metadata$accession), file = "accession_numbers_for_142_selected_genomes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
```

Now, let's use the accession numbers to download the corresponding genome FASTA files, and then run barrnap to find the rRNA genes. 

```bash
# download the genome FASTA files for the 142 selected genomes using the accession numbers
mkdir genome_fastas

while read -r acc; do
  echo "Processing genome: $acc"
  datasets download genome accession "$acc" --filename "genome_fastas/${acc}.zip"
  unzip -o "genome_fastas/${acc}.zip" -d "genome_fastas/${acc}"
  mv genome_fastas/${acc}/ncbi_dataset/data/${acc}/*.fna genome_fastas/${acc}.fna
  rm -rf "genome_fastas/${acc}" "genome_fastas/${acc}.zip"
done < accession_numbers_for_142_selected_genomes.txt

ls genome_fastas/*.fna | wc -l # 127, 15 genomes with no FASTA files

# create folders for the barrnap output and temporary files
mkdir -p barrnap_gff barrnap_16S temp_bed

for f in genome_fastas/*.fna; do
  base=$(basename "$f" .fna)
  echo "Processing $base"

  # 1) find rRNA genes
  barrnap --kingdom arc "$f" > "barrnap_gff/${base}.rrna.gff"

  # 2) keep only 16S and convert to BED with custom names
  awk -v base="$base" 'BEGIN{FS=OFS="\t"; i=0}
    $0 !~ /^#/ && /Name=16S_rRNA/ {
      i++
      print $1, $4-1, $5, base"_16S_copy"i, ".", $7
    }' "barrnap_gff/${base}.rrna.gff" > "temp_bed/${base}.16S.bed"

  # 3) extract the 16S sequence(s)
  if [ -s "temp_bed/${base}.16S.bed" ]; then
    bedtools getfasta \
      -fi "$f" \
      -bed "temp_bed/${base}.16S.bed" \
      -s -name \
      > "barrnap_16S/${base}.16S.fa"
  else
    echo "No 16S found for $base"
  fi
done

# Check how many 16S sequences we have
ls barrnap_16S/*.16S.fa | wc -l # we have 54, from 73 of the 127 genomes have 0 ribosomal RNA features... (M. sp900766745 is represented by only GCA_937975995.1!)

# I will quickly check GCA_937975995.1 with a NCBI BLAST search to see if it is not an obvious contaminant...
cat barrnap_16S/GCA_937975995.1.16S.fa
# Answer = ok, similar to Methanobrevibacter, good sign that it is not a contaminant
```

Now we will:
- combine all Barrnap 16S copies 
- align the full 16S 
- calculate pairwise identity after checking the aligment with MEGA (just to remove any obviously bad sequences that might mess up the distance matrix)

```R
# libraries
library(tidyverse)
library(Biostrings)
library(DECIPHER)

# -----------------------------
# 1) Read metadata and sequences
# -----------------------------
meta <- read.csv("ar53_r226_metadata_for_142_selected_genomes.csv", stringsAsFactors = FALSE)
meta <- meta %>% mutate(acc_clean = gsub("^(RS_|GB_)", "", accession), species_label = species) %>% select(acc_clean, species_label)
meta$acc_simple <- gsub("GC._(.........)..", "\\1", meta$acc_clean)

# Read all barrnap 16S fasta files
fa_files <- list.files("barrnap_16S", pattern = "\\.fa$", full.names = TRUE)
seqs <- do.call(c, lapply(fa_files, readDNAStringSet))

# Parse accession + copy from fasta headers
annot <- tibble(
  old_name = names(seqs),
  acc_clean = sub("^((GCA|GCF)_[0-9]+\\.[0-9]+).*", "\\1", names(seqs)),
  acc_simple = sub("GC._(.........)..", "\\1", acc_clean),
  copy = sub("^.*(_16S_copy[0-9]+).*", "\\1", names(seqs))
) %>%
  left_join(meta %>% select(acc_simple, species_label), by = "acc_simple")

# Create nicer fasta names
annot <- annot %>% mutate(species_tag = gsub("[^A-Za-z0-9]+", "_", species_label), new_name = paste(acc_clean, species_tag, copy, sep = "|"))
names(seqs) <- annot$new_name

# Save combined full-length fasta
writeXStringSet(seqs, "all_16S.full.fa")

# ----------------------------------
# 3) Align with DECIPHER
# ----------------------------------
# Full-length alignment
aln_full <- AlignSeqs(seqs, processors = NULL, verbose = TRUE)
writeXStringSet(aln_full, "all_16S.full.DECIPHER.aln.fa")

# I LOOKED AT THE ALIGNMENT WITH MEGA AND WE NEED TO REMOVE 4 SEQUENCES THAT WERE OBVIOUSLY BAD 
# GCA 000189835.2|Methanocatella smithii| 16S copy1 --> probably a contaminant, very different from the rest of the sequences
# GCA 037250035.1|Methanocatella smithii| 16S copy1 --> too many Ns, especially in the parts we have data for M. sp900766745
# GCA 037389515.1|Methanocatella smithii| 16S copy1 --> too many Ns, especially in the parts we have data for M. sp900766745
# GCA 041537865.1|Methanocatella smithii| 16S copy2 --> probably a contaminant, very different from the rest of the sequences

length(aln_full) # 66 sequences before removing the bad ones

aln_full <- aln_full[!names(aln_full) %in% c(
  "GCA_000189835.2|Methanocatella_smithii|_16S_copy1",
  "GCA_037250035.1|Methanocatella_smithii|_16S_copy1",
  "GCA_037389515.1|Methanocatella_smithii|_16S_copy1",
  "GCA_041537865.1|Methanocatella_smithii|_16S_copy2"
)]

length(aln_full) # 62 sequences after removing the bad ones

# ----------------------------------
# 4) Pairwise identity from DECIPHER
# ----------------------------------
# DistanceMatrix returns pairwise dissimilarities on aligned sequences.
# Convert to percent identity as (1 - distance) * 100.
dist_to_pid <- function(aln) {
  d <- DistanceMatrix(aln)
  pid <- (1 - d) * 100
  diag(pid) <- 100
  pid
}

pid_full <- dist_to_pid(aln_full)

write.csv(as.data.frame(pid_full), "pairwise_identity_full_16S_DECIPHER.csv", row.names = TRUE)

# ----------------------------------
# 5) Summarise within / between species
# ----------------------------------
matrix_to_long <- function(m, region) {
  as.data.frame(as.table(m), stringsAsFactors = FALSE) %>%
    setNames(c("seq1", "seq2", "pid")) %>%
    filter(seq1 < seq2) %>%
    mutate(region = region)
}

long_full <- matrix_to_long(pid_full, "full_16S")

head(as.data.frame(long_full), 2)
#                                                seq1
# 1 GCA_000016525.1|Methanocatella_smithii|_16S_copy1
# 2 GCA_000016525.1|Methanocatella_smithii|_16S_copy1
#                                                seq2       pid   region
# 1 GCA_000016525.1|Methanocatella_smithii|_16S_copy2  98.78706 full_16S
# 2 GCA_000151225.1|Methanocatella_smithii|_16S_copy1 100.00000 full_16S

# extract species from the sequence names (second pipe-delimited field)
long_full$sp1 <- sub(".*\\|(.+)\\|_16S.*", "\\1", long_full$seq1)
long_full$sp2 <- sub(".*\\|(.+)\\|_16S.*", "\\1", long_full$seq2)

# extract accession (first field) to flag within-genome pairs
long_full$acc1 <- sub("\\|.*", "", long_full$seq1)
long_full$acc2 <- sub("\\|.*", "", long_full$seq2)

# sort species labels so sp_a <= sp_b (avoids duplicate A-B / B-A pairs)
long_full$sp_a <- pmin(long_full$sp1, long_full$sp2)
long_full$sp_b <- pmax(long_full$sp1, long_full$sp2)

# remove within-genome copy comparisons (same accession, different 16S copies)
df <- subset(long_full, acc1 != acc2)

# summarise
species_summary_full <- aggregate(pid ~ sp_a + sp_b, data = df, FUN = function(x) {
  c(n      = length(x),
    mean   = round(mean(x),   3),
    median = round(median(x), 3),
    min    = round(min(x),    3),
    max    = round(max(x),    3))
})

as.data.frame(do.call(cbind, list(
  species_summary_full[, c("sp_a","sp_b")],
  as.data.frame(species_summary_full$pid)
)))

#                       sp_a                       sp_b    n   mean median    min
# 1   Methanocatella_smithii     Methanocatella_smithii 1476 99.779 99.932 96.346
# 2   Methanocatella_smithii   Methanocatella_smithii_A  330 99.542 99.593 98.315
# 3 Methanocatella_smithii_A   Methanocatella_smithii_A   14 99.898 99.864 99.796
# 4   Methanocatella_smithii Methanocatella_sp900766745   55 97.077 97.102 95.914
# 5 Methanocatella_smithii_A Methanocatella_sp900766745    6 96.977 96.999 96.898
#       max
# 1 100.000
# 2  99.813
# 3 100.000
# 4  97.772
# 5  97.033
```

Considering that microbiome analyses often rely on 16S amplicons covering only a fraction of the gene, the resolution is likely to be even worse in practice... It is well known that 16S can fail to resolve closely related species, and considering that these 3 are very close, I don't expect that using 16S V4 regions will to do a good job. 

I should try to cut the aligments to get V4 with 515F (5’-GTGYCAGCMGCCGCGGTAA-3’) and 806R (5’-GGACTACNVGGGTWTCTAAT-3’) primers. 

```bash
# cut the full-length 16S alignment to the V4 region using cutadapt with anchored primers and strict parameters to avoid indels and mismatches in the primer regions
cutadapt  -g GTGYCAGCMGCCGCGGTAA...ATTAGAWACCCBNGTAGTCC --discard-untrimmed --no-indels -e 3 -m 150 -M 400 -o all_16S.V4.cutadapt.fa all_16S.full.fa     
# discard sequences where this primer pair is not found in the required arrangement
# do not allow insertions/deletions when matching the primers; mismatches only
# allow up to 3 errors in each primer match
# keep only trimmed sequences with length >= 150 nt, length <= 400 nt

head all_16S.V4.cutadapt.fa # seems good
grep ">" all_16S.V4.cutadapt.fa | wc -l # 61 sequences, 5 sequences were discarded because the primers were not found in the required arrangement

# I checked with MEGA and we have two weird 16S copies that we need to remove, but the rest of the sequences look good, with no obvious misalignments or bad sequences
# GCA 000189835.2|Methanocatella smithii| 16S copy1
# GCA 041537865.1|Methanocatella smithii| 16S copy2
```

```R
# libraries
library(tidyverse)
library(Biostrings)
library(DECIPHER)

# -----------------------------
# 1) Read metadata and sequences
# -----------------------------
meta <- read.csv("ar53_r226_metadata_for_142_selected_genomes.csv", stringsAsFactors = FALSE)
meta <- meta %>% mutate(acc_clean = gsub("^(RS_|GB_)", "", accession), species_label = species) %>% select(acc_clean, species_label)
meta$acc_simple <- gsub("GC._(.........)..", "\\1", meta$acc_clean)

# Read all all_16S.V4.cutadapt.fa
seqs <- readDNAStringSet("all_16S.V4.cutadapt.fa")

# Parse accession + copy from fasta headers
annot <- tibble(
  old_name = names(seqs),
  acc_clean = sub("^((GCA|GCF)_[0-9]+\\.[0-9]+).*", "\\1", names(seqs)),
  acc_simple = sub("GC._(.........)..", "\\1", acc_clean),
  copy = sub("^.*(_16S_copy[0-9]+).*", "\\1", names(seqs))
) %>%
  left_join(meta %>% select(acc_simple, species_label), by = "acc_simple")

# Create nicer fasta names
annot <- annot %>% mutate(species_tag = gsub("[^A-Za-z0-9]+", "_", species_label), new_name = paste(acc_clean, species_tag, copy, sep = "|"))
names(seqs) <- annot$new_name

# ----------------------------------
# 3) Align with DECIPHER
# ----------------------------------
# Alignment
aln_full <- AlignSeqs(seqs, processors = NULL, verbose = TRUE)

#  WE NEED TO REMOVE THE 2 SEQUENCES THAT WERE OBVIOUSLY BAD 
# GCA 000189835.2|Methanocatella smithii| 16S copy1 --> probably a contaminant, very different from the rest of the sequences
# GCA 041537865.1|Methanocatella smithii| 16S copy2 --> probably a contaminant, very different from the rest of the sequences

length(aln_full) # 1 sequences before removing the bad ones

aln_full <- aln_full[!names(aln_full) %in% c(
  "GCA_000189835.2|Methanocatella_smithii|_16S_copy1",
  "GCA_041537865.1|Methanocatella_smithii|_16S_copy2"
)]

length(aln_full) # 59 sequences after removing the bad ones

# ----------------------------------
# 4) Pairwise identity from DECIPHER
# ----------------------------------
# DistanceMatrix returns pairwise dissimilarities on aligned sequences.
# Convert to percent identity as (1 - distance) * 100.
dist_to_pid <- function(aln) {
  d <- DistanceMatrix(aln)
  pid <- (1 - d) * 100
  diag(pid) <- 100
  pid
}

pid_full <- dist_to_pid(aln_full)

write.csv(as.data.frame(pid_full), "pairwise_identity_V4_16S_DECIPHER.csv", row.names = TRUE)

# ----------------------------------
# 5) Summarise within / between species
# ----------------------------------
matrix_to_long <- function(m, region) {
  as.data.frame(as.table(m), stringsAsFactors = FALSE) %>%
    setNames(c("seq1", "seq2", "pid")) %>%
    filter(seq1 < seq2) %>%
    mutate(region = region)
}

long_full <- matrix_to_long(pid_full, "V4_16S")

head(as.data.frame(long_full), 2)
#                                                seq1
# 1 GCA_000016525.1|Methanocatella_smithii|_16S_copy1
# 2 GCA_000016525.1|Methanocatella_smithii|_16S_copy1
#                                                seq2 pid region
# 1 GCA_000016525.1|Methanocatella_smithii|_16S_copy2 100 V4_16S
# 2 GCA_000151225.1|Methanocatella_smithii|_16S_copy1 100 V4_16S

# extract species from the sequence names (second pipe-delimited field)
long_full$sp1 <- sub(".*\\|(.+)\\|_16S.*", "\\1", long_full$seq1)
long_full$sp2 <- sub(".*\\|(.+)\\|_16S.*", "\\1", long_full$seq2)

# extract accession (first field) to flag within-genome pairs
long_full$acc1 <- sub("\\|.*", "", long_full$seq1)
long_full$acc2 <- sub("\\|.*", "", long_full$seq2)

# sort species labels so sp_a <= sp_b (avoids duplicate A-B / B-A pairs)
long_full$sp_a <- pmin(long_full$sp1, long_full$sp2)
long_full$sp_b <- pmax(long_full$sp1, long_full$sp2)

# remove within-genome copy comparisons (same accession, different 16S copies)
df <- subset(long_full, acc1 != acc2)

# summarise
species_summary_full <- aggregate(pid ~ sp_a + sp_b, data = df, FUN = function(x) {
  c(n      = length(x),
    mean   = round(mean(x),   3),
    median = round(median(x), 3),
    min    = round(min(x),    3),
    max    = round(max(x),    3))
})

as.data.frame(do.call(cbind, list(
  species_summary_full[, c("sp_a","sp_b")],
  as.data.frame(species_summary_full$pid)
)))

#                       sp_a                       sp_b    n    mean  median
# 1   Methanocatella_smithii     Methanocatella_smithii 1317  99.894 100.000
# 2   Methanocatella_smithii   Methanocatella_smithii_A  312  99.947 100.000
# 3 Methanocatella_smithii_A   Methanocatella_smithii_A   14 100.000 100.000
# 4   Methanocatella_smithii Methanocatella_sp900766745   52  97.600  97.638
# 5 Methanocatella_smithii_A Methanocatella_sp900766745    6  97.638  97.638
#       min     max
# 1  97.244 100.000
# 2  97.638 100.000
# 3 100.000 100.000
# 4  96.063  97.638
# 5  97.638  97.638
```