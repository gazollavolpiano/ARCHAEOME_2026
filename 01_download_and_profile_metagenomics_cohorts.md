# Downloading samples from public metagenomics cohorts

A total of 14 publicly available metagenomic cohorts with stool samples were selected for inclusion in this study.

Below is an example for the KarlssonFH_2013 cohort:

```R
# libraries
library(tidyverse)
library(curatedMetagenomicData)
packageVersion("curatedMetagenomicData") # 3.6.2

# list samples to download 
metadata <- sampleMetadata %>% 
filter(body_site == "stool") %>% 
filter(study_name == "KarlssonFH_2013") 

dim(metadata) # 145 samples

# save the list of NCBI_accession
accession_list <- unlist(strsplit(metadata$NCBI_accession, split = ";")) # split the strings into accession numbers 
writeLines(accession_list, "ncbi_accessions.txt") # write to a line in a text file
```

Use ncbi_accessions.txt with the SRA Toolkit to download FASTQ files.

```bash
# download the data
fastq-dump --version # 3.0.5
cat ncbi_accessions.txt | parallel -j 10 fastq-dump --split-files --gzip {}

# check the number of files downloaded
ls *1.fastq.gz | wc -l # 145
ls *2.fastq.gz | wc -l # 145
```

# Taxonomic classification of metagenomic reads

To obtain species-level taxonomic profiles from the metagenomes, we used a two-step Kraken2/Bracken workflow. First, raw paired-end reads were screened against a human pangenome (HPRC) database with Kraken2, and non-human reads were extracted using KrakenTools. These filtered reads were then classified against a GTDB r214 database with Kraken2, and species-level abundances were estimated with Bracken. Finally, Bracken reports from all samples were combined into a single phyloseq file for downstream analyses.

```bash
# versions
kraken2 --version # 2.1.3
bracken -v # 2.9

# kraken2 with HPRC indices
kraken2 --db /k2_HPRC_20230810/db/ --paired --gzip-compressed SAMPLE_R1.fastq.gz SAMPLE_R2.fastq.gz --output SAMPLE.HPRC.kraken

# read extraction with KrakenTools (extract_kraken_reads.py)
python extract_kraken_reads.py -k SAMPLE.HPRC.kraken -s1 SAMPLE_R1.fastq.gz -s2 SAMPLE_R2.fastq.gz -t 0 -o nonhuman_SAMPLE.1.fastq -o2 nonhuman_SAMPLE.2.fastq --fastq-output

# kraken2 with GTDB indices (128gb database, confidence 0.1)
kraken2 --db /gtdb_r214/128gb/ --output SAMPLE.kraken --report SAMPLE.kraken.report --paired nonhuman_SAMPLE.1.fastq nonhuman_SAMPLE.2.fastq --confidence 0.1 --use-names

# run bracken 
bracken -d /gtdb_r214/128gb/ -i SAMPLE.kraken.report -o SAMPLE.bracken.report -r 150 -l S -t 25 # read length adjusted based on sequencing strategy!

# combine outputs in a single file
python combine_bracken_outputs_mod.py --files *bracken.report -o combined_bracken.report
```

Use R to create a phyloseq object from the combined bracken output file.

```R
# library
library(phyloseq)
library(tidyverse)

# read merged bracken output result and modify it
merged_bracken_output <- read.table("combined_bracken.report", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
                          select(-matches("frac|taxonomy_id|taxonomy_lvl")) %>% 
                          rename(Species = name) 

merged_bracken_output <- merged_bracken_output %>% rename_all(~ str_remove(., ".bracken.report_num"))

# read the GTDB taxonomy file (available at https://data.gtdb.ecogenomic.org/releases/release214/214.0/)
taxonomy <- read.table("gtdb_v214_taxonomy.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# add the full taxonomy to the merged bracken output 
merged_bracken_output <- merge(taxonomy, merged_bracken_output, all.y = TRUE)

# create otu table and tax table
otu <- merged_bracken_output %>% 
      select(-Domain, -Phylum, -Class, -Order, -Family, -Genus, Species) %>%
      column_to_rownames("Species")

tax <- merged_bracken_output %>% select(Domain, Phylum, Class, Order, Family, Genus, Species) 
rownames(tax) <- tax$Species
 
# create phyloseq objectq
physeq <- phyloseq(otu_table(as.matrix(otu), taxa_are_rows = TRUE), tax_table(as.matrix(tax)))

# write phyloseq object
saveRDS(physeq, "phyloseq.rds")
```