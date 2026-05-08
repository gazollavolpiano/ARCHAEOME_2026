## Study overview

This repository contains the analysis code supporting the manuscript characterising the **human gut archaeome** across **13,716 stool metagenomes** from **15 cohorts** in **12 countries** and three continents.

- **15 cohorts** harmonised through a single pipeline against a custom **GTDB R214** database
- Main cohort: **FINRISK 2002** (*n* = 7,226), a population-based Finnish cohort (used for association with host phenotypes)
- Archaeal communities partitioned into four mutually exclusive **archaeal community types (ACTs)** based on the dominant species
- Cross-domain ecology assessed via **PERMANOVA** on bacterial community composition and **FastSpar** correlations between dominant *Methanobrevibacter* species and bacterial genera
- Host-trait associations in FINRISK 2002 modelled with **multinomial logistic regression** (dietary intake, BMI, smoking, antibiotic use, circulating acetate, liver enzymes, prevalent liver disease)
- Functional basis of ACTs explored via **comparative genomics** (Gapseq) and **community-scale metabolic modelling** (MICOM)

## Bioinformatics pipeline summary

```
Raw reads (Illumina paired-end, stool metagenomes; 15 cohorts)
        │
        ▼
  Kraken2 (HPRC database) → KrakenTools (human read removal)
        │
        ▼
  Kraken2 (custom GTDB R214 database, confidence = 0.1)
        │
        ▼
  Bracken (species-level, read length 150, threshold 25)
        │
        ▼
  R: combine_bracken_outputs → phyloseq object
        │
        ├──► Archaeal detection per sample
        │         │
        │         ├──► Cohort-level prevalence summaries 
        │         │
        │         └──► Mixed-effects logistic regression of archaeal detection
        │
        ├──► ACT assignment per sample (dominant archaeal species)
        │         │
        │         ├──► ComplexHeatmap (CLR-transformed ACT heatmap)
        │         └──► Geographic distribution of ACTs
        │
        ├──► GTDB R214 archaeal reference tree → Microreact
        │         (phylogeny of prevalent gut archaea)
        │
        ├──► Bacterial community composition (genus-level, CLR-transformed)
        │         │
        │         ├──► PCA 
        │         └──► PERMANOVA pairwise across ACTs, adjusted for read depth
        │
        ├──► FastSpar (SparCC) correlations: Methanobrevibacter species vs bacterial genera
        │         │
        │         ├──► 200 inference iterations, 10,000 bootstraps
        │         └──► FDR-corrected per cohort (Figure 4)
        │
        └──► FINRISK 2002 host-trait modelling (n = 3,156)
                  │
                  ├──► Multinomial logistic regression 
                  │    Outcome: ACT (M. smithii–dominated as reference)
                  │    Predictors: diet, BMI, smoking, antibiotic use,
                  │                acetate, liver enzymes, prevalent liver disease
                  │    
                  └──► Adjusted for age, sex, sequencing batch, log10 read depth

Methanobrevibacter representative genomes (M. smithii, M. intestini, M. sp900766745)
        │
        ├──► Gapseq (reaction and transporter prediction)
        │         │
        │         └──► Comparative reaction repertoire 
        │
        └──► MICOM community-scale flux balance analysis (PREDICT 1, matched pairs)
                  │
                  ├──► MatchIt (1:1 matching on age, sex, BMI, sequencing depth)
                  ├──► Community models from UHGG + human gut archaeal genome catalogue
                  └──► LASSO classifier on metabolic flux profiles 
```

## Software and versions

| Tool | Version | Purpose |
|---|---|---|
| NCBI SRA Toolkit | 3.0.5 | Download of raw metagenomic reads from SRA |
| Kraken2 | 2.1.3 | Metagenomic taxonomic classification (HPRC + GTDB R214) |
| KrakenTools | — | Removal of human-classified reads |
| Bracken | 2.9 | Species-level abundance re-estimation (read length 150, min 25 reads) |
| FastSpar | 1.0.0 | SparCC correlation inference with bootstrap p-values |
| Gapseq | 1.4.0 | Genome-scale reaction and transporter prediction |
| HMMER | 3.1b1 | Marker gene identification (archaeal phylogeny) |
| Prodigal | 2.6.3 | Gene calling (archaeal phylogeny) |
| IQ-TREE | 2.1.2 | Maximum-likelihood phylogenetic inference (PMSF model) |
| FastTree | 2.1.10 | Initial guide tree for IQ-TREE |
| MICOM (via QIIME2 wrapper) | 0.16.0 | Community-scale flux balance analysis |
| curatedMetagenomicData | 3.6.2 | Public cohort metadata and SRA accession retrieval |
