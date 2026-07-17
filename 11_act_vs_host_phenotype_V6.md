# Associations between archaeal dominance groups and host factors in FINRISK 2002

- `health_data.RData` contains the host variables
- `phyloseq.rs` to define the groups (M. smithii-, M. intestini-, M. sp900766745-dominanted, diverse)

Covariate blocks:
`BASE`  : age + sex + BMI + smoking + recent antibacterial + batch + log10 read depth 
`ALC`   : + alcohol (part of the base model for GGT, ALAT, prevalent & incident LD) 
`DIET`  : + red meat + fibre 
`ENERGY`: + energy-dense food + Q57X (physical activity) 
`RX`    : + PPI + statin + education 
`DMET`  : + T2D x metformin as one 3-level factor (metformin nested in T2D), entering both as independent binaries is near-collinear 

Exclusion sets:
`-abx`   : recent antibacterial users 
`-t2dmet`: prevalent T2D or metformin
`-pld`   : prevalent liver disease 

Obs: When a perturbagen is excluded it leaves the adjustment set, so exclusion rows carry `RX` but not `DMET`, and drop antibacterial from `BASE`

Family 0-1 — archaeal dominance group as outcome (multinom):
0   : `BASE`                                ---> Base model 
1   : `BASE + DIET`                         ---> Primary diet model
1sa : `BASE + DIET - BMI`                   ---> Diet contingent on adiposity
1sb : `BASE + DIET + DMET + RX`             --->  Medication + SES confounding 
1sc : `BASE + DIET + ENERGY`                ---> Total-consumption / energy proxies 
1sd : `BASE + DIET + ENERGY + DMET + RX`    ---> Everything at once 

Family 2 — Acetate (SCFA):
2   : acetate ~ archaeal dominance group  + `BASE` (lm)       ---> Primary — archaeal dominance group as exposure 
2sa : archaeal dominance group ~ acetate + `BASE` (multinom)  ---> Alternative direction 
2sb : `BASE + DMET + RX`                                      --->  Medication + SES 
2sc : `[-abx -t2dmet -pld] BASE + RX`                         --->  Exclusion design 
2sd : `[-abx -t2dmet -pld] BASE + RX + red meat`              ---> Exclusion + red meat 

Family 3 — GGT:
3   : GGT ~ archaeal dominance group + `BASE + ALC` (lm)         ---> Primary 
3sa : archaeal dominance group ~ GGT + `BASE + ALC` (multinom)   ---> Alternative direction 
3sb : `BASE + ALC + DMET + RX`                                   ---> Medication + SES 
3sc : `[-abx -t2dmet -pld] BASE + ALC + RX`                      ---> Exclusion design 
3sd : `[-abx -t2dmet -pld] BASE + ALC + RX + red meat`           --->  Exclusion + red meat 
3se : `[-abx -t2dmet -pld] BASE - BMI + ALC + RX + red meat`     ---> BMI as mediator

Family 4 — ALAT:
4   : ALAT ~ archaeal dominance group + `BASE + ALC` (lm)         ---> Primary 
4sa : archaeal dominance group ~ ALAT + `BASE + ALC` (multinom)   ---> Alternative direction 
4sb : `BASE + ALC + DMET + RX`                                    ---> Medication + SES 
4sc : `[-abx -t2dmet -pld] BASE + ALC + RX`                       ---> Exclusion design 
4sd : `[-abx -t2dmet -pld] BASE + ALC + RX + red meat`            --->  Exclusion + red meat 
4se : `[-abx -t2dmet -pld] BASE - BMI + ALC + RX + red meat`      ---> BMI as mediator

Family 5 — Prevalent liver disease (glm):
5   : PrevLD ~ archaeal dominance group + `BASE + ALC + DMET`    ---> Primary (LD is outcome, cannot exclude it)
5sa : `BASE + ALC + DMET + RX`                                    ---> Medication + SES 
5sb : `archaeal dominance group(2-level) BASE + ALC + DMET`       ---> Powered contrast with M. smithii+ M. intestini (as one level) vs M. sp900766745, base
5sc : `archaeal dominance group(2-level) age + sex + BMI + ALC`   ---> Powered contrast with M. smithii+ M. intestini (as one level) vs M. sp900766745

Family 6 — Incident liver disease (Cox)
6   : `[-pld] BASE + ALC + DMET`                                   ---> Primary
6sa : `[-pld] BASE + ALC + DMET + RX`                               ---> Medication + SES
6sb : `[-pld] archaeal dominance group(2-level) + BASE + ALC + DMET`                  ---> Powered contrast, base
6sc : `[-pld] archaeal dominance group(2-level) + age + sex + BMI + ALC`              ---> Powered contrast, parsimonious

Exclusion design does not work for liver, we do not have enought cases!

energy-dense food because FINRISK 2002 has NO total energy intake (kcal/day) variable = 
- FR02_98 (number of meals and snacks on weekdays, ordinal 1 to 4)
- BREAD_SLICES_PER_DAY (count/day)
- JUNK_FOOD  (self-rated diet healthiness)

energy expenditure = 
- Q57X (Exercise at free-time, ordinal	1=not much;2=walk etc;3=regular 3hrs or more or competitive sports)

Reference group = M. smithii (first factor level) throughout

One CSV is written per FAMILY (1-6); each CSV holds the family's main model AND  all its sensitivity variants stacked long

```R
# ===========================================================================
# 0. SETUP
# ===========================================================================
library(tidyverse)
library(phyloseq)
library(nnet)
library(broom)
library(survival)

options(width = 220, stringsAsFactors = FALSE)

COHORT <- "FINRISK2002"
OUTDIR <- "/home/camigazo/11_act_vs_host_phenotype"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUTDIR)

load("~/data/health_data.RData")           # -> health$baseline, $covariates, $followup
physeq <- readRDS("~/data/phyloseq.rds")    

# ===========================================================================
# 1. BUILD THE ANALYSIS FRAME
# ===========================================================================

# --- 1.0 Dominance groups from the archaeal sub-community --------------------
FOCAL <- c(
  "Methanobrevibacter_A smithii"     = "M. smithii",
  "Methanobrevibacter_A smithii_A"   = "M. intestini",
  "Methanobrevibacter_A sp900766745" = "M. sp900766745"
)
GROUP_LEVELS <- c("M. smithii", "M. intestini", "M. sp900766745", "Diverse")

get_mat <- function(ps) {
  m <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) m <- t(m)
  m
}

arch_names <- taxa_names(physeq)[as(tax_table(physeq)[, "Domain"], "matrix")[, 1] == "Archaea"]
A  <- get_mat(physeq)[, arch_names, drop = FALSE]
A  <- A[rowSums(A) > 0, , drop = FALSE]
sp <- as(tax_table(physeq)[arch_names, "Species"], "matrix")[, 1]

top_i   <- max.col(A, ties.method = "first")
top_val <- A[cbind(seq_len(nrow(A)), top_i)]
n_tied  <- rowSums(A == top_val)
stopifnot(all(n_tied == 1))

dominant <- tibble(
  SampleID = rownames(A),
  top      = unname(sp[top_i]),
  top_frac = as.numeric(top_val / rowSums(A)),
  n_arch   = as.integer(rowSums(A > 0))
) %>%
  mutate(Group = factor(ifelse(top %in% names(FOCAL), FOCAL[top], "Diverse"),
                        levels = GROUP_LEVELS))

message(COHORT, " archaea-positive samples: ", nrow(dominant))   
# FINRISK2002 archaea-positive samples: 3156

print(table(dominant$Group))
# M. smithii   M. intestini M. sp900766745        Diverse
#       1670            363            541            582

# --- 1.1 Join host variables -------------------------------------------------
# ReadDepth column here, so read depth comes from ReadDepthMillions

data <- dominant %>%
  transmute(Barcode = SampleID, Group, top_frac, n_arch) %>%
  left_join(health$covariates %>%
              select(Barcode, ReadDepthMillions, BL_AGE, MEN, BMI,
                     BL_USE_RX_J01, batch_name), by = "Barcode") %>%
  left_join(health$baseline %>%
              select(Barcode, ALKI2_FR02, CURR_SMOKE, FIBER_TOTAL, RED_MEAT,
                     ACETATE, GGT, ALAT,
                     FR02_98, BREAD_SLICES_PER_DAY, JUNK_FOOD, Q57X,
                     KOULGR), by = "Barcode") %>%
  left_join(health$followup %>%
              select(Barcode, PREVAL_DIAB_T2, PREVAL_LIVERDIS,
                     INCIDENT_LIVERDIS, LIVERDIS_AGEDIFF,
                     BL_USE_RX_A02BC, BL_USE_RX_C10AA, BL_USE_RX_A10BA02),
            by = "Barcode")

dim(data)   # 3156 x 29

# --- 1.2 Transformations -----------------------------------------------------
# INCIDENT_LIVERDIS and LIVERDIS_AGEDIFF are kept numeric for Surv()

data <- data %>%
  mutate(
    Group = factor(Group, levels = GROUP_LEVELS),   # reference = M. smithii
    across(c(MEN, CURR_SMOKE, BL_USE_RX_J01,
             BL_USE_RX_A02BC, BL_USE_RX_C10AA, BL_USE_RX_A10BA02,
             PREVAL_DIAB_T2, PREVAL_LIVERDIS, KOULGR, batch_name), as.factor),
    ReadDepth_log10 = log10(ReadDepthMillions * 1e6),
    ALAT_log10      = log10(ALAT + 1),
    GGT_log10       = log10(GGT + 1),
    ACETATE_rin     = qnorm((rank(ACETATE, na.last = "keep") - 0.5) /
                              sum(!is.na(ACETATE))),
    INCIDENT_LIVERDIS = as.integer(as.character(INCIDENT_LIVERDIS)),
    LIVERDIS_AGEDIFF  = as.numeric(LIVERDIS_AGEDIFF)
  )

# --- 1.3 T2D x metformin as one 3-level factor (for DMET adjustment) ---------
data <- data %>%
  mutate(
    DIAB_MET = case_when(
      as.character(BL_USE_RX_A10BA02) == "1"                 ~ "Metformin",          # 54 (50 + 4)
      as.character(PREVAL_DIAB_T2)    == "1"                 ~ "T2D, no metformin",  # 51
      as.character(PREVAL_DIAB_T2)    == "0"                 ~ "No T2D",             # 2995
      TRUE                                                   ~ NA_character_          # 56
    ),
    DIAB_MET = factor(DIAB_MET, levels = c("No T2D", "T2D, no metformin", "Metformin")))

table(data$DIAB_MET, useNA = "ifany")
#No T2D T2D, no metformin         Metformin              <NA> 
#  2995                51                54                56 

with(data, table(PREVAL_DIAB_T2, BL_USE_RX_A10BA02, useNA = "ifany"))
#              BL_USE_RX_A10BA02
#PREVAL_DIAB_T2    0    1 <NA>
#          0    2995    4    0
#          1      51   50    0
#          <NA>    0    0   56

# --- 1.4 Covariate blocks (character vectors) --------------------------------
BASE       <- c("BL_AGE", "MEN", "BMI", "CURR_SMOKE", "BL_USE_RX_J01", "batch_name", "ReadDepth_log10")
BASE_noabx <- setdiff(BASE, "BL_USE_RX_J01") # antibacterial dropped (excluded)
ALC        <- "ALKI2_FR02"
DIET       <- c("FIBER_TOTAL", "RED_MEAT")
ENERGY     <- c("FR02_98", "BREAD_SLICES_PER_DAY", "JUNK_FOOD", "Q57X")
RX         <- c("BL_USE_RX_A02BC", "BL_USE_RX_C10AA", "KOULGR")   # PPI + statin + education
DMET       <- "DIAB_MET"

# --- 1.5 Exclusion helper ----------------------------------------------------
is1 <- function(x) as.character(x) %in% "1"

apply_excl <- function(dd, abx = FALSE, t2dmet = FALSE, pld = FALSE) {
  if (abx)    dd <- dd[!is1(dd$BL_USE_RX_J01), ]
  if (t2dmet) dd <- dd[!is1(dd$PREVAL_DIAB_T2) & !is1(dd$BL_USE_RX_A10BA02) & !is.na(dd$PREVAL_DIAB_T2), ]
  if (pld)    dd <- dd[!is1(dd$PREVAL_LIVERDIS), ]
  dd
}

# ===========================================================================
# 2. MODELLING ENGINE
# ===========================================================================
# A "spec" is one model: label, engine, data subset, outcome, RHS terms, and a
# `focal` regex selecting the terms to echo to the console. run_spec() fits it,
# tidies ALL coefficients for the CSV, prints the focal rows, and returns the
# tidy frame. run_family() stacks a list of specs and writes one CSV.

mk <- function(label, engine, dd, outcome, rhs, focal, effect) {
  list(label = label, engine = engine, dd = dd, focal = focal, effect = effect,
       formula = as.formula(paste(outcome, "~", paste(rhs, collapse = " + "))))
}

run_spec <- function(s) {
  dd   <- s$dd
  form <- s$formula
  environment(form) <- environment()   # formula carries this frame (where dd lives)
  
  fit <- tryCatch(
    switch(s$engine,
           multinom = do.call(nnet::multinom,   list(form, data = dd, trace = FALSE)),
           lm       = do.call(stats::lm,        list(form, data = dd)),
           glm      = do.call(stats::glm,       list(form, data = dd, family = binomial())),
           coxph    = do.call(survival::coxph,  list(form, data = dd))),
    error = function(e) { message("  [", s$label, "] FIT FAILED: ", conditionMessage(e)); NULL })
  if (is.null(fit)) return(NULL)
  
  expo <- s$engine %in% c("multinom", "glm", "coxph")
  out  <- broom::tidy(fit, exponentiate = expo, conf.int = TRUE)
  if (!"y.level" %in% names(out)) out <- mutate(out, y.level = NA_character_)
  
  n_used <- if (s$engine == "coxph") fit$n     else nobs(fit)
  n_ev   <- if (s$engine == "coxph") fit$nevent else NA_integer_
  
  out <- out %>%
    mutate(model = s$label, engine = s$engine, effect = s$effect,
           n = n_used, events = n_ev, .before = 1)
  
  foc <- out %>% filter(str_detect(term, s$focal))
  cat("\n---- ", s$label, "  [", s$effect, "]  n=", n_used,
      if (s$engine == "coxph") paste0("  events=", n_ev) else "", " ----\n", sep = "")
  if (nrow(foc) == 0) {
    cat("   (focal term not estimated)\n")
  } else {
    foc %>%
      mutate(value = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high),
             p     = formatC(p.value, format = "g", digits = 3)) %>%
      select(any_of("y.level"), term, value, p) %>%
      as.data.frame() %>% print(row.names = FALSE)
  }
  out
}

run_family <- function(specs, csv) {
  cat("\n\n=======================================================\n")
  cat("  ", csv, "\n")
  cat("=======================================================\n")
  res <- purrr::map_dfr(specs, run_spec)
  readr::write_csv(res, csv)
  cat("\n[written] ", csv, "  (", nrow(res), " rows)\n", sep = "")
  invisible(res)
}

# ===========================================================================
# 3. FAMILY 1 : archaeal dominance group as OUTCOME (multinom)
# ===========================================================================
# Focal = the diet / adiposity terms. Effect scale = OR (relative-risk ratios)

fam1 <- list(
  mk("0",   "multinom", data, "Group", BASE,                                  "RED_MEAT|FIBER_TOTAL|BMI", "OR"),
  mk("1",   "multinom", data, "Group", c(BASE, DIET),                         "RED_MEAT|FIBER_TOTAL|BMI", "OR"),
  mk("1sa", "multinom", data, "Group", setdiff(c(BASE, DIET), "BMI"),         "RED_MEAT|FIBER_TOTAL",     "OR"),
  mk("1sb", "multinom", data, "Group", c(BASE, DIET, DMET, RX),               "RED_MEAT|FIBER_TOTAL|BMI", "OR"),
  mk("1sc", "multinom", data, "Group", c(BASE, DIET, ENERGY),                 "RED_MEAT|FIBER_TOTAL|BMI", "OR"),
  mk("1sd", "multinom", data, "Group", c(BASE, DIET, ENERGY, DMET, RX),       "RED_MEAT|FIBER_TOTAL|BMI", "OR")
)
res_fam1 <- run_family(fam1, "family1_group_diet.csv")

# ===========================================================================
# 4. FAMILY 2 : Acetate (SCFA)
# ===========================================================================
# Primary direction: acetate ~ Group (lm), Group as exposure. 2sa flips direction.

d2_excl <- apply_excl(data, abx = TRUE, t2dmet = TRUE, pld = TRUE)

fam2 <- list(
  mk("2",   "lm",       data,    "ACETATE_rin", c("Group", BASE),                          "Group",       "beta"),
  mk("2sa", "multinom", data,    "Group",       c("ACETATE_rin", BASE),                    "ACETATE_rin", "OR"),
  mk("2sb", "lm",       data,    "ACETATE_rin", c("Group", BASE, DMET, RX),                "Group",       "beta"),
  mk("2sc", "lm",       d2_excl, "ACETATE_rin", c("Group", BASE_noabx, RX),                "Group",       "beta"),
  mk("2sd", "lm",       d2_excl, "ACETATE_rin", c("Group", BASE_noabx, RX, "RED_MEAT"),    "Group",       "beta")
)
res_fam2 <- run_family(fam2, "family2_acetate.csv")

# ===========================================================================
# 5. FAMILY 3 : GGT
# ===========================================================================

d3_excl <- apply_excl(data, abx = TRUE, t2dmet = TRUE, pld = TRUE)

fam3 <- list(
  mk("3",   "lm",       data,    "GGT_log10", c("Group", BASE, ALC),                                 "Group",     "beta"),
  mk("3sa", "multinom", data,    "Group",     c("GGT_log10", BASE, ALC),                             "GGT_log10", "OR"),
  mk("3sb", "lm",       data,    "GGT_log10", c("Group", BASE, ALC, DMET, RX),                       "Group",     "beta"),
  mk("3sc", "lm",       d3_excl, "GGT_log10", c("Group", BASE_noabx, ALC, RX),                       "Group",     "beta"),
  mk("3sd", "lm",       d3_excl, "GGT_log10", c("Group", BASE_noabx, ALC, RX, "RED_MEAT"),           "Group",     "beta"),
  mk("3se", "lm",       d3_excl, "GGT_log10", c("Group", setdiff(BASE_noabx, "BMI"), ALC, RX, "RED_MEAT"), "Group", "beta")
)
res_fam3 <- run_family(fam3, "family3_ggt.csv")

# ===========================================================================
# 6. FAMILY 4 : ALAT
# ===========================================================================

d4_excl <- apply_excl(data, abx = TRUE, t2dmet = TRUE, pld = TRUE)

fam4 <- list(
  mk("4",   "lm",       data,    "ALAT_log10", c("Group", BASE, ALC),                                 "Group",      "beta"),
  mk("4sa", "multinom", data,    "Group",      c("ALAT_log10", BASE, ALC),                            "ALAT_log10", "OR"),
  mk("4sb", "lm",       data,    "ALAT_log10", c("Group", BASE, ALC, DMET, RX),                       "Group",      "beta"),
  mk("4sc", "lm",       d4_excl, "ALAT_log10", c("Group", BASE_noabx, ALC, RX),                       "Group",      "beta"),
  mk("4sd", "lm",       d4_excl, "ALAT_log10", c("Group", BASE_noabx, ALC, RX, "RED_MEAT"),           "Group",      "beta"),
  mk("4se", "lm",       d4_excl, "ALAT_log10", c("Group", setdiff(BASE_noabx, "BMI"), ALC, RX, "RED_MEAT"), "Group", "beta")
)
res_fam4 <- run_family(fam4, "family4_alat.csv")


# ===========================================================================
# 7. FAMILY 5 : Prevalent liver disease (glm, logistic)
# ===========================================================================
# LD is the outcome, so it cannot be excluded. Exclusion arms were dropped: with
# only ~21 prevalent cases, -t2dmet leaves too few to model. The four-level model
# (5) separates in the M. intestini cell (0 or near-0 cases) -> run it for
# completeness but report the 2-level contrast (5sb/5sc) as the headline, and do
# NOT put the four-level M. intestini row in a table.

# ---------------------------------------------------------------------------
# Shared: collapse to 2-level exposure (M. smithii + M. intestini vs sp900766745)
# ---------------------------------------------------------------------------
make_2lvl <- function(dd) {
  dd %>%
    filter(Group != "Diverse") %>%
    mutate(Group2 = fct_collapse(droplevels(Group),
                                 "M. smithii/intestini" = c("M. smithii", "M. intestini")),
           Group2 = relevel(Group2, ref = "M. smithii/intestini"))
}

data_2lvl <- make_2lvl(data)

fam5 <- list(
  mk("5",   "glm", data,      "PREVAL_LIVERDIS", c("Group",  BASE, ALC, DMET),             "Group",  "OR"),
  mk("5sa", "glm", data,      "PREVAL_LIVERDIS", c("Group",  BASE, ALC, DMET, RX),         "Group",  "OR"),
  mk("5sb", "glm", data_2lvl, "PREVAL_LIVERDIS", c("Group2", BASE, ALC, DMET),             "Group2", "OR"),
  mk("5sc", "glm", data_2lvl, "PREVAL_LIVERDIS", c("Group2", "BL_AGE", "MEN", "BMI", ALC), "Group2", "OR")
)
res_fam5 <- run_family(fam5, "family5_prevalent_LD.csv")

cat("  full (4-level): ", sum(is1(data$PREVAL_LIVERDIS)), " cases / ", nrow(data), "\n", sep = "")
#full (4-level): 21 cases / 3156

print(with(data,      table(Group,  PREVAL_LIVERDIS)))
#PREVAL_LIVERDIS
#Group               0    1
#M. smithii     1630    6
#M. intestini    356    1
#M. sp900766745  523    8
#Diverse         570    6

print(with(data_2lvl, table(Group2, PREVAL_LIVERDIS)))
#PREVAL_LIVERDIS
#Group2                    0    1
#M. smithii/intestini 1986    7
#M. sp900766745        523    8

# ===========================================================================
# 8. FAMILY 6 : Incident liver disease (Cox)
# ===========================================================================
# Risk set = free of prevalent LD at baseline, valid follow-up time.

d_inc <- data %>%
  filter(!is1(PREVAL_LIVERDIS), !is.na(LIVERDIS_AGEDIFF), LIVERDIS_AGEDIFF > 0,
         !is.na(Group))

# ---- POWER CHECK: read this before trusting any HR ----
print(table(d_inc$Group, d_inc$INCIDENT_LIVERDIS))
#0    1
#M. smithii     1603   27
#M. intestini    347    9
#M. sp900766745  512   11
#Diverse         566    4

cat("total incident events: ", sum(d_inc$INCIDENT_LIVERDIS == 1, na.rm = TRUE), "\n", sep = "")
# total incident events: 51

# 6sb/6sc: two-level contrast on the incident risk set (always -pld by design)
d_inc_2lvl <- make_2lvl(d_inc)

surv_lhs <- "Surv(LIVERDIS_AGEDIFF, INCIDENT_LIVERDIS)"

fam6 <- list(
  mk("6",   "coxph", d_inc,      surv_lhs, c("Group",  BASE, ALC, DMET),              "Group",  "HR"),
  mk("6sa", "coxph", d_inc,      surv_lhs, c("Group",  BASE, ALC, DMET, RX),          "Group",  "HR"),
  mk("6sb", "coxph", d_inc_2lvl, surv_lhs, c("Group2", BASE, ALC, DMET),              "Group2", "HR"),
  mk("6sc", "coxph", d_inc_2lvl, surv_lhs, c("Group2", "BL_AGE", "MEN", "BMI", ALC),  "Group2", "HR")
)
res_fam6 <- run_family(fam6, "family6_incident_LD.csv")

# --- 8.1 Proportional-hazards check (reported model 6sc, two-level) ----------
f_ph <- as.formula(paste(surv_lhs, "~",
                         paste(c("Group2", "BL_AGE", "MEN", "BMI", ALC), collapse = " + ")))
d_ph <- d_inc_2lvl[complete.cases(d_inc_2lvl[, all.vars(f_ph)]), ]
fit_ph <- coxph(f_ph, data = d_ph)
zp <- tryCatch(cox.zph(fit_ph), error = function(e) {
  message("cox.zph failed: ", conditionMessage(e)); NULL })

if (!is.null(zp)) {
  cat("\nPH test (model 6sc):\n"); print(zp$table)
  svg("cox_zph_group2.svg", width = 8, height = 6)
  par(mfrow = c(2, 2)); plot(zp)
  dev.off()
}

#PH test (model 6sc):
#  chisq df         p
#Group2     0.521154688  1 0.4703495
#BL_AGE     2.217230235  1 0.1364777
#MEN        0.578502214  1 0.4469000
#BMI        0.004915217  1 0.9441072
#ALKI2_FR02 0.227388338  1 0.6334671
#GLOBAL     3.790567897  5 0.5799457


# ===========================================================================
# 9. EXTRA: lm diagnostics + outcome histograms (primary enzyme models)
# ===========================================================================
diag_lm <- function(outcome, rhs, tag) {
  f <- as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
  m <- lm(f, data = data)
  svg(paste0("diag_", tag, ".svg"), width = 8, height = 8)
  par(mfrow = c(2, 2)); plot(m)
  dev.off()
}
diag_lm("ACETATE_rin", c("Group", BASE),        "acetate")
diag_lm("GGT_log10",   c("Group", BASE, ALC),   "ggt")
diag_lm("ALAT_log10",  c("Group", BASE, ALC),   "alat")

for (v in c("ACETATE_rin", "GGT_log10", "ALAT_log10")) {
  svg(paste0("hist_", v, ".svg"), width = 2, height = 4)
  print(
    ggplot(data, aes(x = .data[[v]])) +
      geom_histogram(aes(y = after_stat(density)), bins = 40,
                     fill = "grey80", colour = "white") +
      stat_function(fun = dnorm,
                    args = list(mean = mean(data[[v]], na.rm = TRUE),
                                sd   = sd(data[[v]],   na.rm = TRUE)),
                    colour = "red", linewidth = 0.7) +
      labs(x = v, y = "Density") + theme_bw(11)
  )
  dev.off()
}

# CONFIRM NUMBERS: 
bind_rows(res_fam1, res_fam2, res_fam3, res_fam4, res_fam5, res_fam6) %>%
  distinct(model, engine, n, events) %>%
  arrange(model) %>%
  print(n = Inf)
```