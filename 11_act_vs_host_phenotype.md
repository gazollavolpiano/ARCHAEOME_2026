# Associations between ACT and host factors in FINRISK 2002

The file `df4g.RData` contains a dataframe with 3156 samples and 359 host variables. The variable `ArchaealCommunityType` contains the ACTs (M. smithii-dominanted, Ca. M. intestini-dominanted, M. sp900766745-dominanted, Diverse). We will use multinomial logistic regression to test for associations between ACT and host factors (diet, metabolites, liver enzymes, liver disease). 

```R
# load libraries
library(tidyverse)

# increase width of the console
options(width = 300)  

# load health data
load("df4g.RData")
dim(df4g) # 3156 x 359

# use multinomial logistic regression 
library(nnet)

# define the levels (reference = M. smithii-dominanted)
df4g$ArchaealCommunityType <- factor(df4g$ArchaealCommunityType, levels = c("M. smithii-dominanted", "Ca. M. intestini-dominanted", "M. sp900766745-dominanted", "Diverse"))

###############################
# DIET: FIBER_TOTAL and RED_MEAT models
################################

# fit the model
mm <- multinom(ArchaealCommunityType ~ FIBER_TOTAL + RED_MEAT + BL_AGE + MEN + BMI + CURR_SMOKE + BL_USE_RX_J01 + batch_name + ReadDepth_log10, data = df4g)

tidy_mm <- broom::tidy(mm, exponentiate = TRUE, conf.int = TRUE) %>% 
            filter(term != "(Intercept)") %>% 
            mutate(effect = paste0(round(estimate,2), " (", round(conf.low,2), ", ", round(conf.high,2), ")")) %>%
            select(y.level, term, effect, p.value) 

tidy_mm %>% write_csv("diet_multinom_results.csv")

################################
#          ACETATE
################################

# rank inverse normal transformation
df4g$acetate_rin <- qnorm((rank(df4g$ACETATE) - 0.5) / nrow(df4g))

# fit the model
mm <- multinom(ArchaealCommunityType ~ acetate_rin + BL_AGE + MEN + BMI + CURR_SMOKE + BL_USE_RX_J01 + batch_name + ReadDepth_log10, data = df4g)

tidy_mm <- broom::tidy(mm, exponentiate = TRUE, conf.int = TRUE) %>% 
            filter(term != "(Intercept)") %>% 
            mutate(effect = paste0(round(estimate,2), " (", round(conf.low,2), ", ", round(conf.high,2), ")")) %>%
            select(y.level, term, effect, p.value) 

tidy_mm %>% write_csv("acetate_multinom_results.csv")

# also plot the distribution of acetate across ACTs
svg("acetate_rin_ACT.svg", width = 2, height = 3)
ggplot(df4g, aes(x = ArchaealCommunityType, y = acetate_rin, fill = ArchaealCommunityType)) +
  geom_boxplot() +
  scale_fill_manual(values = c("M. smithii-dominanted" = "#7E6C9E", 
                               "Ca. M. intestini-dominanted" = "#181A34", 
                               "M. sp900766745-dominanted" = "#EBE572",  
                               "Diverse" = "#E5E5E0"), 
                    name="Archaeal community types") +
  labs(x = "", y = "Acetate (ranked)") +
  theme_classic(10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
dev.off()

################################
# GGT and ALT
################################

# check correlation
cor.test(df4g$GGT, df4g$ALAT, method = "spearman", use = "complete.obs", exact  = FALSE) # correlated (rho=0.6)

# transform the variables
df4g$log10GGT <- log10(df4g$GGT + 1)
df4g$log10ALAT <- log10(df4g$ALAT + 1)

# log10GGT
mm <- multinom(ArchaealCommunityType ~ log10GGT + BL_AGE + MEN + BMI + CURR_SMOKE + BL_USE_RX_J01 + batch_name + ReadDepth_log10, data = df4g)

tidy_mm <- broom::tidy(mm, exponentiate = TRUE, conf.int = TRUE) %>% 
            filter(term != "(Intercept)") %>% 
            mutate(effect = paste0(round(estimate,2), " (", round(conf.low,2), ", ", round(conf.high,2), ")")) %>%
            select(y.level, term, effect, p.value) 

tidy_mm %>% write_csv("log10GGT_multinom_results.csv")

svg("log10GGT_ACT.svg", width = 2, height = 3)
ggplot(df4g, aes(x = ArchaealCommunityType, y = log10GGT, fill = ArchaealCommunityType)) +
  geom_boxplot() +
  scale_fill_manual(values = c("M. smithii-dominanted" = "#7E6C9E", 
                               "Ca. M. intestini-dominanted" = "#181A34", 
                               "M. sp900766745-dominanted" = "#EBE572",  
                               "Diverse" = "#E5E5E0"), 
                    name="Archaeal community types") +
  labs(x = "", y = "log10(GGT+1)") +
  theme_classic(10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
dev.off()

# log10ALAT
mm <- multinom(ArchaealCommunityType ~ log10ALAT + BL_AGE + MEN + BMI + CURR_SMOKE + BL_USE_RX_J01 + batch_name + ReadDepth_log10, data = df4g)

tidy_mm <- broom::tidy(mm, exponentiate = TRUE, conf.int = TRUE) %>% 
            filter(term != "(Intercept)") %>% 
            mutate(effect = paste0(round(estimate,2), " (", round(conf.low,2), ", ", round(conf.high,2), ")")) %>%
            select(y.level, term, effect, p.value) 

tidy_mm %>% write_csv("log10ALAT_multinom_results.csv")

svg("log10ALAT_ACT.svg", width = 2, height = 3)
ggplot(df4g, aes(x = ArchaealCommunityType, y = log10ALAT, fill = ArchaealCommunityType)) +
  geom_boxplot() +
  scale_fill_manual(values = c("M. smithii-dominanted" = "#7E6C9E", 
                               "Ca. M. intestini-dominanted" = "#181A34", 
                               "M. sp900766745-dominanted" = "#EBE572",  
                               "Diverse" = "#E5E5E0"), 
                    name="Archaeal community types") +
  labs(x = "", y = "log10(ALAT+1)") +
  theme_classic(10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
dev.off()

##############################
#       Liver disease
################################
mm <- multinom(ArchaealCommunityType ~ PREVAL_LIVERDIS + ALKI2_FR02 + PREVAL_DIAB_T2 + BL_AGE + MEN + BMI + CURR_SMOKE + BL_USE_RX_J01 + batch_name + ReadDepth_log10, data = df4g)

tidy_mm <- broom::tidy(mm, exponentiate = TRUE, conf.int = TRUE) %>% 
            filter(term != "(Intercept)") %>% 
            mutate(effect = paste0(round(estimate,2), " (", round(conf.low,2), ", ", round(conf.high,2), ")")) %>%
            select(y.level, term, effect, p.value) 

tidy_mm %>% write_csv("PREVAL_LIVERDIS_multinom_results.csv")
```