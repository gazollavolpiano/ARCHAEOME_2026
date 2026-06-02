# Associations between ACT and host factors in FINRISK 2002

The file `health_data.RData` contains the host variables. The file `FINRISK_2002_archcommunity_types_read_depth.csv` contains the variable `ArchaealCommunityType` with the ACTs (M. smithii-dominanted, Ca. M. intestini-dominanted, M. sp900766745-dominanted, Diverse). 

```
# load libs
library(tidyverse)
library(nnet)
library(broom)

# read/load files of interest
load("health_data.RData")

act <- read.csv("FINRISK_2002_archcommunity_types_read_depth.csv")
head(act, 2)
# SampleID Archaea_detected Species ArchaealCommunityType ReadDepth        Study
#1 .........-.                0    <NA>                  <NA>   9168344 FINRISK 2002
#2 .........-.                0    <NA>                  <NA>  10609149 FINRISK 2002

# prepare dataframe for associations
data <- act %>%
          rename(Barcode = SampleID) %>%
          filter(Archaea_detected == 1) %>%
          select(Barcode, ArchaealCommunityType, ReadDepth) %>%
          left_join(health$covariates %>% select(Barcode, ReadDepthMillions, BL_AGE, MEN, BMI, BL_USE_RX_J01, batch_name)) %>%
          left_join(health$baseline %>% select(Barcode, ALKI2_FR02, CURR_SMOKE, FIBER_TOTAL, RED_MEAT, ACETATE, GGT, ALAT)) %>%
          left_join(health$followup %>% select(Barcode, PREVAL_DIAB_T2, PREVAL_LIVERDIS))

head(data,2)
dim(data) # 3156   18

# prepare data transformations
data$ArchaealCommunityType <- factor(data$ArchaealCommunityType, levels=c("M. smithii-dominanted", "Ca. M. intestini-dominanted", "M. sp900766745-dominanted", "Diverse"))
data$PREVAL_LIVERDIS <- as.factor(data$PREVAL_LIVERDIS)
data$PREVAL_DIAB_T2 <- as.factor(data$PREVAL_DIAB_T2)
data$MEN <- as.factor(data$MEN)
data$ReadDepth_log10 <- log10(data$ReadDepth)
data$ACETATE_rin <- qnorm((rank(data$ACETATE, na.last = "keep") - 0.5) / sum(!is.na(data$ACETATE)))
data$ALAT_log10 <- log10(data$ALAT+1)
data$GGT_log10 <- log10(data$GGT+1)

#################
# Model 0: BASE MODEL
##################

mm <- multinom(ArchaealCommunityType ~ BL_AGE + MEN + BMI + CURR_SMOKE + BL_USE_RX_J01 + batch_name + ReadDepth_log10, data) 
broom::tidy(mm, exponentiate = TRUE, conf.int = TRUE) %>% write_csv("basic_result.csv")
nobs(mm) #3067

#################
# RED MEAT and FIBER: CHECK CORRELATION
##################

cor.test(data$RED_MEAT, data$FIBER_TOTAL, method="spearman", exact=FALSE)
#Spearman's rank correlation rho
#
#data:  data$RED_MEAT and data$FIBER_TOTAL
#S = 4282279512, p-value = 0.006789
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#       rho 
#-0.0502319 
sum(complete.cases(data$RED_MEAT, data$FIBER_TOTAL))
# 2903

#################
# MODEL 1: RED MEAT and FIBER (PRIMARY DIET)
##################

mm <- multinom(ArchaealCommunityType ~ FIBER_TOTAL + RED_MEAT + BL_AGE + MEN + BMI + CURR_SMOKE + BL_USE_RX_J01 + batch_name + ReadDepth_log10, data) 
broom::tidy(mm, exponentiate = TRUE, conf.int = TRUE) %>% write_csv("diet_fiber_meat_result.csv")
nobs(mm) #2841

#################
# MODEL 1s: RED MEAT and FIBER, NO BMI
##################

mm <- multinom(ArchaealCommunityType ~ FIBER_TOTAL + RED_MEAT + BL_AGE + MEN + CURR_SMOKE + BL_USE_RX_J01 + batch_name + ReadDepth_log10, data) 
broom::tidy(mm, exponentiate = TRUE, conf.int = TRUE) %>% write_csv("diet_fiber_meat_nobmi_result.csv")
nobs(mm) #2841

#################
# MODEL 2: ACETATE (ArchaealCommunityType as the outcome)
##################

mm <- multinom(ArchaealCommunityType ~ ACETATE_rin + BL_AGE + MEN + BMI + CURR_SMOKE + BL_USE_RX_J01 + batch_name + ReadDepth_log10, data) 
broom::tidy(mm, exponentiate = TRUE, conf.int = TRUE) %>% write_csv("acetate_result_multinom.csv")
nobs(mm) #2805

#################
# MODEL 2s: ACETATE as outcome (lm)
##################

mlm <- lm(ACETATE_rin ~ ArchaealCommunityType + BL_AGE + MEN + BMI + CURR_SMOKE + BL_USE_RX_J01 + batch_name + ReadDepth_log10, data) 
broom::tidy(mlm, conf.int = TRUE) %>% write_csv("acetate_result_acetate_outcome_lm.csv")
nobs(mlm) #2805

svg("acetate_result_acetate_outcome_lm_diagnostics.svg", width = 8, height = 8)
par(mfrow=c(2,2)); plot(mlm) 
dev.off()

svg("hist_acetate_rin.svg", width = 2, height = 4)
ggplot(data, aes(x=ACETATE_rin))+
  geom_histogram(aes(y=after_stat(density)), bins=40,
                fill="grey80", colour="white")+
  stat_function(fun=dnorm,
                args = list(mean=mean(data$ACETATE_rin, na.rm=TRUE),
                            sd=sd(data$ACETATE_rin, na.rm=TRUE)),
                colour="red", linewidth=0.7)+
  labs(x="Acetate (ranked", y="Density")+
  theme_bw(11)
dev.off()

#################
# MODEL 3: GGT (ArchaealCommunityType as the outcome + ALKI2_FR02)
##################

mm <- multinom(ArchaealCommunityType ~ GGT_log10 + BL_AGE + MEN + BMI + CURR_SMOKE + BL_USE_RX_J01 + ALKI2_FR02 + batch_name + ReadDepth_log10, data) 
broom::tidy(mm, exponentiate = TRUE, conf.int = TRUE) %>% write_csv("GGT_result_multinom.csv")
nobs(mm) #2946

#################
# MODEL 3s: GGT as outcome (lm) + ALKI2_FR02
##################

mlm <- lm(GGT_log10 ~ ArchaealCommunityType + BL_AGE + MEN + BMI + CURR_SMOKE + BL_USE_RX_J01 + ALKI2_FR02 + batch_name + ReadDepth_log10, data) 
broom::tidy(mlm, conf.int = TRUE) %>% write_csv("GGT_result_GGT_outcome_lm.csv")
nobs(mlm) #2946

svg("GGT_result_GGT_outcome_lm_diagnostics.svg", width = 8, height = 8)
par(mfrow=c(2,2)); plot(mlm) 
dev.off()

svg("hist_GGT_log10.svg", width = 2, height = 4)
ggplot(data, aes(x=GGT_log10))+
  geom_histogram(aes(y=after_stat(density)), bins=40,
                 fill="grey80", colour="white")+
  stat_function(fun=dnorm,
                args = list(mean=mean(data$GGT_log10, na.rm=TRUE),
                            sd=sd(data$GGT_log10, na.rm=TRUE)),
                colour="red", linewidth=0.7)+
  labs(x="log10(GGT+1)", y="Density")+
  theme_bw(11)
dev.off()

#################
# MODEL 4: ALAT (ArchaealCommunityType as the outcome + ALKI2_FR02)
##################

mm <- multinom(ArchaealCommunityType ~ ALAT_log10 + BL_AGE + MEN + BMI + CURR_SMOKE + BL_USE_RX_J01 + ALKI2_FR02 + batch_name + ReadDepth_log10, data) 
broom::tidy(mm, exponentiate = TRUE, conf.int = TRUE) %>% write_csv("ALAT_result_multinom.csv")
nobs(mm) #2938

#################
# MODEL 4s: ALAT as outcome (lm) + ALKI2_FR02
##################

mlm <- lm(ALAT_log10 ~ ArchaealCommunityType + BL_AGE + MEN + BMI + CURR_SMOKE + BL_USE_RX_J01 + ALKI2_FR02 + batch_name + ReadDepth_log10, data) 
broom::tidy(mlm, conf.int = TRUE) %>% write_csv("ALAT_result_ALAT_outcome_lm.csv")
nobs(mlm) #2938

svg("ALAT_result_ALAT_outcome_lm_diagnostics.svg", width = 8, height = 8)
par(mfrow=c(2,2)); plot(mlm) 
dev.off()

svg("hist_ALAT_log10.svg", width = 2, height = 4)
ggplot(data, aes(x=ALAT_log10))+
  geom_histogram(aes(y=after_stat(density)), bins=40,
                 fill="grey80", colour="white")+
  stat_function(fun=dnorm,
                args = list(mean=mean(data$ALAT_log10, na.rm=TRUE),
                            sd=sd(data$ALAT_log10, na.rm=TRUE)),
                colour="red", linewidth=0.7)+
  labs(x="log10(ALAT+1)", y="Density")+
  theme_bw(11)
dev.off()

#################
# MODEL 5: LD as the outcome, GLM model, + ALKI2_FR02 + PREVAL_DIAB_T2
##################

mglm <- glm(PREVAL_LIVERDIS ~ ArchaealCommunityType + BL_AGE + MEN + BMI + CURR_SMOKE + BL_USE_RX_J01 + 
              ALKI2_FR02 + PREVAL_DIAB_T2 + batch_name + ReadDepth_log10, data,
               family=binomial()) 
write.csv(file="LD_result_glm.csv", broom::tidy(mglm, exponentiate = TRUE, conf.int = TRUE))
nobs(mglm) #2946

#################
# MODEL 5s: LD as the outcome, GLM model, + ALKI2_FR02 + PREVAL_DIAB_T2, NO BATCH INCLUDED
##################

mglm <- glm(PREVAL_LIVERDIS ~ ArchaealCommunityType + BL_AGE + MEN + BMI + CURR_SMOKE + BL_USE_RX_J01 + 
              ALKI2_FR02 + PREVAL_DIAB_T2  + ReadDepth_log10, data,
            family=binomial()) 
write.csv(file="LD_result_glm_no_batch.csv", broom::tidy(mglm, exponentiate = TRUE, conf.int = TRUE))
nobs(mglm) #2946
```
