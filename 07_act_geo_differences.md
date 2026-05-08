# Tests claim that geography affect ACT

We can test this by looking at the distribution of ACTs across regions. We can use a chi-squared test to see if there is an association between region and ACT, and then look at the standardized residuals to see which ACTs drive the association. Finally, we can calculate Cramer's V to measure the strength of the association.

However, we should be cautious in interpreting these results, as there may be confounding factors.

```R
# libraries
library(tidyverse)

df <- read.csv("proportions_FR02_CMD_13716_ACTs.csv")

region_map <- c(
  "JieZ"="Asia","Kailuan"="Asia","QinN"="Asia","YachidaS"="Asia",
  "HMP2"="North America","MLVS"="North America",
  "500FG"="Europe","BackhedF"="Europe","BBS"="Europe","DIABIMMUNE"="Europe",
  "FINRISK_2002"="Europe","KarlssonFH"="Europe","MetaCardis"="Europe","PREDICT_1"="Europe"
)

df <- df %>% mutate(region = region_map[study])

tab <- df %>%
  group_by(region, archaeal_community_type) %>%
  summarise(n = sum(count), .groups="drop") %>%
  pivot_wider(names_from = archaeal_community_type, values_from = n, values_fill = 0) %>%
  tibble::column_to_rownames("region") %>%
  as.matrix()

tab
#               Diverse Methano_intestinii Methano_smithii Methano_sp900766745
# Asia              140                 26             278                 302
# Europe           1011                838            2846                1389
# North America      42                 23             172                  46

chisq.test(tab) # global association
#         Pearson's Chi-squared test

# data:  tab
# X-squared = 182.43, df = 6, p-value < 2.2e-16

chisq.test(tab)$stdres  # which ACTs drive it
#                  Diverse Methano_intestinii Methano_smithii Methano_sp900766745
# Asia           1.5412267          -7.851030       -5.252218           10.793561
# Europe        -0.8494090           8.092450        1.812582           -7.588433
# North America -0.8873395          -2.256718        4.971247           -3.263050

cs <- chisq.test(tab)
N <- sum(tab)
k <- min(nrow(tab)-1, ncol(tab)-1)
V <- sqrt(unname(cs$statistic) / (N * k))
V # [1] 0.1132424
```