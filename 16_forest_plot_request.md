# Forest plot (collaborator request)

We are going to create a plot for  family 1: multinomial, archaeal dominance group as outcome (n = 2,841)

```R
library(ggplot2)
library(dplyr)

setwd("/Users/camilagazollavolpiano/Desktop/ARCHAEA/V6_corrections")
infile  <- "data_for_forrest_plot.csv"
outfile <- "fig_forest_family1.svg"

## -- 1. Read; drop terms not interpretable as effects on this plot ----------
d <- read.csv(infile, stringsAsFactors = FALSE)
d <- d %>% filter(term != "(Intercept)", !grepl("^batch_name", term))

## -- 2. Labels and display order (first listed = top of panel) --------------
cont_unit <- c(BL_AGE = "per year", BMI = "per kg/m\u00B2", FIBER_TOTAL = "per unit", RED_MEAT = "per unit")

term_map <- c(
  BMI             = sprintf("BMI (%s)", cont_unit["BMI"]),
  RED_MEAT        = sprintf("Red meat intake (%s)", cont_unit["RED_MEAT"]),
  FIBER_TOTAL     = sprintf("Total fibre intake (%s)", cont_unit["FIBER_TOTAL"]),
  CURR_SMOKE1     = "Current smoking",
  BL_USE_RX_J011  = "Antibiotic use",
  BL_AGE          = sprintf("Age (%s)", cont_unit["BL_AGE"]),
  MEN1            = "Male sex",
  ReadDepth_log10 = "Sequencing depth (log10)"
)

d <- d %>%
  mutate(
    term_lab = factor(unname(term_map[term]), levels = rev(unname(term_map))),
    y.level  = factor(y.level, levels = c("M. intestini", "M. sp900766745", "Diverse")),
    sig      = p.value < 0.05
  )

## Italic species names in strip labels
strip_lab <- as_labeller(c(
  "M. intestini"   = "italic('M. intestini')",
  "M. sp900766745" = "italic('M.')~'sp900766745'",
  "Diverse"        = "'Diverse'"
), label_parsed)

## -- 4. theme: ~7 pt sans, hairline rules, no gridlines --------
bs <- 7
theme_nature <- theme_classic(base_size = bs, base_family = "sans") +
  theme(
    axis.line.x       = element_line(linewidth = 0.25, colour = "black"),
    axis.line.y       = element_blank(),
    axis.ticks.x      = element_line(linewidth = 0.25, colour = "black"),
    axis.ticks.y      = element_blank(),
    axis.ticks.length = unit(1.5, "pt"),
    axis.text         = element_text(colour = "black", size = bs),
    axis.title        = element_text(colour = "black", size = bs),
    strip.background  = element_blank(),
    strip.text        = element_text(size = bs, hjust = 0),
    panel.spacing.x   = unit(6, "pt"),
    legend.position   = "none",
    plot.margin       = margin(4, 6, 4, 4)
  )

## Replace with the Figure 2 community-type palette for cross-figure consistency
pal <- c("M. intestini" = "#4C72B0", "M. sp900766745" = "#C44E52", "Diverse" = "#8C8C8C")

## -- 5. Plot ----------------------------------------------------------------
p <- ggplot(d, aes(x = estimate, y = term_lab, colour = y.level)) +
  geom_vline(xintercept = 1, linewidth = 0.25, linetype = "dashed", colour = "grey60") +
  geom_linerange(aes(xmin = conf.low, xmax = conf.high), linewidth = 0.4) +
  geom_point(aes(shape = sig), size = 1.4, fill = "white", stroke = 0.4) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21)) +
  scale_x_continuous(
    trans  = "log10",
    breaks = c(0.5, 0.7, 1, 1.5, 2, 3),
    labels = c("0.5", "0.7", "1", "1.5", "2", "3"),
    expand = expansion(mult = 0.04)
  ) +
  scale_colour_manual(values = pal) +
  facet_wrap(~ y.level, nrow = 1, labeller = strip_lab) +
  labs(x = expression("Odds ratio (95% CI) vs "*italic("M. smithii")*"-dominated"),
       y = NULL) +
  theme_nature

## -- 6. Write (svg / print / dev.off) ---------------------------------------
svg(outfile, width = 7.2, height = 2.3)   # 183 mm double-column
print(p)
dev.off()
cat("wrote", outfile, "\n")
```
