### TECHNICAL QUALITY CONTROL ###
# batch-to-batch precision (pooled QC) and accuracy (NIST SRM 1950) of all targeted metabolites;
# metabolites with a pooled QC CV > 20% are excluded from all further analyses
# outputs: Figure S2, Figure S3, Figure S4, Table S3

library(openxlsx)
library(tidyverse)
library(patchwork)

rm(list = ls())

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

qc_file <- "01_input/QC NIST-Plasma.xlsx"

aa_not_detected <- c("gamma-aminobutyrate", "carnosine")
ac_not_detected <- c("hydroxyoctadecanoylcarnitine", "hydroxyhexadecanoylcarnitine", "hexenoylcarnitine", "sebacoylcarnitine",
                     "suberoylcarnitine", "hydroxyhexanoylcarnitine", "succinylcarnitine", "malonylcarnitine")



#_______________________________________________________________________________
# mean, sd and CV per metabolite across all repeated QC injections

read_qc <- function(sheet, qc_type, not_detected) {
  read.xlsx(qc_file, sheet = sheet) %>%
    filter(Filename != "Filename") %>%
    filter(str_detect(Filename, regex(qc_type, ignore_case = TRUE))) %>%
    select(!any_of(not_detected)) %>%
    t() %>%
    as.data.frame() %>%
    setNames(.[1,]) %>%
    slice(2:nrow(.)) %>%
    mutate_all(~ gsub("NF", NA, .)) %>%
    mutate_all(~ as.numeric(.)) %>%
    mutate(class = sheet,
           mean = rowMeans(across(starts_with("CoLaus")), na.rm = TRUE),
           sd = apply(across(starts_with("CoLaus")), 1, function(x) sd(x, na.rm = TRUE)),
           CV = sd/mean * 100) %>%
    rownames_to_column("metabolite") %>%
    arrange(CV) %>%
    select(metabolite, class, mean, sd, CV)
}

pool_full <- bind_rows(read_qc("AA", "pool", aa_not_detected),
                       read_qc("AC", "pool", ac_not_detected))

NIST_AA <- read_qc("AA", "NIST", aa_not_detected)
NIST_AC <- read_qc("AC", "NIST", ac_not_detected)



#_______________________________________________________________________________
# accuracy: measured vs. certified NIST concentrations (amino acids only)

AA_conc_NIST <- read.xlsx("01_input/NIST reference concentration.xlsx") %>%
  left_join(NIST_AA, by = "metabolite") %>%
  mutate(upper_measured = mean + 1.96*sd,
         lower_measured = mean - 1.96*sd) %>%
  mutate(upper_certified = NISTconc + NISTrange,
         lower_certified = NISTconc - NISTrange) %>%
  arrange(mean) %>%
  mutate(accuracy = mean / NISTconc * 100)

accuracy <- ggplot(AA_conc_NIST, mapping = aes(x = mean, y = NISTconc)) +
  geom_errorbar(aes(ymin = lower_certified, ymax = upper_certified), color = "gray") +
  geom_errorbarh(aes(xmin = lower_measured, xmax = upper_measured), color = "gray") +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_text(size = 2.5,
            aes(label = metabolite,
                vjust = c(+0.5, +1.00, +0.0, +0.5, +0.5, +0.50, +0.90, +0.5, +0.5, +0.5),       # label positions, in order of increasing concentration
                hjust = c(-0.2, -0.12, -0.3, +1.4, -0.3, -0.60, +1.45, -1.3, +1.6, +1.7))) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
  labs(x = "NIST - measured concentrations [µM]",
       y = "NIST - certified concentrations [µM]")



#_______________________________________________________________________________
# precision: pooled QC CVs and exclusion of metabolites with CV > 20%

exclude <- pool_full %>%
  filter(CV > 20)

pool_retained <- pool_full %>%
  filter(CV < 20)

precision <- ggplot(pool_retained, aes(x = class, y = CV)) +
  geom_violin(aes(fill = class), width = 1) +
  scale_fill_brewer(palette = "Pastel2") +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0)) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red") +
  labs(y = "CV [%]") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position = "none"
  ) +
  ylim(0, 25)

figure_S2 <- precision + accuracy +
  plot_annotation(tag_levels = "A")

ggsave("02_output/Figure S2 - Precision and Accuracy.tif", figure_S2, width = 176, height = 85, units = "mm", dpi = 600)



#_______________________________________________________________________________
# comparison of NIST and pooled QC CVs (retained metabolites only)

NIST_full <- bind_rows(NIST_AA, NIST_AC) %>%
  left_join(select(AA_conc_NIST, metabolite, NISTconc, NISTrange, accuracy), by = "metabolite") %>%
  filter(metabolite %in% pool_full$metabolite)

comparison <- pool_retained %>%
  mutate(qc = "pooled QC") %>%
  bind_rows(filter(NIST_full, metabolite %in% pool_retained$metabolite)) %>%
  mutate(qc = ifelse(is.na(qc), "NIST QC", qc),
         qc = as.factor(qc))

ggplot(data = comparison,
       mapping = aes(x = qc, y = CV)) +
  geom_violin(fill = "lightgray") +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(x = as.numeric(qc), y = CV), method = "lm", se = FALSE, linewidth = 2, color = "red") +
  geom_line(aes(group = metabolite), alpha = 0.3) +
  facet_wrap(~ class) +
  theme_light() +
  theme(axis.title.x = element_blank())
ggsave("02_output/Figure S3 - CVs NIST vs Pooled QC.png", plot = last_plot(), width = 18, height = 15, units = "cm")

# the CV difference between both QCs follows their concentration difference
relative_comparison <- comparison %>%
  pivot_wider(
    id_cols = c(metabolite, class),
    names_from = qc,
    values_from = c(mean, sd, CV),
    names_glue = "{qc}_{.value}"
  ) %>%
  mutate(conc_ratio = `NIST QC_mean` / `pooled QC_mean`,
         cv_ratio = `NIST QC_CV` / `pooled QC_CV`)

corr_df <- relative_comparison %>%
  group_by(class) %>%
  reframe(
    cor = cor(conc_ratio, cv_ratio, method = "pearson"),
    p   = cor.test(conc_ratio, cv_ratio, method = "pearson")$p.value
  ) %>%
  mutate(
    label = paste0("r = ", round(cor, 3),
                   "\np = ", signif(p, 3)),
    x = Inf,
    y = Inf
  )

ggplot(data = relative_comparison,
       mapping = aes(x = conc_ratio, y = cv_ratio)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "black") +
  facet_wrap(~ class) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_text(data = corr_df,
            aes(x = x, y = y, label = label),
            hjust = 1.1, vjust = 1.2,
            size = 4) +
  theme_light() +
  ylab("CV ratio (NIST/pool)") +
  xlab("concentration ratio (NIST/pool)")
ggsave("02_output/Figure S4 - Concentration Effect on CVs.png", plot = last_plot(), width = 18, height = 15, units = "cm")



#_______________________________________________________________________________
# export (Table S3: pooled QC and NIST QC CVs, NIST accuracy)

write.xlsx(pool_full, "02_output/Pooled QC CVs.xlsx")
write.xlsx(NIST_full, "02_output/NIST QC CVs.xlsx")
write.xlsx(exclude, "02_output/excluded metabolites.xlsx")
