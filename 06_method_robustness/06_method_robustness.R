### ROBUSTNESS OF THE RI CALCULATION METHODS ###
# RIs from 100 random subsets (n = 1200, drawn without replacement) are compared with the RIs of the full dataset
# (same method) by the Jaccard index; differences between methods are quantified with Cliff's delta
# outputs: Figure 5, Jaccard indices, Cliff's deltas

library(openxlsx)
library(tidyverse)
library(effsize)
library(RColorBrewer)

rm(list = ls())

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

source("../ri_methods.R")

# partitions with < 20 participants are removed, as they are not reported
full_data <- read.xlsx("../02_preprocessing/02_output/Untransformed IQR Data.xlsx") %>%
  add_partitions() %>%
  filter(!(
    (sex == "Female" & menopausal_status == "postmenopausal" & age_groups == "40-45") |
    (sex == "Female" & menopausal_status == "premenopausal"  & age_groups == "55-60") |
    (sex == "Female" & menopausal_status == "premenopausal"  & age_groups == "60-65")
  ))

metabolites <- colnames(full_data)[7:56]

partitions <- count_partitions(full_data)

all_RIs <- function(dat) {
  bind_rows(RIs_partitioned(dat, partitions, metabolites, limits_NP, "non-parametric"),
            RIs_partitioned(dat, partitions, metabolites, limits_BC, "Box-Cox"),
            RIs_LM(LM_constants(dat, metabolites), partitions))
}



#_______________________________________________________________________________
# full data RIs ("big") vs. subset RIs ("small")

RIs_big <- all_RIs(full_data) %>%
  rename(golden_lower_border = lower_border, golden_upper_border = upper_border)

set.seed(121)

jaccard_indices <- map_dfr(1:100, function(iteration) {

  subsampled_data <- full_data %>%
    sample_n(size = 1200, replace = FALSE)

  print(paste("iteration", iteration))

  all_RIs(subsampled_data) %>%
    inner_join(RIs_big, by = c("method", "partition", "metabolite")) %>%
    mutate(iteration = iteration,
           jaccard = (pmin(upper_border, golden_upper_border) - pmax(lower_border, golden_lower_border)) /
                     (pmax(upper_border, golden_upper_border) - pmin(lower_border, golden_lower_border)))
})

write.csv(jaccard_indices, "02_output/Jaccard Indices.csv", row.names = FALSE)



#_______________________________________________________________________________
# Figure 5

jaccard_boxplot <- jaccard_indices %>%
  mutate(method = ifelse(method == "lm", "linear model", method),
         method = factor(method, levels = c("non-parametric", "Box-Cox", "linear model")))

colors <- brewer.pal(9, "Greens")[c(3,5,7)]

ggplot(jaccard_boxplot, aes(y = jaccard, x = method, fill = method)) +
  geom_boxplot(outlier.size = 0.1) +
  scale_fill_manual(values = colors) +
  labs(x = "Reference interval calculation method", y = "Jaccard index") +
  ylim(c(0.5, 1)) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6))

ggsave("02_output/Figure 5 - Jaccard Indices.tif", width = 85, height = 85, units = "mm", dpi = 600)
ggsave("02_output/Figure 5 - Jaccard Indices.png", width = 10, height = 10, units = "cm")



#_______________________________________________________________________________
# Cliff's delta between methods

jaccard_stat <- jaccard_boxplot %>%
  pivot_wider(id_cols = c("partition", "metabolite", "iteration"), names_from = "method", values_from = "jaccard")

cliff <- function(a, b) {
  d <- cliff.delta(jaccard_stat[[a]], jaccard_stat[[b]])
  tibble(comparison = paste(a, "vs.", b), delta = d$estimate, ci_low = d$conf.int[[1]], ci_high = d$conf.int[[2]],
         magnitude = as.character(d$magnitude))
}

cliffs_deltas <- bind_rows(cliff("linear model", "non-parametric"),
                           cliff("linear model", "Box-Cox"),
                           cliff("Box-Cox", "non-parametric"))
print(cliffs_deltas)

write.xlsx(cliffs_deltas, "02_output/Cliffs Deltas.xlsx")
