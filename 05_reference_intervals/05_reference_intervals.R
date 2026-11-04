### REFERENCE INTERVALS ###
# RIs of all metabolites with the non-parametric (NP), Box-Cox (BC) and linear model (LM) approach;
# LM-based RIs are reported at the center of each age group
# outputs: Table S1 (partition counts), Table S6 (LM), Table S7 (NP), Table S8 (BC), LM constants for the Shiny app

library(openxlsx)
library(tidyverse)

rm(list = ls())

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

source("../ri_methods.R")

full_data <- read.xlsx("../02_preprocessing/02_output/Untransformed IQR Data.xlsx") %>%
  add_partitions()

metabolites <- colnames(full_data)[7:56]

strata_count <- count_partitions(full_data)



#_______________________________________________________________________________
# RI calculation

NP_summary <- RIs_partitioned(full_data, strata_count, metabolites, limits_NP, "non-parametric")
BC_summary <- RIs_partitioned(full_data, strata_count, metabolites, limits_BC, "Box-Cox")

LM_constants_full <- LM_constants(full_data, metabolites)
LM_summary <- RIs_LM(LM_constants_full, strata_count)



#_______________________________________________________________________________
# export: one column per partition; partitions with < 20 participants are not reported

small_partitions <- c("Female 40-45 postmenopausal", "Female 55-60 premenopausal", "Female 60-65 premenopausal")

to_table <- function(summary) {
  summary %>%
    mutate(lower_border = round(lower_border, digits = 2),
           upper_border = round(upper_border, digits = 2)) %>%
    separate(partition, into = c("sex", "age_group", "menopausal_status"), sep = " ", remove = FALSE) %>%
    arrange(age_group, match(menopausal_status, c("Male", "premenopausal", "postmenopausal"))) %>%
    mutate(reference_interval = paste(lower_border, "-", upper_border)) %>%
    select(metabolite, partition, reference_interval) %>%
    pivot_wider(names_from = partition, values_from = reference_interval) %>%
    mutate(class = ifelse(grepl("carnitine", metabolite, ignore.case = TRUE), "acylcarnitine", "amino acid"),
           unit  = ifelse(class == "acylcarnitine", "nM", "µM")) %>%
    relocate(class, unit, .before = metabolite) %>%
    arrange(class, metabolite) %>%
    select(!all_of(small_partitions))
}

write.xlsx(strata_count, "02_output/Table S1 - Partition Counts.xlsx")
write.xlsx(to_table(LM_summary), "02_output/Table S6 - Linear Model RIs.xlsx")
write.xlsx(to_table(NP_summary), "02_output/Table S7 - Non-Parametric RIs.xlsx")
write.xlsx(to_table(BC_summary), "02_output/Table S8 - Box-Cox RIs.xlsx")

# coefficients and residual limits of the LM approach, also read by the Shiny app
LM_constants_full <- select(LM_constants_full, !lambda_YJ)
write.xlsx(LM_constants_full, "02_output/LM Constants.xlsx")
write.csv(LM_constants_full, "../07_shiny_app/lm_constants.csv", row.names = FALSE)
