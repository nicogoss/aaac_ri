### CONCENTRATION RANGES ###
# mean concentration of each metabolite in the total cohort, males, premenopausal and postmenopausal females
# outputs: Figure 1 (amino acids), Figure 2 (acylcarnitines)

library(openxlsx)
library(tidyverse)

rm(list = ls())

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

full_data    <- read.xlsx("../02_preprocessing/02_output/Untransformed IQR Data.xlsx")
abbreviation <- read.xlsx("../02_preprocessing/01_input/Metabolite Abbreviations.xlsx") %>%
  select(metabolite, abbreviation)

metabolites <- colnames(full_data)[5:54]

group_colors <- c("male" = "steelblue", "premenopausal" = "darkorange", "postmenopausal" = "firebrick")



#_______________________________________________________________________________
# average concentrations per stratum

long_data <- full_data %>%
  pivot_longer(all_of(metabolites), names_to = "metabolite", values_to = "conc")

total_average <- long_data %>%
  group_by(metabolite) %>%
  summarize(mean_conc = mean(conc, na.rm = TRUE), .groups = "drop") %>%
  mutate(stratum = "total")

strata_average <- long_data %>%
  mutate(menopausal_status = ifelse(menopausal_status == "Male", "male", menopausal_status)) %>%
  filter(!menopausal_status == "unresolved") %>%
  group_by(menopausal_status, metabolite) %>%
  summarize(mean_conc = mean(conc, na.rm = TRUE), .groups = "drop") %>%
  rename("stratum" = "menopausal_status")

average <- bind_rows(total_average, strata_average) %>%
  left_join(abbreviation, by = "metabolite") %>%
  mutate(class = ifelse(grepl("carnitine", metabolite, ignore.case = TRUE), "AC", "AA"))



#_______________________________________________________________________________
# plot: one bar per stratum, metabolites connected across strata

# spreads label positions by a minimum gap on the log axis while staying as close to the true values as possible
spread_labels <- function(x, min_gap) {
  ord <- order(x)
  steps <- (seq_along(x) - 1) * min_gap
  spread <- isoreg(log10(x[ord]) - steps)$yf + steps
  out <- numeric(length(x))
  out[ord] <- 10^spread
  out
}

plot_bars <- function(dat, x_title, x_breaks, label_gap) {

  plot_data <- dat %>%
    mutate(stratum = factor(stratum, levels = c("total", "male", "premenopausal", "postmenopausal")),
           y = as.numeric(stratum)) %>%
    arrange(y) # geom_path connects in data order

  ranges <- plot_data %>%
    group_by(stratum, y) %>%
    summarize(min_conc = min(mean_conc), max_conc = max(mean_conc), .groups = "drop")

  labels_top    <- filter(plot_data, stratum == "postmenopausal") %>%
    mutate(label_x = spread_labels(mean_conc, label_gap))
  labels_bottom <- filter(plot_data, stratum == "total") %>%
    mutate(label_x = spread_labels(mean_conc, label_gap))

  ggplot(plot_data, aes(x = mean_conc, y = y)) +
    geom_rect(data = ranges, aes(xmin = min_conc, xmax = max_conc, ymin = y - 0.07, ymax = y + 0.07, fill = stratum),
              alpha = 0.5, inherit.aes = FALSE) +
    scale_fill_manual(values = c("total" = "grey60", group_colors)) +
    geom_path(aes(group = metabolite), linetype = "longdash", color = "grey40", linewidth = 0.2, alpha = 0.7) +
    geom_point(shape = 124, size = 2.2, stroke = 0.8) +
    geom_segment(data = labels_top, aes(x = mean_conc, xend = label_x, y = y + 0.09, yend = y + 0.25),
                 linewidth = 0.15, color = "grey50") +
    geom_segment(data = labels_bottom, aes(x = mean_conc, xend = label_x, y = y - 0.09, yend = y - 0.25),
                 linewidth = 0.15, color = "grey50") +
    geom_text(data = labels_top, aes(x = label_x, y = y + 0.27, label = abbreviation), size = 2.3, angle = 90, hjust = 0) +
    geom_text(data = labels_bottom, aes(x = label_x, y = y - 0.27, label = abbreviation), size = 2.3, angle = 90, hjust = 1) +
    scale_x_log10(breaks = x_breaks, labels = scales::label_number(big.mark = ","),
                  minor_breaks = scales::minor_breaks_log()) +
    scale_y_continuous(breaks = 1:4, labels = levels(plot_data$stratum)) +
    coord_cartesian(ylim = c(0.15, 4.85), clip = "off") +
    labs(x = x_title, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(axis.text = element_text(size = 7, color = "black"),
          axis.title.x = element_text(size = 8),
          legend.position = "none",
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(color = "grey70", linetype = "dotted", linewidth = 0.3),
          panel.grid.minor.x = element_line(color = "grey85", linetype = "dotted", linewidth = 0.2))
}

bars_AA <- plot_bars(filter(average, class == "AA"), "AA concentration [µM]", c(2, 5, 10, 20, 50, 100, 200, 500), label_gap = 0.045)
bars_AC <- plot_bars(filter(average, class == "AC"), "AC concentration [nM]", c(20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000), label_gap = 0.055)

ggsave("02_output/Figure 1 - Concentration Ranges AA.tif", plot = bars_AA, height = 90, width = 176, unit = "mm", dpi = 600, compression = "lzw")
ggsave("02_output/Figure 2 - Concentration Ranges AC.tif", plot = bars_AC, height = 90, width = 176, unit = "mm", dpi = 600, compression = "lzw")
ggsave("02_output/Figure 1 - Concentration Ranges AA.png", plot = bars_AA, height = 90, width = 176, unit = "mm", dpi = 300)
ggsave("02_output/Figure 2 - Concentration Ranges AC.png", plot = bars_AC, height = 90, width = 176, unit = "mm", dpi = 300)
