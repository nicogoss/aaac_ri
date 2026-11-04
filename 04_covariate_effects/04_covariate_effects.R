### EFFECTS OF AGE, SEX AND MENOPAUSAL STATUS ###
# (1) effects of age, sex, menopause and their interactions on log-transformed concentrations of all 50 metabolites
# (2) calibration of the LM-based reference intervals in 5-fold cross-validation
# outputs: Figure 3, Figure 4, Figure S5, Table S4, Table S5, age distribution of pre- and postmenopausal females

library(openxlsx)
library(tidyverse)
library(car)
library(VGAM)
library(patchwork)

rm(list = ls())

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

data <- read.xlsx("../02_preprocessing/02_output/Untransformed IQR Data.xlsx") %>%
  filter(!menopausal_status == "unresolved")

abbreviation <- read.xlsx("../02_preprocessing/01_input/Metabolite Abbreviations.xlsx") %>%
  select(metabolite, abbreviation)

metabolites <- names(data)[which(names(data) == "creatinine"):ncol(data)]

age_group_breaks <- c(40, 45, 50, 55, 60, 65, 80)
age_group_labels <- c("40-45", "45-50", "50-55", "55-60", "60-65", "65-80")



#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||(0) AGE OVERLAP OF PRE- AND POSTMENOPAUSAL FEMALES||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

age_summary <- data %>%
  group_by(menopausal_status) %>%
  summarize(n = n(), min_age = min(age), q1 = quantile(age, 0.25), median = median(age),
            q3 = quantile(age, 0.75), max_age = max(age), .groups = "drop")

overlap_window <- data %>%
  filter(sex == "Female") %>%
  mutate(in_window = age >= 45 & age < 55) %>%
  count(menopausal_status, in_window)



#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||(1) EFFECTS AS PERCENT DIFFERENCES||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# metabolites log-transformed, age in years; effects reported as % difference (100 * (exp(beta) - 1)),
# age and interaction effects per decade
# female = 1 for all females, post = 1 for postmenopausal females (reference: males and premenopausal females), so that
# sex compares premenopausal females with males and menopause compares postmenopausal with premenopausal females
# main effects from model 1 (no interactions), age:sex from model 2, age:menopause from model 3

stopifnot(all(unlist(data[metabolites]) > 0, na.rm = TRUE))

data_scaled <- data %>%
  mutate(across(all_of(metabolites), log))

get_term <- function(fit, term) {
  ct <- summary(fit)$coefficients
  ci <- confint(fit)
  tibble(estimate = ct[term, "Estimate"], ci_low = ci[term, 1], ci_high = ci[term, 2],
         p = ct[term, "Pr(>|t|)"])
}

effects <- map_dfr(metabolites, function(metabolite) {

  d <- data_scaled %>% mutate(y = .data[[metabolite]], female = sex == "Female", post = menopausal_status == "postmenopausal")

  fit_main     <- lm(y ~ age + female + post, data = d)                         # model 1
  fit_age_sex  <- lm(y ~ female * age, data = d)                                # model 2
  fit_age_meno <- lm(y ~ female * age + post * age, data = d)                   # model 3

  bind_rows(
    get_term(fit_main, "age")                 %>% mutate(predictor = "age", across(estimate:ci_high, ~ . * 10)),
    get_term(fit_main, "femaleTRUE")          %>% mutate(predictor = "sex"),
    get_term(fit_main, "postTRUE")            %>% mutate(predictor = "menopause"),
    get_term(fit_age_sex, "femaleTRUE:age")   %>% mutate(predictor = "age:sex", across(estimate:ci_high, ~ . * 10)),
    get_term(fit_age_meno, "age:postTRUE")    %>% mutate(predictor = "age:menopause", across(estimate:ci_high, ~ . * 10))
  ) %>%
    mutate(metabolite = metabolite)
}) %>%
  rename(log_estimate = estimate, log_ci_low = ci_low, log_ci_high = ci_high) %>%
  mutate(estimate = 100 * (exp(log_estimate) - 1),
         ci_low   = 100 * (exp(log_ci_low) - 1),
         ci_high  = 100 * (exp(log_ci_high) - 1),
         p_bonferroni = pmin(p * length(metabolites), 1),
         significant = p_bonferroni < 0.05,
         class = ifelse(grepl("carnitine", metabolite, ignore.case = TRUE), "Acylcarnitines", "Amino acids"),
         predictor = factor(predictor, levels = c("age", "sex", "menopause", "age:sex", "age:menopause")))

effect_counts <- effects %>%
  group_by(predictor, class) %>%
  summarize(n_significant = sum(significant), n_up = sum(significant & estimate > 0),
            n_down = sum(significant & estimate < 0), .groups = "drop")

# Figure 3: heatmap
metabolite_order <- effects %>%
  filter(predictor == "age") %>%
  arrange(class, estimate) %>%
  pull(metabolite)

# colors on the log scale (so that e.g. -50% and +100% are equally intense), legend labelled in %
limit <- max(abs(effects$log_estimate))

ggplot(effects %>% mutate(metabolite = factor(metabolite, levels = metabolite_order)),
       aes(x = predictor, y = metabolite, fill = log_estimate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(significant, "*", "")), vjust = 0.75, size = 4) +
  facet_grid(class ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-limit, limit),
                       breaks = log(c(0.5, 0.75, 1, 1.25, 1.5)),
                       labels = function(x) sprintf("%+.0f%%", 100 * (exp(x) - 1)),
                       name = "Difference") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

ggsave("02_output/Figure 3 - Effect Heatmap.tif", plot = last_plot(), height = 240, width = 120, units = "mm", dpi = 600, compression = "lzw")
ggsave("02_output/Figure 3 - Effect Heatmap.png", plot = last_plot(), height = 240, width = 120, units = "mm", dpi = 300)

# Table S4: % difference (95% CI) and Bonferroni-adjusted p-values
table_effects <- effects %>%
  mutate(effect = sprintf("%.1f (%.1f, %.1f)", estimate, ci_low, ci_high),
         p_adj = ifelse(significant, formatC(p_bonferroni, format = "e", digits = 2), "n.s.")) %>%
  select(class, metabolite, predictor, effect, p_adj) %>%
  pivot_wider(names_from = predictor, values_from = c(effect, p_adj), names_vary = "slowest") %>%
  left_join(abbreviation, by = "metabolite") %>%
  mutate(class = ifelse(class == "Amino acids", "amino acid", "acylcarnitine")) %>%
  relocate(abbreviation, .after = class) %>%
  arrange(class, metabolite)



#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||(1b) FIGURE 4: AGE TRAJECTORIES||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

trajectory_data <- data %>%
  mutate(group = ifelse(sex == "Male", "male", menopausal_status),
         group = factor(group, levels = c("premenopausal", "postmenopausal", "male")))

plot_trajectory <- function(metabolite, y_title) {
  ggplot(trajectory_data, aes(x = age, y = .data[[metabolite]])) +
    geom_point(size = 0.5, color = "grey") +
    stat_smooth(method = "lm", formula = y ~ x, color = "#A50F15", fill = "#FB6A4A", linewidth = 0.5) +
    facet_wrap(~ group, strip.position = "bottom") +
    labs(y = y_title, x = "age") +
    theme_light()
}

figure_4 <- plot_trajectory("ornithine", "Ornithine [µM]") / plot_trajectory("carnitine", "Free Carnitine [nM]") +
  plot_annotation(tag_levels = "A")

ggsave("02_output/Figure 4 - Ornithine and Carnitine.tif", plot = figure_4, height = 220, width = 176, units = "mm", dpi = 600, compression = "lzw")
ggsave("02_output/Figure 4 - Ornithine and Carnitine.png", plot = figure_4, height = 220, width = 176, units = "mm", dpi = 300)



#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||(2) CALIBRATION OF LM-BASED REFERENCE INTERVALS|||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# same model and residual handling as the LM approach (see ri_methods.R)

full_data <- data %>%
  mutate(age_groups = cut(age, breaks = age_group_breaks, labels = age_group_labels),
         stratum = ifelse(sex == "Male", "Male", menopausal_status))

lm_ri_fit <- function(train, metabolite) {

  fit <- lm(get(metabolite) ~ age + (sex == "Female") + (menopausal_status == "postmenopausal") +
              (sex == "Female") * age + (menopausal_status == "postmenopausal") * age,
            data = train, na.action = "na.exclude")
  res <- na.omit(residuals(fit))
  lambda <- powerTransform(res, family = "yjPower")$lambda
  res_yj <- yeo.johnson(res, lambda)
  list(fit = fit,
       lower = yeo.johnson(mean(res_yj) - 1.96 * sd(res_yj), lambda, inverse = TRUE),
       upper = yeo.johnson(mean(res_yj) + 1.96 * sd(res_yj), lambda, inverse = TRUE))
}

lm_ri_flag <- function(ri, test, metabolite) {

  pred <- predict(ri$fit, newdata = test)
  obs <- test[[metabolite]]
  tibble(Sample = test$Sample, stratum = test$stratum, age_groups = test$age_groups, metabolite = metabolite,
         below = obs < pred + ri$lower, above = obs > pred + ri$upper) %>%
    filter(!is.na(obs))
}

# in-sample
coverage_in <- map_dfr(metabolites, function(metabolite) {
  lm_ri_flag(lm_ri_fit(full_data, metabolite), full_data, metabolite)
}) %>% mutate(type = "in-sample")

# 5-fold cross-validation
set.seed(121)
folds <- sample(rep(1:5, length.out = nrow(full_data)))

coverage_cv <- map_dfr(1:5, function(k) {
  map_dfr(metabolites, function(metabolite) {
    lm_ri_flag(lm_ri_fit(full_data[folds != k, ], metabolite), full_data[folds == k, ], metabolite)
  })
}) %>% mutate(type = "5-fold CV")

coverage <- bind_rows(coverage_in, coverage_cv)

coverage_overall <- coverage %>%
  group_by(type) %>%
  summarize(pct_below = 100 * mean(below), pct_above = 100 * mean(above),
            pct_outside = 100 * mean(below | above), .groups = "drop")

coverage_stratum <- coverage %>%
  group_by(type, stratum) %>%
  summarize(n_obs = n(), pct_below = 100 * mean(below), pct_above = 100 * mean(above),
            pct_outside = 100 * mean(below | above), .groups = "drop")

coverage_stratum_age <- coverage %>%
  group_by(type, stratum, age_groups) %>%
  summarize(n_obs = n(), pct_below = 100 * mean(below), pct_above = 100 * mean(above),
            pct_outside = 100 * mean(below | above), .groups = "drop")

# Table S5: cross-validated calibration per metabolite and group, with a summary row over all metabolites
calibration_summary <- function(dat) {
  dat %>%
    summarize(outside = 100 * mean(below | above), below = 100 * mean(below), above = 100 * mean(above),
              .groups = "drop") %>%
    relocate(outside, .after = above)
}

coverage_cv_groups <- bind_rows(coverage_cv %>% mutate(stratum = "all"), coverage_cv) %>%
  mutate(stratum = factor(stratum, levels = c("all", "Male", "premenopausal", "postmenopausal")))

table_calibration <- bind_rows(
  coverage_cv_groups %>% group_by(stratum) %>% calibration_summary() %>% mutate(metabolite = "all metabolites"),
  coverage_cv_groups %>% group_by(metabolite, stratum) %>% calibration_summary()
) %>%
  mutate(across(c(below, above, outside), ~ round(., 2))) %>%
  pivot_wider(names_from = stratum, values_from = c(below, above, outside), names_vary = "slowest") %>%
  left_join(abbreviation, by = "metabolite") %>%
  mutate(class = case_when(metabolite == "all metabolites" ~ "",
                           grepl("carnitine", metabolite, ignore.case = TRUE) ~ "acylcarnitine",
                           TRUE ~ "amino acid")) %>%
  relocate(class, abbreviation, .before = metabolite) %>%
  arrange(metabolite != "all metabolites", class, metabolite)

# Figure S5: calibration per age group
ggplot(coverage_stratum_age %>% filter(type == "5-fold CV"),
       aes(x = age_groups, y = pct_outside, color = stratum, group = stratum)) +
  geom_hline(yintercept = 5, linetype = "dashed") +
  geom_point(aes(size = n_obs)) +
  geom_line() +
  scale_color_manual(values = c("Male" = "#4575B4", "premenopausal" = "#FDAE61", "postmenopausal" = "#D73027")) +
  labs(x = "Age group [years]", y = "Observations outside LM-based RI [%]", size = "n (obs.)", color = NULL) +
  theme_light()

ggsave("02_output/Figure S5 - LM RI Calibration.png", plot = last_plot(), height = 100, width = 176, units = "mm", dpi = 300)



#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||EXPORT||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

write.xlsx(list("Age Summary" = age_summary,
                "Overlap 45-55" = overlap_window,
                "Effects" = effects,
                "Effect Counts" = effect_counts,
                "Coverage Overall" = coverage_overall,
                "Coverage Stratum" = coverage_stratum,
                "Coverage Stratum Age" = coverage_stratum_age,
                "Table S4 Effects" = table_effects,
                "Table S5 Calibration" = table_calibration),
           "02_output/Effects and LM RI Calibration.xlsx")
