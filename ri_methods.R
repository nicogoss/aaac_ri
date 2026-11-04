### REFERENCE INTERVAL METHODS ###
# shared by 05_reference_intervals and 06_method_robustness:
# non-parametric (NP) and Box-Cox (BC) approach per partition, linear model (LM) approach on the full data

library(MASS, include.only = "boxcox")
library(car, include.only = "powerTransform")
library(VGAM, include.only = "yeo.johnson")



#_______________________________________________________________________________
# partitions: 5-year age groups (65-80 as one group), sex, menopausal status

age_group_breaks  <- c(40, 45, 50, 55, 60, 65, 80)
age_group_labels  <- c("40-45", "45-50", "50-55", "55-60", "60-65", "65-80")
age_group_centers <- c(42.5, 47.5, 52.5, 57.5, 62.5, 72.5)                      # ages at which the LM-based RIs are reported

add_partitions <- function(dat) {
  dat %>%
    mutate(age_groups       = cut(age, breaks = age_group_breaks, labels = age_group_labels),
           age_group_center = cut(age, breaks = age_group_breaks, labels = age_group_centers)) %>%
    relocate(age_groups, age_group_center, .after = age) %>%
    mutate(sex = as.factor(sex),
           menopausal_status = as.factor(menopausal_status)) %>%
    filter(!menopausal_status == "unresolved")
}

count_partitions <- function(dat) {
  dat %>%
    group_by(age_groups, age_group_center, menopausal_status, sex) %>%
    summarize(n = n(), .groups = "drop") %>%
    arrange(age_groups, menopausal_status, sex)
}



#_______________________________________________________________________________
# RI limits of a single concentration vector

limits_NP <- function(x) {
  unname(quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
}

limits_BC <- function(x) {

  current_boxcox <- boxcox(x ~ 1, lambda = seq(-5, 5, 1/100), plotit = FALSE)  # find best lambda
  lambda <- current_boxcox$x[which.max(current_boxcox$y)]

  transformed <- if (lambda == 0) log(x + 1) else ((x^lambda) - 1)/lambda
  limits <- mean(transformed, na.rm = TRUE) + c(-1.96, 1.96) * sd(transformed, na.rm = TRUE)

  if (lambda == 0) exp(limits) - 1 else (lambda * limits + 1)^(1/lambda)         # back-transform
}



#_______________________________________________________________________________
# NP and BC approach: limits per partition and metabolite

RIs_partitioned <- function(dat, partitions, metabolites, limit_function, method) {

  pmap_dfr(partitions, function(sex, menopausal_status, age_groups, ...) {

    partition <- dat[dat$sex == sex & dat$menopausal_status == menopausal_status & dat$age_groups == age_groups, ]
    limits <- sapply(partition[metabolites], limit_function)

    tibble(method = method,
           partition = paste(sex, age_groups, menopausal_status),
           metabolite = metabolites,
           lower_border = limits[1,],
           upper_border = limits[2,])
  })
}



#_______________________________________________________________________________
# LM approach: model of age, sex, menopausal status and their interactions on untransformed concentrations;
# Yeo-Johnson-transformed residuals give the limits (mean +/- 1.96 SD), which are back-transformed and added to the prediction

LM_constants <- function(dat, metabolites) {

  map_dfr(metabolites, function(metabolite) {

    d <- dat %>% mutate(y = .data[[metabolite]])
    fit <- lm(y ~ age +
                (sex == "Female") +
                (menopausal_status == "postmenopausal") +
                (sex == "Female") * age +
                (menopausal_status == "postmenopausal") * age,
              data = d, na.action = "na.exclude")

    coefficients <- unname(coef(fit))
    residuals <- fit$residuals

    lambda <- powerTransform(residuals, family = "yjPower")$lambda[[1]]
    residuals_YJ <- yeo.johnson(residuals, lambda)
    limits_YJ <- mean(residuals_YJ) + c(-1.96, 1.96) * sd(residuals_YJ)

    tibble(metabolite = metabolite,
           Intercept = coefficients[1],
           coeff_age = coefficients[2],
           coeff_sex = coefficients[3],
           coeff_menopausal_status = coefficients[4],
           coeff_age_sex = coefficients[5],
           coeff_age_menopausal_status = coefficients[6],
           lambda_YJ = lambda,
           lower_border = yeo.johnson(limits_YJ[1], lambda, inverse = TRUE),
           upper_border = yeo.johnson(limits_YJ[2], lambda, inverse = TRUE))
  })
}

RIs_LM <- function(constants, partitions) {

  expand_grid(partitions %>% select(age_group_center, sex, menopausal_status), constants) %>%
    mutate(age = as.numeric(as.character(age_group_center)),
           female = sex == "Female",
           post = menopausal_status == "postmenopausal",
           prediction = Intercept + age * coeff_age + female * coeff_sex + post * coeff_menopausal_status +
             female * age * coeff_age_sex + post * age * coeff_age_menopausal_status,
           lower_border = round(prediction + lower_border, digits = 2),
           upper_border = round(prediction + upper_border, digits = 2),
           lower_border = ifelse(lower_border < 0, 0, lower_border),
           method = "lm",
           partition = paste(sex, cut(age, breaks = age_group_breaks, labels = age_group_labels), menopausal_status)) %>%
    select(method, partition, metabolite, lower_border, upper_border)
}
