### DATA PRE-PROCESSING ###
# removes non-detected and imprecise metabolites, averages duplicate injections, sets IQR outliers to missing
# and merges the concentrations with age, sex and menopausal status
# output: analysis dataset used by all downstream scripts

library(openxlsx)
library(dplyr)

rm(list = ls())

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}



#_______________________________________________________________________________
# loading concentrations and metadata

aa_data <- read.xlsx("01_input/Resultat CoLaus AA-AC.xlsx", sheet = "AA", startRow = 2, colNames = TRUE) %>%
  rename("Sample" = "Filename")

ac_data <- read.xlsx("01_input/Resultat CoLaus AA-AC.xlsx", sheet = "AC", startRow = 2, colNames = TRUE) %>%
  rename("Sample" = "Filename")

metadata <- read.csv("01_input/colaus_metadata_harmonized.csv") %>%
  filter(fu == 1) %>%                                                           # first follow-up
  select(pt, sex, age, mp) %>%
  mutate(sex = ifelse(sex == 0, "Female", "Male")) %>%
  mutate(mp = ifelse(mp == "male", "Male", ifelse(mp == "pre", "premenopausal", ifelse(mp == "post", "postmenopausal", mp)))) %>%
  rename("Sample" = "pt",
         "menopausal_status" = "mp")



#_______________________________________________________________________________
# removing metabolites which were not detected or too imprecise (pooled QC CV > 20%, see 01_technical_qc)

excluded <- read.xlsx("../01_technical_qc/02_output/excluded metabolites.xlsx")

aa_not_detected <- c("gamma-aminobutyrate", "carnosine")
ac_not_detected <- c("hydroxyoctadecanoylcarnitine", "hydroxyhexadecanoylcarnitine", "hexenoylcarnitine", "sebacoylcarnitine",
                     "suberoylcarnitine", "hydroxyhexanoylcarnitine", "succinylcarnitine", "malonylcarnitine")

aa_found <- aa_data[, !names(aa_data) %in% aa_not_detected] %>%
  select(!intersect(colnames(.), excluded$metabolite))

ac_found <- ac_data[, !names(ac_data) %in% ac_not_detected] %>%
  select(!intersect(colnames(.), excluded$metabolite))



#_______________________________________________________________________________
# reducing sample names to participant IDs (e.g. "Colaus_AA&AC_0520220_100_P1_64" -> "100")

to_participant_id <- function(sample) {
  sample <- gsub("_P(.*)", "", sample)                                          # "Colaus_AA&AC_0520220_100_P1_64" -> "Colaus_AA&AC_0520220_100"
  sample <- gsub("-A(.*)", "", sample)                                          # "Colaus_AA&AC_0520220_100-AX" -> "Colaus_AA&AC_0520220_100" (duplicates)
  gsub("Colaus(.*)_", "", sample)                                               # "Colaus_AA&AC_0520220_100" -> "100"
}

aa_found$Sample <- to_participant_id(aa_found$Sample)
ac_found$Sample <- to_participant_id(ac_found$Sample)



#_______________________________________________________________________________
# averaging duplicate injections (duplicates follow directly on the first injection)

average_duplicates <- function(dat) {

  duplicates <- sort(which(duplicated(dat$Sample)), decreasing = TRUE)         # from bottom to top to avoid shifts in indices
  dat[] <- lapply(dat, as.numeric)

  for (i in duplicates) {
    dat[(i-1),] <- colMeans(rbind(dat[(i-1),], dat[i,]))
    dat <- dat[-i,]
  }
  dat
}

aa_clean <- average_duplicates(aa_found)
ac_clean <- average_duplicates(ac_found)

ac_clean$carnitine <- ac_clean$carnitine * 1000                                 # µM -> nM, so that all ACs are in nM



#_______________________________________________________________________________
# merging AAs and ACs, setting outliers (outside Q1 - 1.5 IQR, Q3 + 1.5 IQR) to missing

aaac_data <- merge(aa_clean, ac_clean, by = "Sample", all = TRUE)
aaac_data$Sample <- as.character(aaac_data$Sample)

aaac_clean_IQR <- aaac_data[aaac_data$Sample != 5106,]                          # sample 5106: almost no metabolites measured

for (i in 2:ncol(aaac_clean_IQR)) {

  x <- aaac_clean_IQR[,i]
  Q1 <- quantile(x, .25, na.rm = TRUE)[["25%"]]
  Q3 <- quantile(x, .75, na.rm = TRUE)[["75%"]]
  IQR <- IQR(x, na.rm = TRUE)

  aaac_clean_IQR[which(x > Q3 + 1.5*IQR | x < Q1 - 1.5*IQR), i] <- NA
}



#_______________________________________________________________________________
# merging with metadata and export

aaac_raw <- merge(metadata, aaac_clean_IQR, by = "Sample")

write.xlsx(aaac_raw, "02_output/Untransformed IQR Data.xlsx", rowNames = FALSE)
