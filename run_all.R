### RUNS THE COMPLETE ANALYSIS ###
# each script is run in its own folder and R session, in the order of the analysis

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

scripts <- c("01_technical_qc/01_technical_qc.R",
             "02_preprocessing/02_preprocessing.R",
             "03_concentration_ranges/03_concentration_ranges.R",
             "04_covariate_effects/04_covariate_effects.R",
             "05_reference_intervals/05_reference_intervals.R",
             "06_method_robustness/06_method_robustness.R")              # takes about half an hour

rscript <- file.path(R.home("bin"), "Rscript")

for (script in scripts) {
  message("running ", script)
  owd <- setwd(dirname(script))
  status <- system2(rscript, shQuote(basename(script)))
  setwd(owd)
  if (status != 0) stop(script, " failed")
}
