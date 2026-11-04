### DEPLAOYS THE SHINY APP ###

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

rsconnect::deployApp("07_shiny_app")
