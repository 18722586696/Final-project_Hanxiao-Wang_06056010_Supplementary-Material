# 00_setup.R
# Project setup and reproducibility checks.

required_packages <- c("ape", "caper", "dplyr", "ggplot2")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install required packages before continuing: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(ape)
library(dplyr)
library(ggplot2)

required_dirs <- c(
  "data_raw",
  "data_raw/shark_trees_100",
  "outputs/100_tree_results",
  "outputs/EvoHeritage_100trees_final",
  "outputs/tables",
  "outputs/figures",
  "outputs/logs",
  "EvoHeritage"
)

invisible(lapply(
  required_dirs,
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

tools_file <- "EvoHeritage/EvoHeritage Tools.R"

if (!file.exists(tools_file)) {
  stop(
    "Missing ", tools_file, ". ",
    "Copy the original EvoHeritage Tools.R into the EvoHeritage folder."
  )
}

source(tools_file)

normalise_species <- function(x) {
  x <- gsub("_", " ", as.character(x), fixed = TRUE)
  trimws(gsub("[[:space:]]+", " ", x))
}

capture.output(
  sessionInfo(),
  file = "outputs/logs/R_session_info_start.txt"
)

message("Setup complete.")
