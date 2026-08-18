# Configuration for 100 Chondrichthyan trees

tree.file.paths <- list.files(
  "data_raw/shark_trees_100",
  pattern = "\\.tre$",
  full.names = TRUE
)

stopifnot(length(tree.file.paths) == 100)

results.file.path <- "outputs/100_tree_results/"
graphs.file.path <- "outputs/figures/"
csv.file.path <- "outputs/tables/"

dir.create(
  results.file.path,
  recursive = TRUE,
  showWarnings = FALSE
)
