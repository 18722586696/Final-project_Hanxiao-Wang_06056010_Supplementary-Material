# 02_ED_analysis.R
# Calculate ED independently on each of the 100 complete phylogenetic trees,
# then summarise species-level scores and ranks across trees.

ed_output_dir <- "outputs/100_tree_results"

for (j in seq_along(tree.file.paths)) {

  out_file <- sprintf(
    "%s/LivingFossil%04dED.csv",
    ed_output_dir,
    j
  )

  if (file.exists(out_file)) {
    message("ED tree ", j, " already completed — skipping.")
    next
  }

  message("Calculating ED: tree ", j, " of ", length(tree.file.paths))

  current_tree <- ape::read.tree(tree.file.paths[j])
  current_tree <- ape::reorder.phylo(current_tree, order = "cladewise")

  result <- Calculate_ED(input.tree = current_tree)

  write.csv(
    result,
    out_file,
    row.names = FALSE
  )
}

ed_files <- sort(list.files(
  ed_output_dir,
  pattern = "^LivingFossil[0-9]{4}ED\\.csv$",
  full.names = TRUE
))

if (length(ed_files) != 100) {
  stop("Expected 100 ED result files, found ", length(ed_files))
}

all_ed <- lapply(seq_along(ed_files), function(i) {
  x <- read.csv(ed_files[i], check.names = FALSE)
  x$tree_id <- i
  x
}) |>
  dplyr::bind_rows()

required_ed_cols <- c("tip.label", "partitioned.EvoHeritage", "lf.rank")
if (!all(required_ed_cols %in% names(all_ed))) {
  stop(
    "ED output is missing required columns: ",
    paste(setdiff(required_ed_cols, names(all_ed)), collapse = ", ")
  )
}

ed_summary_100 <- all_ed |>
  dplyr::group_by(tip.label) |>
  dplyr::summarise(
    mean_ED = mean(partitioned.EvoHeritage, na.rm = TRUE),
    median_ED = median(partitioned.EvoHeritage, na.rm = TRUE),
    sd_ED = sd(partitioned.EvoHeritage, na.rm = TRUE),
    mean_rank = mean(lf.rank, na.rm = TRUE),
    median_rank = median(lf.rank, na.rm = TRUE),
    top20_count = sum(lf.rank <= 20, na.rm = TRUE),
    top20_frequency = top20_count / 100,
    .groups = "drop"
  ) |>
  dplyr::arrange(mean_rank) |>
  dplyr::mutate(species = normalise_species(tip.label))

top20_ED_final <- ed_summary_100 |>
  dplyr::arrange(mean_rank) |>
  dplyr::slice(1:20)

write.csv(
  ed_summary_100,
  "outputs/tables/ED_100tree_species_summary.csv",
  row.names = FALSE
)

write.csv(
  top20_ED_final,
  "outputs/tables/ED_100tree_top20.csv",
  row.names = FALSE
)

message("ED analysis and 100-tree summary complete.")
