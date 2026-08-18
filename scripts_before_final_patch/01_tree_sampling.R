# 01_tree_sampling.R
# Reproducibly sample 100 trees from the 10,000-tree VertLife distribution.

tree_set_file <- "data_raw/Chond.10Cal.10kTreeSet.tre"

if (!file.exists(tree_set_file)) {
  stop("Missing phylogenetic tree set: ", tree_set_file)
}

sampled_tree_dir <- "data_raw/shark_trees_100"

existing_sampled <- list.files(
  sampled_tree_dir,
  pattern = "^shark_tree_[0-9]{3}\\.tre$",
  full.names = TRUE
)

if (length(existing_sampled) == 100) {
  message("100 sampled tree files already exist — tree sampling skipped.")
} else {
  message("Reading the 10,000-tree distribution...")
  all_trees <- ape::read.tree(tree_set_file)

  if (length(all_trees) < 100) {
    stop("Tree distribution contains fewer than 100 trees.")
  }

  set.seed(123)
  selected_index <- sample(seq_along(all_trees), 100, replace = FALSE)
  selected_trees <- all_trees[selected_index]

  write.csv(
    data.frame(
      sampled_tree_id = seq_along(selected_index),
      original_tree_index = selected_index
    ),
    "outputs/tables/selected_tree_indices.csv",
    row.names = FALSE
  )

  for (i in seq_along(selected_trees)) {
    ape::write.tree(
      selected_trees[[i]],
      file = sprintf(
        "%s/shark_tree_%03d.tre",
        sampled_tree_dir,
        i
      )
    )
  }

  rm(all_trees, selected_trees)
  gc()
}

tree.file.paths <- list.files(
  sampled_tree_dir,
  pattern = "^shark_tree_[0-9]{3}\\.tre$",
  full.names = TRUE
)

tree.file.paths <- sort(tree.file.paths)

stopifnot(length(tree.file.paths) == 100)

message("Tree sampling step complete: ", length(tree.file.paths), " trees.")
