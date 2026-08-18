# 04_EvoHeritage_analysis.R
# Formal 100-tree, 10,000-repeat EvoHeritage analyses at 66 Ma and 145 Ma.
# Restart-safe: completed result files are skipped.

eh_output_dir <- "outputs/EvoHeritage_100trees_final"

run_eh_one_tree <- function(
  tree,
  min_age,
  rho = 0.01,
  repeats = 10000
) {

  eh_tree <- make.ancestral.evoheritage.tree(
    tree,
    rho = rho,
    lambda = 1,
    min.age = min_age,
    max.age = 4025
  )

  Partitioned.EvoHeritage(
    input.tree = eh_tree,
    num.repeats = repeats
  )
}

for (j in seq_along(tree.file.paths)) {

  message("EvoHeritage tree ", j, " of ", length(tree.file.paths))

  current_tree <- ape::read.tree(tree.file.paths[j])
  current_tree <- ape::reorder.phylo(current_tree, order = "cladewise")

  file66 <- sprintf(
    "%s/EH66_tree_%03d.csv",
    eh_output_dir,
    j
  )

  if (!file.exists(file66)) {
    message("  Running 66 Ma...")
    set.seed(66000 + j)

    result_66 <- run_eh_one_tree(
      current_tree,
      min_age = 66,
      rho = 0.01,
      repeats = 10000
    )

    write.csv(
      result_66,
      file66,
      row.names = FALSE
    )
  } else {
    message("  66 Ma already completed — skipping.")
  }

  file145 <- sprintf(
    "%s/EH145_tree_%03d.csv",
    eh_output_dir,
    j
  )

  if (!file.exists(file145)) {
    message("  Running 145 Ma...")
    set.seed(145000 + j)

    result_145 <- run_eh_one_tree(
      current_tree,
      min_age = 145,
      rho = 0.01,
      repeats = 10000
    )

    write.csv(
      result_145,
      file145,
      row.names = FALSE
    )
  } else {
    message("  145 Ma already completed — skipping.")
  }
}

n66 <- length(list.files(eh_output_dir, pattern = "^EH66_tree_[0-9]{3}\\.csv$"))
n145 <- length(list.files(eh_output_dir, pattern = "^EH145_tree_[0-9]{3}\\.csv$"))

if (n66 != 100 || n145 != 100) {
  stop(
    "EvoHeritage run incomplete. 66 Ma files: ", n66,
    "; 145 Ma files: ", n145
  )
}

message("Formal EvoHeritage analyses complete.")
