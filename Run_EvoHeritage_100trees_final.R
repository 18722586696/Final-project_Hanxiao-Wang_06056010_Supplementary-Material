library(ape)

source("EvoHeritage/EvoHeritage Tools.R")
source("EvoHeritage/Living Fossils/Living Fossil Configuration.R")

dir.create(
  "outputs/EvoHeritage_100trees_final",
  recursive = TRUE,
  showWarnings = FALSE
)

run_eh_one_tree <- function(tree, min_age, rho = 0.01, repeats = 10000) {
  
  eh_tree <- make.ancestral.evoheritage.tree(
    tree,
    rho = rho,
    lambda = 1,
    min.age = min_age,
    max.age = 4025
  )
  
  result <- Partitioned.EvoHeritage(
    input.tree = eh_tree,
    num.repeats = repeats
  )
  
  result
}

for (j in seq_along(tree.file.paths)) {
  
  cat("\nTree", j, "of", length(tree.file.paths), "\n")
  
  current_tree <- ape::read.tree(tree.file.paths[j])
  
  # ---------- 66 Ma ----------
  file66 <- sprintf(
    "outputs/EvoHeritage_100trees_final/EH66_tree_%03d.csv",
    j
  )
  
  if (!file.exists(file66)) {
    
    cat("Running 66 Ma...\n")
    
    set.seed(66000 + j)
    
    result_66 <- run_eh_one_tree(
      tree = current_tree,
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
    
    cat("66 Ma already completed — skipping\n")
  }
  
  # ---------- 145 Ma ----------
  file145 <- sprintf(
    "outputs/EvoHeritage_100trees_final/EH145_tree_%03d.csv",
    j
  )
  
  if (!file.exists(file145)) {
    
    cat("Running 145 Ma...\n")
    
    set.seed(145000 + j)
    
    result_145 <- run_eh_one_tree(
      tree = current_tree,
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
    
    cat("145 Ma already completed — skipping\n")
  }
}