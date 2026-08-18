# 03_trait_analysis.R
# Match the final 100-tree ED summary to the imputed trait dataset.
# Phylogenetic metrics are calculated before trait matching.

trait_file <- "data_raw/shark_missforest (1).RDA"

if (!file.exists(trait_file)) {
  stop("Missing trait data: ", trait_file)
}

trait_env <- new.env(parent = emptyenv())
loaded_names <- load(trait_file, envir = trait_env)

if (!"trait.mat" %in% loaded_names) {
  stop(
    "The trait RDA must contain an object called trait.mat. Found: ",
    paste(loaded_names, collapse = ", ")
  )
}

trait.mat <- trait_env$trait.mat

if (!"species" %in% names(trait.mat)) {
  if (is.null(rownames(trait.mat))) {
    stop("trait.mat contains neither a species column nor row names.")
  }
  trait.mat$species <- rownames(trait.mat)
}

trait.mat$species <- normalise_species(trait.mat$species)

merged_100tree <- ed_summary_100 |>
  dplyr::left_join(trait.mat, by = "species") |>
  dplyr::mutate(
    top20_ED_100tree = ifelse(
      species %in% top20_ED_final$species,
      "Top 20 ED",
      "Other species"
    )
  )

stopifnot(sum(merged_100tree$top20_ED_100tree == "Top 20 ED") == 20)

continuous_traits <- c(
  "length",
  "trophic.level",
  "max.depth",
  "growth.ratio",
  "shape.pc1",
  "shape.pc2"
)

continuous_traits <- continuous_traits[
  continuous_traits %in% names(merged_100tree)
]

wilcox_results <- lapply(continuous_traits, function(trait) {

  dat <- merged_100tree[
    !is.na(merged_100tree[[trait]]) &
      !is.na(merged_100tree$top20_ED_100tree),
  ]

  test <- wilcox.test(
    dat[[trait]] ~ dat$top20_ED_100tree,
    exact = FALSE
  )

  data.frame(
    trait = trait,
    W = unname(test$statistic),
    p_value = test$p.value,
    n_top20 = sum(
      dat$top20_ED_100tree == "Top 20 ED"
    ),
    n_other = sum(
      dat$top20_ED_100tree == "Other species"
    )
  )
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    p_adjusted_BH = p.adjust(p_value, method = "BH")
  )

write.csv(
  wilcox_results,
  "outputs/tables/trait_continuous_wilcoxon_BH.csv",
  row.names = FALSE
)

trait_summary <- merged_100tree |>
  dplyr::group_by(top20_ED_100tree) |>
  dplyr::summarise(
    n_species = dplyr::n(),
    mean_ED = mean(mean_ED, na.rm = TRUE),
    median_ED = median(mean_ED, na.rm = TRUE),
    dplyr::across(
      dplyr::all_of(continuous_traits),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        median = ~median(.x, na.rm = TRUE)
      )
    ),
    .groups = "drop"
  )

write.csv(
  trait_summary,
  "outputs/tables/trait_group_summary.csv",
  row.names = FALSE
)

categorical_traits <- c(
  "reproductive.guild",
  "vertical.position"
)

categorical_traits <- categorical_traits[
  categorical_traits %in% names(merged_100tree)
]

categorical_results <- list()

for (trait in categorical_traits) {

  tab <- table(
    merged_100tree$top20_ED_100tree,
    merged_100tree[[trait]],
    useNA = "no"
  )

  ft <- fisher.test(tab)

  categorical_results[[trait]] <- data.frame(
    trait = trait,
    p_value = ft$p.value
  )

  write.csv(
    as.data.frame.matrix(tab),
    sprintf(
      "outputs/tables/%s_contingency_table.csv",
      gsub("\\.", "_", trait)
    )
  )
}

if (length(categorical_results) > 0) {
  categorical_test_table <- dplyr::bind_rows(categorical_results)

  write.csv(
    categorical_test_table,
    "outputs/tables/trait_categorical_fisher.csv",
    row.names = FALSE
  )
}

write.csv(
  merged_100tree,
  "outputs/tables/ED_100tree_with_traits.csv",
  row.names = FALSE
)

message("Trait analysis complete.")
