# 03b_trait_sensitivity_observed.R
# Sensitivity analysis using original observed (non-imputed) shark traits.

observed_file <- "data_raw/shark_db.RDA"
if (!file.exists(observed_file)) stop("Missing observed trait data: ", observed_file)

obs_env <- new.env(parent = emptyenv())
obs_names <- load(observed_file, envir = obs_env)
if (!"shark.db" %in% obs_names) stop("shark_db.RDA must contain shark.db")
shark.db <- obs_env$shark.db

shark.db$species <- normalise_species(shark.db$species)
if (anyDuplicated(shark.db$species)) stop("Duplicated species names in shark.db")

# Convert the two columns known to be stored as character, while checking that
# conversion does not create new missing values from non-missing strings.
for (v in c("trophic.level", "max.depth")) {
  old <- shark.db[[v]]
  new <- suppressWarnings(as.numeric(as.character(old)))
  bad <- !is.na(old) & nzchar(trimws(as.character(old))) & is.na(new)
  if (any(bad)) stop("Non-numeric values found in ", v, ": ", paste(unique(old[bad]), collapse = ", "))
  shark.db[[v]] <- new
}

observed_merged <- ed_summary_100 |>
  dplyr::left_join(shark.db, by = "species") |>
  dplyr::mutate(
    top20_ED_100tree = ifelse(species %in% top20_ED_final$species,
                              "Top 20 ED", "Other species")
  )

matching_summary <- data.frame(
  metric = c("original_records", "matched_tree_species", "ED_top20_with_observed_record"),
  n = c(nrow(shark.db),
        sum(!is.na(observed_merged$spec.code)),
        sum(observed_merged$top20_ED_100tree == "Top 20 ED" & !is.na(observed_merged$spec.code)))
)
write.csv(matching_summary, "outputs/tables/trait_observed_matching_summary.csv", row.names = FALSE)

continuous_traits_obs <- c("length", "trophic.level", "max.depth", "growth.ratio", "shape.pc1", "shape.pc2")
observed_wilcox_results <- lapply(continuous_traits_obs, function(trait) {
  dat <- observed_merged[!is.na(observed_merged[[trait]]) & !is.na(observed_merged$top20_ED_100tree), ]
  test <- wilcox.test(dat[[trait]] ~ dat$top20_ED_100tree, exact = FALSE)
  data.frame(
    trait = trait,
    W = unname(test$statistic),
    p_value = test$p.value,
    n_top20 = sum(dat$top20_ED_100tree == "Top 20 ED"),
    n_other = sum(dat$top20_ED_100tree == "Other species")
  )
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    p_adjusted_BH = p.adjust(p_value, method = "BH"),
    significant_BH = p_adjusted_BH < 0.05
  )
write.csv(observed_wilcox_results, "outputs/tables/trait_observed_wilcoxon_BH.csv", row.names = FALSE)

observed_fisher_results <- lapply(c("reproductive.guild", "vertical.position"), function(trait) {
  dat <- observed_merged[!is.na(observed_merged[[trait]]) & !is.na(observed_merged$top20_ED_100tree), ]
  tab <- table(dat$top20_ED_100tree, dat[[trait]])
  ft <- fisher.test(tab)
  data.frame(trait = trait, p_value = ft$p.value, n_species = nrow(dat))
}) |> dplyr::bind_rows()
write.csv(observed_fisher_results, "outputs/tables/trait_observed_fisher.csv", row.names = FALSE)

message("Observed-only trait sensitivity analysis complete.")
