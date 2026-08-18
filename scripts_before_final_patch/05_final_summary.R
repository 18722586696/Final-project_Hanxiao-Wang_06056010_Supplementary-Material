# 05_final_summary.R
# Summarise formal EvoHeritage results and compare ED, 66 Ma and 145 Ma.

eh66_files_final <- sort(list.files(
  "outputs/EvoHeritage_100trees_final",
  pattern = "^EH66_tree_[0-9]{3}\\.csv$",
  full.names = TRUE
))

all_eh66_final <- lapply(seq_along(eh66_files_final), function(i) {
  x <- read.csv(eh66_files_final[i], check.names = FALSE)
  x$tree_id <- i
  x
}) |>
  dplyr::bind_rows()

eh66_summary_final <- all_eh66_final |>
  dplyr::group_by(tip.label) |>
  dplyr::summarise(
    mean_EH_66 = mean(partitioned.EvoHeritage, na.rm = TRUE),
    median_EH_66 = median(partitioned.EvoHeritage, na.rm = TRUE),
    sd_EH_66 = sd(partitioned.EvoHeritage, na.rm = TRUE),
    mean_rank_66 = mean(lf.rank, na.rm = TRUE),
    median_rank_66 = median(lf.rank, na.rm = TRUE),
    top20_count_66 = sum(lf.rank <= 20, na.rm = TRUE),
    top20_frequency_66 = top20_count_66 / 100,
    .groups = "drop"
  ) |>
  dplyr::arrange(mean_rank_66)

eh145_files_final <- sort(list.files(
  "outputs/EvoHeritage_100trees_final",
  pattern = "^EH145_tree_[0-9]{3}\\.csv$",
  full.names = TRUE
))

all_eh145_final <- lapply(seq_along(eh145_files_final), function(i) {
  x <- read.csv(eh145_files_final[i], check.names = FALSE)
  x$tree_id <- i
  x
}) |>
  dplyr::bind_rows()

eh145_summary_final <- all_eh145_final |>
  dplyr::group_by(tip.label) |>
  dplyr::summarise(
    mean_EH_145 = mean(partitioned.EvoHeritage, na.rm = TRUE),
    median_EH_145 = median(partitioned.EvoHeritage, na.rm = TRUE),
    sd_EH_145 = sd(partitioned.EvoHeritage, na.rm = TRUE),
    mean_rank_145 = mean(lf.rank, na.rm = TRUE),
    median_rank_145 = median(lf.rank, na.rm = TRUE),
    top20_count_145 = sum(lf.rank <= 20, na.rm = TRUE),
    top20_frequency_145 = top20_count_145 / 100,
    .groups = "drop"
  ) |>
  dplyr::arrange(mean_rank_145)

comparison_final <- ed_summary_100 |>
  dplyr::select(
    tip.label,
    species,
    mean_ED,
    median_ED,
    sd_ED,
    mean_rank,
    median_rank,
    top20_count,
    top20_frequency
  ) |>
  dplyr::left_join(
    eh66_summary_final,
    by = "tip.label"
  ) |>
  dplyr::left_join(
    eh145_summary_final,
    by = "tip.label"
  )

top20_ED_final <- comparison_final |>
  dplyr::arrange(mean_rank) |>
  dplyr::slice(1:20)

top20_EH66_final <- comparison_final |>
  dplyr::arrange(mean_rank_66) |>
  dplyr::slice(1:20)

top20_EH145_final <- comparison_final |>
  dplyr::arrange(mean_rank_145) |>
  dplyr::slice(1:20)

ED_species_final <- top20_ED_final$tip.label
EH66_species_final <- top20_EH66_final$tip.label
EH145_species_final <- top20_EH145_final$tip.label

overlap_ED_66 <- length(intersect(
  ED_species_final,
  EH66_species_final
))

overlap_ED_145 <- length(intersect(
  ED_species_final,
  EH145_species_final
))

overlap_66_145 <- length(intersect(
  EH66_species_final,
  EH145_species_final
))

common_all_three_final <- Reduce(
  intersect,
  list(
    ED_species_final,
    EH66_species_final,
    EH145_species_final
  )
)

overlap_summary <- data.frame(
  comparison = c(
    "ED vs 66 Ma",
    "ED vs 145 Ma",
    "66 Ma vs 145 Ma",
    "All three"
  ),
  shared_species = c(
    overlap_ED_66,
    overlap_ED_145,
    overlap_66_145,
    length(common_all_three_final)
  ),
  denominator = c(20, 20, 20, 20),
  proportion = c(
    overlap_ED_66 / 20,
    overlap_ED_145 / 20,
    overlap_66_145 / 20,
    length(common_all_three_final) / 20
  )
)

write.csv(
  comparison_final,
  "outputs/tables/ED_EvoHeritage_66_145_FINAL_100trees.csv",
  row.names = FALSE
)

write.csv(
  top20_ED_final,
  "outputs/tables/FINAL_top20_ED.csv",
  row.names = FALSE
)

write.csv(
  top20_EH66_final,
  "outputs/tables/FINAL_top20_EH66.csv",
  row.names = FALSE
)

write.csv(
  top20_EH145_final,
  "outputs/tables/FINAL_top20_EH145.csv",
  row.names = FALSE
)

write.csv(
  overlap_summary,
  "outputs/tables/FINAL_top20_overlap_summary.csv",
  row.names = FALSE
)

write.csv(
  data.frame(tip.label = common_all_three_final),
  "outputs/tables/FINAL_species_shared_by_all_three.csv",
  row.names = FALSE
)

print(overlap_summary)
print(common_all_three_final)

message("Final ED/EvoHeritage comparison complete.")
