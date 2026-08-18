# 05_final_summary.R
# Summarise formal EvoHeritage results and compare ED, 66 Ma and 145 Ma.

eh66_files_final <- sort(list.files("outputs/EvoHeritage_100trees_final", pattern = "^EH66_tree_[0-9]{3}\\.csv$", full.names = TRUE))
eh145_files_final <- sort(list.files("outputs/EvoHeritage_100trees_final", pattern = "^EH145_tree_[0-9]{3}\\.csv$", full.names = TRUE))
if (length(eh66_files_final) != 100 || length(eh145_files_final) != 100) stop("Expected 100 result files for each EvoHeritage boundary.")

all_eh66_final <- lapply(seq_along(eh66_files_final), function(i) { x <- read.csv(eh66_files_final[i], check.names=FALSE); x$tree_id <- i; x }) |> dplyr::bind_rows()
all_eh145_final <- lapply(seq_along(eh145_files_final), function(i) { x <- read.csv(eh145_files_final[i], check.names=FALSE); x$tree_id <- i; x }) |> dplyr::bind_rows()

eh66_summary_final <- all_eh66_final |> dplyr::group_by(tip.label) |> dplyr::summarise(
  mean_EH_66=mean(partitioned.EvoHeritage,na.rm=TRUE), median_EH_66=median(partitioned.EvoHeritage,na.rm=TRUE), sd_EH_66=sd(partitioned.EvoHeritage,na.rm=TRUE),
  mean_rank_66=mean(lf.rank,na.rm=TRUE), median_rank_66=median(lf.rank,na.rm=TRUE), top20_count_66=sum(lf.rank<=20,na.rm=TRUE), top20_frequency_66=top20_count_66/100, .groups="drop")
eh145_summary_final <- all_eh145_final |> dplyr::group_by(tip.label) |> dplyr::summarise(
  mean_EH_145=mean(partitioned.EvoHeritage,na.rm=TRUE), median_EH_145=median(partitioned.EvoHeritage,na.rm=TRUE), sd_EH_145=sd(partitioned.EvoHeritage,na.rm=TRUE),
  mean_rank_145=mean(lf.rank,na.rm=TRUE), median_rank_145=median(lf.rank,na.rm=TRUE), top20_count_145=sum(lf.rank<=20,na.rm=TRUE), top20_frequency_145=top20_count_145/100, .groups="drop")

comparison_final <- ed_summary_100 |>
  dplyr::select(tip.label,species,mean_ED,median_ED,sd_ED,mean_rank,median_rank,top20_count,top20_frequency) |>
  dplyr::left_join(eh66_summary_final,by="tip.label") |>
  dplyr::left_join(eh145_summary_final,by="tip.label")

top20_ED_final <- comparison_final |> dplyr::arrange(mean_rank) |> dplyr::slice(1:20)
top20_EH66_final <- comparison_final |> dplyr::arrange(mean_rank_66) |> dplyr::slice(1:20)
top20_EH145_final <- comparison_final |> dplyr::arrange(mean_rank_145) |> dplyr::slice(1:20)

ED_species_final <- top20_ED_final$tip.label; EH66_species_final <- top20_EH66_final$tip.label; EH145_species_final <- top20_EH145_final$tip.label
common_all_three_final <- Reduce(intersect,list(ED_species_final,EH66_species_final,EH145_species_final))
overlap_summary <- data.frame(
  comparison=c("ED vs 66 Ma","ED vs 145 Ma","66 Ma vs 145 Ma","All three"),
  shared_species=c(length(intersect(ED_species_final,EH66_species_final)),length(intersect(ED_species_final,EH145_species_final)),length(intersect(EH66_species_final,EH145_species_final)),length(common_all_three_final)),
  denominator=20
) |> dplyr::mutate(proportion=shared_species/denominator)

rank_correlations <- data.frame(
  comparison=c("ED vs 66 Ma","ED vs 145 Ma","66 Ma vs 145 Ma"),
  spearman_rho=c(
    cor(comparison_final$mean_rank,comparison_final$mean_rank_66,method="spearman",use="complete.obs"),
    cor(comparison_final$mean_rank,comparison_final$mean_rank_145,method="spearman",use="complete.obs"),
    cor(comparison_final$mean_rank_66,comparison_final$mean_rank_145,method="spearman",use="complete.obs")
  )
)

# Major-group mapping must contain exactly: species, major_group.
group_file <- "data_processed/chondrichthyan_major_groups.csv"
if (!file.exists(group_file)) stop("Missing ", group_file, ". Create it from the verified 1192-species mapping (500 shark, 600 ray, 49 chimaera, 43 unknown).")
major_group_mapping <- read.csv(group_file, stringsAsFactors=FALSE) |> dplyr::select(species,major_group)
major_group_mapping$species <- normalise_species(major_group_mapping$species)
if (nrow(major_group_mapping)!=1192 || anyDuplicated(major_group_mapping$species)) stop("Major-group mapping must contain 1192 unique species.")
comparison_with_groups <- comparison_final |> dplyr::left_join(major_group_mapping,by="species")
if (any(is.na(comparison_with_groups$major_group))) stop("Some phylogeny species lack major_group mapping.")

group_summary_final <- comparison_with_groups |>
  dplyr::filter(major_group %in% c("shark","ray","chimaera")) |>
  dplyr::group_by(major_group) |>
  dplyr::summarise(
    n_species=dplyr::n(),
    group_mean_ED=mean(mean_ED,na.rm=TRUE), group_median_ED=median(mean_ED,na.rm=TRUE),
    group_mean_EH66=mean(mean_EH_66,na.rm=TRUE), group_median_EH66=median(mean_EH_66,na.rm=TRUE),
    group_mean_EH145=mean(mean_EH_145,na.rm=TRUE), group_median_EH145=median(mean_EH_145,na.rm=TRUE), .groups="drop")

top20_group_counts_final <- dplyr::bind_rows(
  top20_ED_final |> dplyr::left_join(major_group_mapping,by="species") |> dplyr::count(major_group,name="n") |> dplyr::mutate(method="ED"),
  top20_EH66_final |> dplyr::left_join(major_group_mapping,by="species") |> dplyr::count(major_group,name="n") |> dplyr::mutate(method="66 Ma"),
  top20_EH145_final |> dplyr::left_join(major_group_mapping,by="species") |> dplyr::count(major_group,name="n") |> dplyr::mutate(method="145 Ma")
)

write.csv(comparison_final,"outputs/tables/ED_EvoHeritage_66_145_FINAL_100trees.csv",row.names=FALSE)
write.csv(top20_ED_final,"outputs/tables/FINAL_top20_ED.csv",row.names=FALSE)
write.csv(top20_EH66_final,"outputs/tables/FINAL_top20_EH66.csv",row.names=FALSE)
write.csv(top20_EH145_final,"outputs/tables/FINAL_top20_EH145.csv",row.names=FALSE)
write.csv(overlap_summary,"outputs/tables/FINAL_top20_overlap_summary.csv",row.names=FALSE)
write.csv(data.frame(tip.label=common_all_three_final),"outputs/tables/FINAL_species_shared_by_all_three.csv",row.names=FALSE)
write.csv(rank_correlations,"outputs/tables/FINAL_rank_correlations.csv",row.names=FALSE)
write.csv(group_summary_final,"outputs/tables/FINAL_major_group_summary.csv",row.names=FALSE)
write.csv(top20_group_counts_final,"outputs/tables/FINAL_top20_group_counts.csv",row.names=FALSE)

print(overlap_summary); print(rank_correlations); print(group_summary_final); print(top20_group_counts_final)
message("Final ED/EvoHeritage, correlation and major-group summaries complete.")
