# 06_figures_tables.R
# Generate compact final figures from the formal summaries.

# Figure 1: ED Top 20 by mean ED
p_ed_top20 <- top20_ED_final |>
  dplyr::mutate(
    species_label = normalise_species(tip.label),
    species_label = reorder(species_label, mean_ED)
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(x = species_label, y = mean_ED)
  ) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = NULL,
    y = "Mean evolutionary distinctiveness",
    title = "Overall ED Top 20 across 100 phylogenetic trees"
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(
  "outputs/figures/FINAL_ED_top20.png",
  p_ed_top20,
  width = 8,
  height = 7,
  dpi = 300
)

# Figure 2: Top-20 frequency of the overall ED Top 20
p_stability <- top20_ED_final |>
  dplyr::mutate(
    species_label = normalise_species(tip.label),
    species_label = reorder(species_label, top20_frequency)
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(x = species_label, y = top20_frequency)
  ) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Proportion of trees in ED Top 20",
    title = "Phylogenetic stability of overall ED Top-20 candidates"
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(
  "outputs/figures/FINAL_ED_top20_stability.png",
  p_stability,
  width = 8,
  height = 7,
  dpi = 300
)

# Figure 3: Pairwise Top-20 overlap
plot_overlap <- overlap_summary[
  overlap_summary$comparison != "All three",
]

p_overlap <- ggplot2::ggplot(
  plot_overlap,
  ggplot2::aes(x = comparison, y = proportion)
) +
  ggplot2::geom_col() +
  ggplot2::geom_text(
    ggplot2::aes(
      label = paste0(shared_species, "/20")
    ),
    vjust = -0.4
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Proportion of shared Top-20 species",
    title = "Top-20 overlap among ED and time-bounded EvoHeritage"
  ) +
  ggplot2::theme_minimal(base_size = 11)

ggplot2::ggsave(
  "outputs/figures/FINAL_top20_overlap.png",
  p_overlap,
  width = 7,
  height = 5,
  dpi = 300
)

capture.output(
  sessionInfo(),
  file = "outputs/logs/R_session_info_final.txt"
)

message("Final figures and session information written.")
