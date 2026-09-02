#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else getwd()
input_csv <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(
    project_root,
    "AMRC-repo-files",
    "Streptococcus pneumoniae analysis",
    "Single_subs_all_S.pneumo.csv"
  )
}
output_dir <- if (length(args) >= 3) {
  args[[3]]
} else {
  file.path(project_root, "manuscript", "source_data")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

single_subs <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c(
  "Relative_Comparison", "PBP", "Loci", "AA_1", "AA_2",
  "median_phenotypic_distance", "sd", "n", "se"
)
missing_columns <- setdiff(required_columns, colnames(single_subs))
if (length(missing_columns) > 0) {
  stop(
    sprintf("Missing required input columns: %s", paste(missing_columns, collapse = ", ")),
    call. = FALSE
  )
}

# `median_phenotypic_distance` is a legacy column name. Script 17 calculates it
# with mean(Phenotypic_distance), and `se` is therefore the standard error of
# that mean. Give both quantities their correct names in the figure source data.
figure_data <- single_subs %>%
  mutate(
    position = as.integer(Loci),
    amino_acid_contrast = paste0(AA_1, position, AA_2),
    mean_pairwise_distance_log2_mic = as.numeric(median_phenotypic_distance),
    sd_pairwise_distance = as.numeric(sd),
    n_pairwise_isolate_comparisons = as.integer(n),
    se_mean_pairwise_distance = as.numeric(se)
  ) %>%
  group_by(PBP, amino_acid_contrast) %>%
  mutate(genetic_background_index = row_number()) %>%
  ungroup() %>%
  mutate(
    display_label = if_else(
      genetic_background_index == 1L,
      amino_acid_contrast,
      paste0(amino_acid_contrast, " (", genetic_background_index, ")")
    )
  ) %>%
  filter(mean_pairwise_distance_log2_mic > 0.5) %>%
  arrange(PBP, position, amino_acid_contrast, genetic_background_index) %>%
  mutate(
    plot_key = paste(PBP, display_label, Relative_Comparison, sep = "::"),
    plot_key = factor(plot_key, levels = unique(plot_key))
  )

source_data <- figure_data %>%
  transmute(
    relative_comparison = Relative_Comparison,
    pbp = PBP,
    position,
    amino_acid_contrast,
    genetic_background_index,
    display_label,
    mean_pairwise_distance_log2_mic,
    sd_pairwise_distance,
    n_pairwise_isolate_comparisons,
    se_mean_pairwise_distance
  )

source_csv <- file.path(output_dir, "FigureS14_single_substitution_comparisons.csv")
write.csv(source_data, source_csv, row.names = FALSE)

red <- "#E41A1C"
figure_s14 <- ggplot(
  figure_data,
  aes(x = mean_pairwise_distance_log2_mic, y = plot_key)
) +
  geom_vline(xintercept = 1, linewidth = 0.65, colour = "black") +
  geom_segment(
    data = filter(figure_data, !is.na(se_mean_pairwise_distance)),
    aes(
      x = pmax(0, mean_pairwise_distance_log2_mic - se_mean_pairwise_distance),
      xend = mean_pairwise_distance_log2_mic + se_mean_pairwise_distance,
      yend = plot_key
    ),
    linewidth = 0.65,
    colour = "black"
  ) +
  geom_point(shape = 21, size = 3, stroke = 0.35, colour = red, fill = red, alpha = 0.70) +
  facet_wrap(~PBP, scales = "free_y", nrow = 1) +
  scale_x_continuous(limits = c(0, 5), breaks = 0:5, expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_discrete(labels = function(x) sub("^[^:]+::([^:]+)::.*$", "\\1", x)) +
  labs(
    x = expression("Mean pairwise phenotypic distance (log"[2] * " MIC units)"),
    y = "Amino-acid contrast (PBP TPD region)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_line(colour = "grey80", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey85", colour = "grey40"),
    strip.text = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.position = "none"
  )

png_path <- file.path(output_dir, "FigureS14_single_substitution_comparisons.png")
pdf_path <- file.path(output_dir, "FigureS14_single_substitution_comparisons.pdf")
ggsave(png_path, plot = figure_s14, width = 8.2, height = 5.1, dpi = 600, bg = "white")
ggsave(pdf_path, plot = figure_s14, width = 8.2, height = 5.1, bg = "white")

message(sprintf("Input comparisons: %d", nrow(single_subs)))
message(sprintf("Displayed comparisons (>0.5): %d", nrow(figure_data)))
message(sprintf("Displayed by PBP: %s", paste(names(table(figure_data$PBP)), table(figure_data$PBP), collapse = ", ")))
message(sprintf("Wrote %s", source_csv))
message(sprintf("Wrote %s", png_path))
message(sprintf("Wrote %s", pdf_path))
