#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else getwd()
data_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(project_root, "AMRC-repo-files", "AMR-cartography-results", "data")
}
summary_csv <- if (length(args) >= 3) {
  args[[3]]
} else {
  file.path(
    project_root,
    "AMRC-repo-files",
    "Streptococcus pneumoniae analysis",
    "Single_subs_all_S.pneumo.csv"
  )
}
output_dir <- if (length(args) >= 4) {
  args[[4]]
} else {
  file.path(project_root, "manuscript", "source_data")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(project_root, "analysis", "config.R"))
load(file.path(data_dir, "tablemic_pneumo_gen_3628.RData"))
load(file.path(data_dir, "Spneumo_3628_PCA_start_2D_METRIC.RData"))
metadata <- read.csv(
  file.path(data_dir, "meta_data_Spneumoniae.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

excluded_isolates <- c(
  "20156696", "20162849", "20151885", "20153985", "20154509",
  "2013224047", "2013218247", "2014200662", "5869-99", "2513-99"
)

theta <- AMRC_MAP_ROTATION_DEGREES * pi / 180
rotation <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), ncol = 2)
rotated_coordinates <- torg_met$conf %*% rotation

map_coordinates <- bind_cols(as.data.frame(rotated_coordinates), metadata) %>%
  filter(!LABID %in% excluded_isolates) %>%
  rename(D1 = V1, D2 = V2) %>%
  mutate(
    D1 = D1 / AMRC_PHEN_DILATION_SLOPE,
    D2 = D2 / AMRC_PHEN_DILATION_SLOPE
  )

if (nrow(map_coordinates) != nrow(PBPseq)) {
  stop(
    sprintf(
      "Filtered map rows (%d) do not match PBP sequence rows (%d)",
      nrow(map_coordinates), nrow(PBPseq)
    ),
    call. = FALSE
  )
}

# Reproduce the PBP-type numbering used by Script 17. PBP2B_N656 is excluded
# because its apparent variation is an amino-acid calling error.
pbp_type <- PBPseq %>%
  select(PBP1A_T371:PBP2X_V587, -PBP2B_N656) %>%
  unite("PBP_type_sequence", everything(), remove = TRUE) %>%
  group_by(PBP_type_sequence) %>%
  group_indices()

pbp_data <- PBPseq %>% mutate(PBP_type = as.character(pbp_type))
map_data <- bind_cols(map_coordinates, pbp_data)

single_subs <- read.csv(summary_csv, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c(
  "Relative_Comparison", "Reference_PBP", "Comparison_PBP",
  "number_of_isolates_of_reference", "number_of_isolates_of_comparison_group",
  "Location", "PBP", "Loci", "AA_1", "AA_2",
  "median_phenotypic_distance", "sd", "n", "se"
)
missing_columns <- setdiff(required_columns, colnames(single_subs))
if (length(missing_columns) > 0) {
  stop(
    sprintf("Missing required summary columns: %s", paste(missing_columns, collapse = ", ")),
    call. = FALSE
  )
}

# The legacy `median_phenotypic_distance` field is produced with mean() in
# Script 17. Rename it in the figure data and select the same >=1-unit panels.
comparison_summary <- single_subs %>%
  mutate(
    pbp_order = factor(PBP, levels = c("PBP1A", "PBP2B", "PBP2X")),
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
  filter(mean_pairwise_distance_log2_mic >= 1) %>%
  arrange(pbp_order, position, amino_acid_contrast, genetic_background_index) %>%
  mutate(
    display_contrast = if_else(
      genetic_background_index == 1L,
      amino_acid_contrast,
      paste0(amino_acid_contrast, " (", genetic_background_index, ")")
    ),
    panel_id = sprintf("panel_%02d", row_number()),
    panel_label = sprintf(
      "%s %s - Mean: %.3f",
      sub("^PBP", "", PBP),
      display_contrast,
      mean_pairwise_distance_log2_mic
    )
  )

if (nrow(comparison_summary) == 0) {
  stop("No comparisons met the >=1 log2 MIC-unit display threshold", call. = FALSE)
}

build_panel_data <- function(i) {
  comparison <- comparison_summary[i, , drop = FALSE]
  location <- comparison$Location[[1]]
  reference_type <- as.character(comparison$Reference_PBP[[1]])
  comparison_type <- as.character(comparison$Comparison_PBP[[1]])

  if (!location %in% colnames(map_data)) {
    stop(sprintf("PBP position column not found: %s", location), call. = FALSE)
  }

  points <- map_data %>%
    filter(PBP_type %in% c(reference_type, comparison_type)) %>%
    transmute(
      panel_id = comparison$panel_id[[1]],
      panel_label = comparison$panel_label[[1]],
      relative_comparison = comparison$Relative_Comparison[[1]],
      pbp = comparison$PBP[[1]],
      position = comparison$position[[1]],
      amino_acid_contrast = comparison$amino_acid_contrast[[1]],
      mean_pairwise_distance_log2_mic = comparison$mean_pairwise_distance_log2_mic[[1]],
      isolate_id = LABID,
      pbp_type = PBP_type,
      comparison_group = if_else(PBP_type == reference_type, "reference", "comparison"),
      amino_acid = .data[[location]],
      D1,
      D2
    )

  observed_amino_acids <- sort(unique(points$amino_acid))
  expected_amino_acids <- sort(unique(c(comparison$AA_1[[1]], comparison$AA_2[[1]])))
  if (!identical(observed_amino_acids, expected_amino_acids)) {
    stop(
      sprintf(
        "%s: observed amino acids [%s], expected [%s]",
        comparison$Relative_Comparison[[1]],
        paste(observed_amino_acids, collapse = ","),
        paste(expected_amino_acids, collapse = ",")
      ),
      call. = FALSE
    )
  }

  reference_points <- points %>%
    filter(comparison_group == "reference") %>%
    select(reference_isolate_id = isolate_id, D1_reference = D1, D2_reference = D2)
  comparison_points <- points %>%
    filter(comparison_group == "comparison") %>%
    select(comparison_isolate_id = isolate_id, D1_comparison = D1, D2_comparison = D2)

  segments <- crossing(reference_points, comparison_points) %>%
    mutate(
      panel_id = comparison$panel_id[[1]],
      panel_label = comparison$panel_label[[1]],
      relative_comparison = comparison$Relative_Comparison[[1]],
      .before = 1
    )

  list(points = points, segments = segments)
}

panel_data <- lapply(seq_len(nrow(comparison_summary)), build_panel_data)
highlight_points <- bind_rows(lapply(panel_data, `[[`, "points"))
pairwise_segments <- bind_rows(lapply(panel_data, `[[`, "segments"))

panel_levels <- comparison_summary$panel_label
comparison_summary$panel_label <- factor(comparison_summary$panel_label, levels = panel_levels)
highlight_points$panel_label <- factor(highlight_points$panel_label, levels = panel_levels)
pairwise_segments$panel_label <- factor(pairwise_segments$panel_label, levels = panel_levels)

background_points <- crossing(
  comparison_summary %>%
    transmute(
      panel_id,
      panel_label,
      reference_type = as.character(Reference_PBP),
      comparison_type = as.character(Comparison_PBP)
    ),
  map_data %>% select(isolate_id = LABID, pbp_type = PBP_type, D1, D2)
) %>%
  filter(!pbp_type %in% c(reference_type, comparison_type))

summary_source <- comparison_summary %>%
  transmute(
    relative_comparison = Relative_Comparison,
    pbp = PBP,
    position,
    amino_acid_contrast,
    genetic_background_index,
    display_contrast,
    reference_pbp_type = Reference_PBP,
    comparison_pbp_type = Comparison_PBP,
    number_of_isolates_in_reference_group = number_of_isolates_of_reference,
    number_of_isolates_in_comparison_group = number_of_isolates_of_comparison_group,
    mean_pairwise_distance_log2_mic,
    sd_pairwise_distance,
    n_pairwise_isolate_comparisons,
    se_mean_pairwise_distance,
    panel_label
  )

write.csv(
  summary_source,
  file.path(output_dir, "FigureS15_comparison_summary.csv"),
  row.names = FALSE
)
write.csv(
  highlight_points,
  file.path(output_dir, "FigureS15_isolate_points.csv"),
  row.names = FALSE
)
write.csv(
  pairwise_segments,
  file.path(output_dir, "FigureS15_pairwise_segments.csv"),
  row.names = FALSE
)

amino_acid_colours <- c(
  "A" = "#E41A1C", "D" = "#377EB8", "E" = "#4DAF4A",
  "G" = "#984EA3", "I" = "#FFFF33", "K" = "#FF7F00",
  "L" = "#66A61E", "M" = "#A65628", "N" = "#F781BF",
  "P" = "#A6761D", "S" = "#E7298A", "T" = "#1B9E77",
  "V" = "#7570B3"
)
amino_acids_in_figure <- sort(unique(highlight_points$amino_acid))
missing_colours <- setdiff(amino_acids_in_figure, names(amino_acid_colours))
if (length(missing_colours) > 0) {
  stop(
    sprintf("No plotting colour defined for amino acids: %s", paste(missing_colours, collapse = ", ")),
    call. = FALSE
  )
}

x_limits <- range(map_data$D1, na.rm = TRUE)
y_limits <- range(map_data$D2, na.rm = TRUE)

figure_s15 <- ggplot() +
  geom_point(
    data = background_points,
    aes(x = D1, y = D2),
    shape = 16,
    size = 0.28,
    colour = "grey62",
    alpha = 0.75
  ) +
  geom_segment(
    data = pairwise_segments,
    aes(
      x = D1_reference,
      y = D2_reference,
      xend = D1_comparison,
      yend = D2_comparison
    ),
    linewidth = 0.18,
    colour = "grey15",
    alpha = 0.65
  ) +
  geom_point(
    data = highlight_points,
    aes(x = D1, y = D2, fill = amino_acid),
    shape = 21,
    size = 1.15,
    stroke = 0.28,
    colour = "black"
  ) +
  facet_wrap(~panel_label, ncol = 5) +
  scale_fill_manual(
    values = amino_acid_colours,
    breaks = amino_acids_in_figure,
    name = "Amino acid"
  ) +
  scale_x_continuous(limits = x_limits, breaks = seq(floor(x_limits[[1]]), ceiling(x_limits[[2]]), 1)) +
  scale_y_continuous(limits = y_limits, breaks = seq(floor(y_limits[[1]]), ceiling(y_limits[[2]]), 1)) +
  coord_fixed() +
  theme_linedraw(base_size = 9) +
  theme(
    panel.grid.major = element_line(colour = "grey72", linewidth = 0.28),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(size = 7.2, colour = "black", face = "plain"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.45, "cm"),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.25),
    panel.spacing = unit(0.12, "cm")
  ) +
  guides(fill = guide_legend(override.aes = list(size = 2.4)))

png_path <- file.path(output_dir, "FigureS15_single_substitution_maps.png")
pdf_path <- file.path(output_dir, "FigureS15_single_substitution_maps.pdf")
ggsave(png_path, plot = figure_s15, width = 10.5, height = 6.3, dpi = 600, bg = "white")
ggsave(pdf_path, plot = figure_s15, width = 10.5, height = 6.3, bg = "white")

message(sprintf("Displayed comparisons (>=1): %d", nrow(comparison_summary)))
message(sprintf("Highlighted isolate rows: %d", nrow(highlight_points)))
message(sprintf("Pairwise segments: %d", nrow(pairwise_segments)))
message(sprintf("Wrote %s", png_path))
message(sprintf("Wrote %s", pdf_path))
