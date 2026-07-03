#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(jsonlite)
  library(patchwork)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) args[[1]] else getwd()
results_root <- if (length(args) >= 2) args[[2]] else file.path(project_root, "analysis_outputs", "recomputed_170_thresholds")
out_dir <- if (length(args) >= 3) args[[3]] else file.path(results_root, "manuscript_outputs", "supplement_figures", "original_style")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

thresholds <- fromJSON(file.path(results_root, "thresholds", "recomputed_thresholds.json"))
additive_threshold <- thresholds$additive$adjusted_threshold
epistasis_threshold <- thresholds$epistasis$epistasis_threshold
additive_meff <- thresholds$additive$galwey_meff

red <- "#E41A1C"
pbp_colours <- c("1A" = "#FDBB84", "2B" = "#EF6548", "2X" = "#B30000")

require_exact_package <- function(package, figure) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      sprintf(
        "%s requires the original plotting package `%s`; refusing to draw a non-identical substitute.",
        figure, package
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

parse_pbp <- function(marker) {
  str_match(marker, "^PBP(1A|2B|2X)_")[, 2]
}

first_position <- function(marker) {
  as.integer(str_match(marker, "^[^0-9]*([0-9]+)")[, 2])
}

p_to_numeric <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, "^<", "")
  suppressWarnings(as.numeric(x))
}

read_csv_base <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

save_plot <- function(plot, stem, width, height) {
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  ggsave(png_path, plot = plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(pdf_path, plot = plot, width = width, height = height, bg = "white")
  invisible(c(png = png_path, pdf = pdf_path))
}

save_current_device_plot <- function(plot_fun, stem, width, height) {
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  png(png_path, width = width, height = height, units = "in", res = 300, bg = "white")
  plot_fun()
  dev.off()
  pdf(pdf_path, width = width, height = height, bg = "white")
  plot_fun()
  dev.off()
  invisible(c(png = png_path, pdf = pdf_path))
}

convert_svg_outputs <- function(svg_path, stem) {
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))

  if (nzchar(Sys.which("inkscape"))) {
    pdf_status <- system2(
      "inkscape",
      c(svg_path, "--export-type=pdf", paste0("--export-filename=", pdf_path))
    )
    png_status <- system2(
      "inkscape",
      c(svg_path, "--export-type=png", "--export-dpi=300", paste0("--export-filename=", png_path))
    )
    if (!identical(pdf_status, 0L) || !identical(png_status, 0L)) {
      stop("Inkscape failed while converting the nVennR SVG output.", call. = FALSE)
    }
  } else if (nzchar(Sys.which("convert"))) {
    pdf_status <- system2("convert", c(svg_path, pdf_path))
    png_status <- system2("convert", c("-density", "300", svg_path, png_path))
    if (!identical(pdf_status, 0L) || !identical(png_status, 0L)) {
      stop("ImageMagick failed while converting the nVennR SVG output.", call. = FALSE)
    }
  } else {
    stop("Need `inkscape` or ImageMagick `convert` to convert nVennR SVG output.", call. = FALSE)
  }

  invisible(c(svg = svg_path, png = png_path, pdf = pdf_path))
}

pad_nvennr_svg_canvas <- function(svg_path, width = 1100, height = 500) {
  svg <- readLines(svg_path, warn = FALSE)
  svg[1] <- sub('width="[0-9.]+"', sprintf('width="%s"', width), svg[1])
  svg[1] <- sub('height="[0-9.]+"', sprintf('height="%s"', height), svg[1])
  if (!grepl("viewBox=", svg[1], fixed = TRUE)) {
    svg[1] <- sub(">$", sprintf(' viewBox="0 0 %s %s">', width, height), svg[1])
  }
  if (!any(grepl('id="nvennr-white-background"', svg, fixed = TRUE))) {
    svg <- append(svg, '<rect id="nvennr-white-background" width="100%" height="100%" fill="white"/>', after = 1)
  }
  writeLines(svg, svg_path)
  invisible(svg_path)
}

split_compound_marker <- function(marker) {
  parts <- unlist(str_split(marker, "_PBP", simplify = FALSE), use.names = FALSE)
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) {
    return(character(0))
  }
  parts <- ifelse(str_starts(parts, "PBP"), parts, paste0("PBP", parts))
  parts
}

marker_to_position_label <- function(marker) {
  parts <- split_compound_marker(marker)
  parts <- str_replace(parts, "_[A-Z*?]$", "")
  str_remove(parts, "^PBP")
}

position_number_from_label <- function(label) {
  as.integer(str_match(label, "_[A-Z]?([0-9]+)")[, 2])
}

load_tpd_position_labels <- function() {
  pbp_sequence_path <- file.path(project_root, "AMRC-repo-files", "AMR-cartography-results", "data", "PBP_Sequence_dataset2.csv")
  pbp_columns <- colnames(read_csv_base(pbp_sequence_path))
  pbp_columns <- pbp_columns[str_starts(pbp_columns, "PBP")]
  str_remove(pbp_columns, "^PBP")
}

load_additive_effects <- function() {
  pvals <- read_csv_base(file.path(results_root, "additive", "merged", "mvLMM_p_values_normal_pneumo_low_freq_vars.csv"))
  effects <- read_csv_base(file.path(results_root, "additive", "merged", "mvLMM_effect_sizes_normal_pneumo_low_freq_vars.csv"))

  candidates <- effects %>%
    filter(effect_type == "candidate") %>%
    select(marker, env, effsize, effsize_se) %>%
    pivot_wider(names_from = env, values_from = c(effsize, effsize_se))

  pvals %>%
    mutate(pv20_adj_galwey = pmin(pv20 * additive_meff, 1)) %>%
    select(marker, pv20, pv20_adj_galwey) %>%
    left_join(candidates, by = "marker") %>%
    mutate(
      effsize_D1 = effsize_env1_D1,
      effsize_D2 = effsize_env1_D2,
      effsize_se_D1 = effsize_se_env1_D1,
      effsize_se_D2 = effsize_se_env1_D2,
      Joint_effsize = sqrt(effsize_D1^2 + effsize_D2^2),
      LMM = "mvLMM",
      PBP = parse_pbp(marker)
    )
}

load_epistasis_effects <- function() {
  epi <- read_csv_base(file.path(results_root, "epistasis", "merged", "corrected_epistasis_p_values.csv"))
  interaction_counts <- read_csv_base(file.path(results_root, "interactions", "corrected_epistasis_interactions.csv")) %>%
    select(interaction, min_cell_count)

  epi %>%
    left_join(interaction_counts, by = "interaction") %>%
    mutate(
      effsize_D1 = effsize_env1_D1,
      effsize_D2 = effsize_env1_D2,
      Joint_effsize = joint_effect_size,
      LMM = "Epistatic",
      neg_log10_adj_p = -log10(pmax(pv20_adj_galwey, .Machine$double.xmin))
    )
}

make_s17 <- function() {
  require_exact_package("ggraph", "Supplementary Figure S17B")
  require_exact_package("tidygraph", "Supplementary Figure S17B")

  epi <- load_epistasis_effects()
  marker_support <- read_csv_base(file.path(results_root, "epistasis", "merged", "corrected_epistasis_marker_support.csv")) %>%
    mutate(PBP = parse_pbp(marker), position = first_position(marker))

  s17_scatter_source <- epi %>%
    select(interaction, parent_a, parent_b, pv20_adj_galwey, Joint_effsize, effsize_D1, effsize_D2,
           joint_effect_size_se, min_cell_count, passes_epistasis_threshold,
           passes_epistasis_effect_filter, epistasis_support)
  write.csv(s17_scatter_source, file.path(out_dir, "recomputed_S17_epistasis_scatter_source.csv"), row.names = FALSE)
  write.csv(marker_support, file.path(out_dir, "recomputed_S17_marker_support_source.csv"), row.names = FALSE)

  p_scatter <- ggplot(
    filter(epi),
    aes(x = joint_effect_size, y = -log10(pv20_adj_galwey), size = min_cell_count)
  ) +
    geom_vline(xintercept = 1) +
    geom_hline(yintercept = -log10(epistasis_threshold)) +
    geom_point(shape = 16, colour = red, alpha = 0.25, position = position_jitter(width = 0., height = 0.)) +
    theme_linedraw() +
    scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1)) +
    labs(title = "", x = "Effect size (MIC units)", y = "-log10 adj. p-value") +
    guides(size = guide_legend(title = "No. isolates (min)")) +
    theme_bw() +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16)) +
    annotate("label", x = 4.5, y = -log10(epistasis_threshold), label = signif(epistasis_threshold, 2)) +
    scale_colour_manual(values = c(red, "#377EB8"))

  save_plot(p_scatter, "recomputed_S17A_spneumo_effect_size_plot", width = 7.0, height = 5.0)

  all_position_labels <- load_tpd_position_labels()
  node_tbl <- data.frame(label = all_position_labels, stringsAsFactors = FALSE) %>%
    mutate(
      PBP = str_match(label, "^(1A|2B|2X)_")[, 2],
      position = position_number_from_label(label)
    ) %>%
    arrange(PBP, position, label) %>%
    mutate(id = row_number(), .before = 1)
  write.csv(node_tbl, file.path(out_dir, "recomputed_S17_network_node_source.csv"), row.names = FALSE)

  interaction_network <- epi %>%
    filter(pv20_adj_galwey <= epistasis_threshold) %>%
    rowwise() %>%
    mutate(
      AA1_list = list(marker_to_position_label(parent_a)),
      AA2_list = list(marker_to_position_label(parent_b))
    ) %>%
    ungroup() %>%
    select(AA1_list, AA2_list, pv20_adj_galwey) %>%
    tidyr::unnest_longer(AA1_list, values_to = "AA1") %>%
    tidyr::unnest_longer(AA2_list, values_to = "AA2") %>%
    mutate(Significance = -log10(pv20_adj_galwey))

  edges <- interaction_network %>%
    left_join(node_tbl, by = c("AA1" = "label")) %>%
    rename(from = id) %>%
    left_join(node_tbl, by = c("AA2" = "label")) %>%
    rename(to = id) %>%
    filter(!is.na(from), !is.na(to)) %>%
    select(from, to, Significance)

  routes_tidy <- tidygraph::tbl_graph(nodes = node_tbl, edges = edges, directed = TRUE) %>%
    tidygraph::activate(edges) %>%
    arrange(desc(Significance))

  p_network <- ggraph::ggraph(routes_tidy, layout = "linear", circular = TRUE) +
    ggraph::geom_edge_arc(aes(width = Significance), width = 0.05, alpha = .85, colour = "#B30000") +
    ggraph::scale_edge_width(range = c(.05, .1)) +
    ggraph::geom_node_point(aes(colour = PBP), shape = 16) +
    scale_size(range = c(1, 4), guide = "none") +
    labs(edge_width = "-log10 adj. p-value") +
    ggraph::geom_node_label(aes(label = ifelse(!is.na(position) & position %% 50 == 0, position, NA)), size = 5, hjust = 0.5, vjust = 0.5, na.rm = TRUE) +
    ggraph::theme_graph() +
    coord_fixed(ratio = 1, clip = "off") +
    scale_colour_manual(values = pbp_colours) +
    theme(legend.position = "none", plot.margin = margin(12, 24, 12, 24)) +
    annotate("label", x = .9, y = .9, size = 6, label = "PBP1A") +
    annotate("label", x = 0, y = -1.15, size = 6, label = "PBP2B") +
    annotate("label", x = -.9, y = .9, size = 6, label = "PBP2X")

  save_plot(p_network, "recomputed_S17B_spneumo_epistatic_interaction_network", width = 7.0, height = 7.0)

  a_mean <- marker_support %>%
    summarize(median_PBP_count = median(num_sig_interactions), mean_PBP_count = mean(num_sig_interactions))

  p_hist <- ggplot(filter(marker_support), aes(x = num_sig_interactions)) +
    labs(title = "C", x = "Total number of significant interactions", y = "Frequency") +
    geom_histogram(fill = red, alpha = 0.65, colour = "black", bins = 50) +
    geom_vline(data = a_mean, aes(xintercept = median_PBP_count), color = red, linetype = "dashed", linewidth = 1) +
    scale_x_continuous(breaks = seq(0, 70, 10)) +
    scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 5)) +
    theme_bw() +
    coord_fixed(ratio = 0.4) +
    geom_label(data = a_mean, inherit.aes = FALSE, aes(x = 30, y = 12.5, label = paste("Median =", median_PBP_count)))

  save_plot(p_hist, "recomputed_S17C_spneumo_epistatic_histogram", width = 7.0, height = 5.0)

  combined <- ((p_scatter | p_network) + plot_layout(widths = c(1, 1.25))) / p_hist +
    plot_layout(heights = c(3.2, 1.0))
  save_plot(combined, "recomputed_S17_epistasis_original_style", width = 10.5, height = 6.9)

  list(
    tested_interactions = nrow(epi),
    p_threshold_only_interactions = sum(epi$pv20_adj_galwey <= epistasis_threshold, na.rm = TRUE),
    supported_interactions = sum(as.character(epi$epistasis_support) == "True", na.rm = TRUE),
    markers_with_support = nrow(marker_support),
    markers_gt_40_supported_interactions = sum(marker_support$num_sig_interactions > 40),
    epistasis_threshold = epistasis_threshold
  )
}

make_s18 <- function() {
  additive <- load_additive_effects() %>%
    filter(pv20_adj_galwey <= additive_threshold) %>%
    select(marker, PBP, pv20_adj_galwey, effsize_D1, effsize_D2, Joint_effsize, LMM)

  epistasis <- load_epistasis_effects() %>%
    filter(pv20_adj_galwey <= epistasis_threshold) %>%
    transmute(marker = interaction, PBP = NA_character_, pv20_adj_galwey, effsize_D1, effsize_D2, Joint_effsize, LMM)

  write.csv(additive, file.path(out_dir, "recomputed_S18_additive_source.csv"), row.names = FALSE)
  write.csv(epistasis, file.path(out_dir, "recomputed_S18_epistasis_source.csv"), row.names = FALSE)

  effect_sizes <- bind_rows(additive, epistasis)
  effect_sizes <- effect_sizes %>%
    mutate(
      category_V1 = case_when(
        (effsize_D1 < -1) ~ "-",
        (effsize_D1 < 1 & effsize_D1 > -1) ~ "/",
        (effsize_D1 > 1) ~ "+"
      ),
      category_V2 = case_when(
        (effsize_D2 < -1) ~ "-",
        (effsize_D2 < 1 & effsize_D2 > -1) ~ "/",
        (effsize_D2 > 1) ~ "+"
      )
    )

  test_n <- as.data.frame(t(as.data.frame.matrix(table(effect_sizes$category_V1, effect_sizes$category_V2))))
  test_prop <- as.data.frame(t(as.data.frame.matrix(prop.table(table(effect_sizes$category_V1, effect_sizes$category_V2)) * 100)))
  test_prop <- round(test_prop, 2)

  summary_rows <- expand.grid(category_V2 = c("-", "/", "+"), category_V1 = c("-", "/", "+"), stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(count = as.integer(test_n[category_V2, category_V1]), percent = as.numeric(test_prop[category_V2, category_V1])) %>%
    ungroup()
  write.csv(summary_rows, file.path(out_dir, "recomputed_S18_axis_category_summary.csv"), row.names = FALSE)

  label_text <- function(v2, v1) paste(test_n[v2, v1], " - ", test_prop[v2, v1], "%")

  p <- ggplot(filter(effect_sizes, LMM == "mvLMM"), aes(x = effsize_D1, y = effsize_D2)) +
    geom_vline(xintercept = 1, size = 0.2) +
    geom_vline(xintercept = -1, size = 0.2) +
    geom_hline(yintercept = 1, size = 0.2) +
    geom_hline(yintercept = -1, size = 0.2) +
    geom_point(data = filter(effect_sizes, LMM != "mvLMM"), shape = 16, alpha = 0.25, colour = red, size = 1.5) +
    geom_point(shape = 21, size = 2.5, alpha = 0.8, colour = "black", fill = red) +
    guides() +
    theme(legend.position = "none") +
    theme_linedraw() +
    scale_x_continuous(limits = c(-4, 5), breaks = seq(-4, 5, 1)) +
    scale_y_continuous(limits = c(-4, 3), breaks = seq(-4, 3, 1)) +
    labs(title = "", x = "Effect size (Map axis 1 - Cephalosporins)", y = "Effect size (Map axis 2 - Penicillin)") +
    theme_bw() +
    coord_fixed() +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +
    annotate("label", x = -3.5, y = -4, label = label_text("-", "-"), size = 3.5) +
    annotate("label", x = -3.5, y = -0.75, label = label_text("/", "-"), size = 3.5) +
    annotate("label", x = -3.5, y = 3, label = label_text("+", "-"), size = 3.5) +
    annotate("label", x = 0, y = -4, label = label_text("-", "/"), size = 3.5) +
    annotate("label", x = 0, y = -0.75, label = label_text("/", "/"), size = 3.5) +
    annotate("label", x = 0, y = 3, label = label_text("+", "/"), size = 3.5) +
    annotate("label", x = 4.5, y = -4, label = label_text("-", "+"), size = 3.5) +
    annotate("label", x = 4.5, y = -0.75, label = label_text("/", "+"), size = 3.5) +
    annotate("label", x = 4.5, y = 3, label = label_text("+", "+"), size = 3.5)

  save_plot(p, "recomputed_S18_effect_size_axes_original_style", width = 7.2, height = 6.4)

  list(
    additive_marker_rows = nrow(additive),
    epistasis_p_threshold_rows = nrow(epistasis),
    total_plotted_effects = nrow(effect_sizes)
  )
}

make_s19 <- function() {
  require_exact_package("nVennR", "Supplementary Figure S19")

  sf1 <- read_csv_base(file.path(results_root, "manuscript_outputs", "Supplementary_File_1.csv")) %>%
    mutate(
      rank = case_when(
        Evidence == "Very Strong" ~ 4,
        Evidence == "Strong" ~ 3,
        Evidence == "Moderate" ~ 2,
        Evidence == "Weak" ~ 1,
        TRUE ~ 0
      ),
      Location = str_extract(Substitution, "[0-9]+"),
      PBP_position = paste(PBP, Location, sep = "_")
    ) %>%
    filter(!is.na(Location)) %>%
    arrange(rank) %>%
    group_by(PBP_position) %>%
    slice_tail(n = 1) %>%
    ungroup()

  prior <- read_csv_base(file.path(project_root, "AMRC-repo-files", "Streptococcus pneumoniae analysis", "CDC_GWAS_overlap_TPD.csv")) %>%
    mutate(
      PBP = recode(PBP, "1a" = "1A", "2b" = "2B", "2x" = "2X"),
      Location = as.character(Position),
      PBP_position = paste(PBP, Location, sep = "_"),
      GWAS = tolower(GWAS) == "yes",
      CDC = tolower(CDC) == "yes",
      Laboratory = tolower(Laboratory) == "yes"
    )

  sets <- list(
    `Cartography (Weak)` = sf1 %>% filter(rank >= 1) %>% pull(PBP_position),
    `Cartography (Moderate+)` = sf1 %>% filter(rank >= 2) %>% pull(PBP_position),
    GWAS = prior %>% filter(GWAS) %>% pull(PBP_position),
    CDC = prior %>% filter(CDC) %>% pull(PBP_position),
    Laboratory = prior %>% filter(Laboratory) %>% pull(PBP_position)
  )
  sets <- lapply(sets, unique)

  all_positions <- sort(unique(unlist(sets)))
  membership <- data.frame(PBP_position = all_positions, stringsAsFactors = FALSE)
  for (set_name in names(sets)) {
    membership[[set_name]] <- membership$PBP_position %in% sets[[set_name]]
  }
  write.csv(membership, file.path(out_dir, "recomputed_S19_overlap_membership.csv"), row.names = FALSE)
  set_sizes <- data.frame(set = names(sets), size = vapply(sets, length, integer(1)), stringsAsFactors = FALSE)
  write.csv(set_sizes, file.path(out_dir, "recomputed_S19_overlap_set_sizes.csv"), row.names = FALSE)

  any_previous <- unique(c(sets$GWAS, sets$CDC, sets$Laboratory))
  lab_only <- setdiff(sets$Laboratory, union(sets$GWAS, sets$CDC))
  gwas_only <- setdiff(sets$GWAS, union(sets$Laboratory, sets$CDC))
  summary <- data.frame(
    cartography_weak_plus_positions = length(sets$`Cartography (Weak)`),
    cartography_moderate_plus_positions = length(sets$`Cartography (Moderate+)`),
    previous_any_source_positions = length(any_previous),
    weak_plus_overlap_any_previous = length(intersect(sets$`Cartography (Weak)`, any_previous)),
    weak_plus_not_previous = length(setdiff(sets$`Cartography (Weak)`, any_previous)),
    weak_plus_lab_only_recovered = length(intersect(sets$`Cartography (Weak)`, lab_only)),
    previous_not_recovered_by_weak_plus = length(setdiff(any_previous, sets$`Cartography (Weak)`)),
    previous_not_recovered_gwas_only = length(setdiff(gwas_only, sets$`Cartography (Weak)`)),
    previous_not_recovered_lab_only = length(setdiff(lab_only, sets$`Cartography (Weak)`)),
    moderate_plus_overlap_any_previous = length(intersect(sets$`Cartography (Moderate+)`, any_previous)),
    moderate_plus_not_previous = length(setdiff(sets$`Cartography (Moderate+)`, any_previous)),
    moderate_plus_lab_only_recovered = length(intersect(sets$`Cartography (Moderate+)`, lab_only)),
    previous_not_recovered_by_moderate_plus = length(setdiff(any_previous, sets$`Cartography (Moderate+)`))
  )
  write.csv(summary, file.path(out_dir, "recomputed_S19_overlap_summary.csv"), row.names = FALSE)

  stem <- "recomputed_S19_overlap_venn_original_style"
  svg_path <- file.path(out_dir, paste0(stem, ".svg"))
  nVennR::plotVenn(
    list(
      `Cartography (Weak)` = sets$`Cartography (Weak)`,
      `Cartography (Moderate+)` = sets$`Cartography (Moderate+)`,
      GWAS = sets$GWAS,
      CDC = sets$CDC,
      Laboratory = sets$Laboratory
    ),
    nCycles = 70000,
    showPlot = TRUE,
    opacity = 0.4,
    borderWidth = 2,
    labelRegions = FALSE,
    outFile = svg_path,
    systemShow = FALSE,
    fontScale = 2,
    setColors = c(
      "#FBB4AE",
      red,
      "#4DAF4A",
      "#FFFF33",
      "#377EB8"
    )
  )
  pad_nvennr_svg_canvas(svg_path)
  convert_svg_outputs(svg_path, stem)
  as.list(summary[1, ])
}

requested_figures <- str_split(Sys.getenv("FIGURES", "S17,S18,S19"), ",")[[1]]
requested_figures <- toupper(str_trim(requested_figures))
requested_figures <- requested_figures[nzchar(requested_figures)]

s17 <- if ("S17" %in% requested_figures) make_s17() else NULL
s18 <- if ("S18" %in% requested_figures) make_s18() else NULL
s19 <- if ("S19" %in% requested_figures) make_s19() else NULL

summary_lines <- c(
  "# Recomputed Supplementary Figures - Original Style",
  "",
  "Generated by `generate_recomputed_supplement_figures_original_style.R`.",
  "",
  "This script mirrors the original R/Rmd figure styles while reading recomputed 170-marker outputs.",
  "",
  "Note: this script now refuses to draw non-identical substitutes for the package-specific original plots. S17B requires `ggraph`/`tidygraph`; S19 requires the original `nVennR::plotVenn()` implementation."
)

if (!is.null(s17)) {
  summary_lines <- c(
    summary_lines,
    "",
    "## S17",
    "",
    sprintf("- Output: `%s` and `%s`.", file.path(out_dir, "recomputed_S17_epistasis_original_style.png"), file.path(out_dir, "recomputed_S17_epistasis_original_style.pdf")),
    sprintf("- Panel A: `%s` and `%s`.", file.path(out_dir, "recomputed_S17A_spneumo_effect_size_plot.png"), file.path(out_dir, "recomputed_S17A_spneumo_effect_size_plot.pdf")),
    sprintf("- Panel B: `%s` and `%s`.", file.path(out_dir, "recomputed_S17B_spneumo_epistatic_interaction_network.png"), file.path(out_dir, "recomputed_S17B_spneumo_epistatic_interaction_network.pdf")),
    sprintf("- Panel C: `%s` and `%s`.", file.path(out_dir, "recomputed_S17C_spneumo_epistatic_histogram.png"), file.path(out_dir, "recomputed_S17C_spneumo_epistatic_histogram.pdf")),
    sprintf("- Tested interactions: `%s`.", s17$tested_interactions),
    sprintf("- P-threshold interactions: `%s`.", s17$p_threshold_only_interactions),
    sprintf("- Supported interactions after lower-bound effect filter: `%s`.", s17$supported_interactions),
    sprintf("- Markers with support: `%s`.", s17$markers_with_support),
    sprintf("- Markers with >40 supported interactions: `%s`.", s17$markers_gt_40_supported_interactions),
    sprintf("- Adjusted epistasis threshold: `%s`.", signif(s17$epistasis_threshold, 6))
  )
}

if (!is.null(s18)) {
  summary_lines <- c(
    summary_lines,
    "",
    "## S18",
    "",
    sprintf("- Output: `%s` and `%s`.", file.path(out_dir, "recomputed_S18_effect_size_axes_original_style.png"), file.path(out_dir, "recomputed_S18_effect_size_axes_original_style.pdf")),
    sprintf("- Additive mvLMM marker rows plotted: `%s`.", s18$additive_marker_rows),
    sprintf("- Epistasis p-threshold rows plotted: `%s`.", s18$epistasis_p_threshold_rows),
    sprintf("- Total plotted effects: `%s`.", s18$total_plotted_effects)
  )
}

if (!is.null(s19)) {
  summary_lines <- c(
    summary_lines,
    "",
    "## S19",
    "",
    sprintf("- Output: `%s`, `%s`, and `%s`.", file.path(out_dir, "recomputed_S19_overlap_venn_original_style.png"), file.path(out_dir, "recomputed_S19_overlap_venn_original_style.pdf"), file.path(out_dir, "recomputed_S19_overlap_venn_original_style.svg")),
    sprintf("- Cartography Weak+ positions: `%s`.", s19$cartography_weak_plus_positions),
    sprintf("- Cartography Moderate+ positions: `%s`.", s19$cartography_moderate_plus_positions),
    sprintf("- Previous-source positions: `%s`.", s19$previous_any_source_positions),
    sprintf("- Weak+ overlap with previous sources: `%s`.", s19$weak_plus_overlap_any_previous),
    sprintf("- Weak+ not previous: `%s`.", s19$weak_plus_not_previous)
  )
}
writeLines(summary_lines, file.path(out_dir, "recomputed_supplement_figures_original_style_summary.md"))

cat(paste(summary_lines, collapse = "\n"), "\n")
