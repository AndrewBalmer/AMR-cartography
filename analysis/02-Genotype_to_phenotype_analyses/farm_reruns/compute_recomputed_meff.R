#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(poolr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: compute_recomputed_meff.R DATA_DIR INTERACTION_FILE OUT_JSON", call. = FALSE)
}

data_dir <- args[[1]]
interaction_file <- args[[2]]
out_json <- args[[3]]

markers <- read.csv(file.path(data_dir, "S.pneumo_map_dummy_gen_test_markers.csv"), check.names = FALSE)
interactions <- read.csv(interaction_file, check.names = FALSE)

additive_cor <- cor(markers, method = "pearson")
additive_meff <- unname(meff(additive_cor, method = "galwey"))

interaction_matrix <- matrix(0, nrow = nrow(markers), ncol = nrow(interactions))
for (i in seq_len(nrow(interactions))) {
  interaction_matrix[, i] <- markers[[interactions$parent_a[[i]]]] * markers[[interactions$parent_b[[i]]]]
}
colnames(interaction_matrix) <- interactions$interaction
epistasis_cor <- cor(interaction_matrix, method = "pearson")
epistasis_meff <- unname(meff(epistasis_cor, method = "galwey"))

json <- sprintf(
  paste0(
    "{\n",
    "  \"method\": \"poolr::meff(cor(...), method = \\\"galwey\\\")\",\n",
    "  \"additive_marker_count\": %d,\n",
    "  \"additive_meff\": %.17g,\n",
    "  \"additive_meff_floor\": %d,\n",
    "  \"epistasis_interaction_count\": %d,\n",
    "  \"epistasis_meff\": %.17g,\n",
    "  \"epistasis_meff_floor\": %d\n",
    "}\n"
  ),
  ncol(markers),
  additive_meff,
  floor(additive_meff),
  nrow(interactions),
  epistasis_meff,
  floor(epistasis_meff)
)

dir.create(dirname(out_json), recursive = TRUE, showWarnings = FALSE)
writeLines(json, out_json)
cat(json)
