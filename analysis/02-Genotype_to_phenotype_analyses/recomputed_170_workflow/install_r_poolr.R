#!/usr/bin/env Rscript

repo <- Sys.getenv("CRAN_REPO", "https://cloud.r-project.org")
user_lib <- Sys.getenv("R_LIBS_USER")
if (!nzchar(user_lib)) {
  user_lib <- file.path(Sys.getenv("HOME"), "R", paste0(R.version$platform, "-library"), paste(R.version$major, R.version$minor, sep = "."))
}

dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(user_lib, .libPaths())))

if (!requireNamespace("poolr", quietly = TRUE)) {
  install.packages("poolr", lib = user_lib, repos = repo)
}

suppressPackageStartupMessages(library(poolr))
cat("poolr installed and loadable from:\n")
cat(system.file(package = "poolr"), "\n")
