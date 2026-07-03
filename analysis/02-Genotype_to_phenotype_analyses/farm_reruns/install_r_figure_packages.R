#!/usr/bin/env Rscript

lib <- Sys.getenv(
  "R_LIBS_USER",
  file.path(Sys.getenv("HOME"), "R", "x86_64-pc-linux-gnu-library", paste0(R.version$major, ".", R.version$minor))
)
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

options(repos = c(CRAN = "https://cloud.r-project.org"))

install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, lib = lib, dependencies = NA)
  }
}

for (package in c("ggraph", "tidygraph")) {
  install_if_missing(package)
}

if (!requireNamespace("nVennR", quietly = TRUE)) {
  install_if_missing("remotes")
  remotes::install_github("vqf/nVennR", lib = lib, upgrade = "never")
}

missing <- c("ggraph", "tidygraph", "nVennR")[
  !vapply(c("ggraph", "tidygraph", "nVennR"), requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing) > 0) {
  stop("Missing required original figure packages after installation: ", paste(missing, collapse = ", "))
}

cat("Original figure packages available in:\n", paste(.libPaths(), collapse = "\n"), "\n")
