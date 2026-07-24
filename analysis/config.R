# =============================================================================
# AMR-cartography — shared analysis configuration
# -----------------------------------------------------------------------------
# Single source of truth for (a) the input data directory and (b) the canonical
# map transform constants used across the numbered analysis scripts. Added
# 2026-07-01 to remove hardcoded machine-specific paths and duplicated constants
# (audit findings #1 and #3); it does not change any analysis logic.
#
# The numbered scripts load this file with:
#     source(Sys.getenv("AMRC_CONFIG",
#                        unset = path.expand("~/AMR-cartography/analysis/config.R")))
# which resolves on both the original Mac (/Users/ajb306/AMR-cartography) and the
# farm (~/AMR-cartography) without depending on the working directory. Override
# the location with the AMRC_CONFIG environment variable if the repo lives
# elsewhere.
# =============================================================================

# ---- Input data directory ---------------------------------------------------
# Defaults to the original author path so existing local workflows are
# unchanged. Point elsewhere (farm, another machine) by exporting AMRC_DATA_DIR
# before launching R, e.g.
#     export AMRC_DATA_DIR=/path/to/AMR-cartography-results/data
# or from R: Sys.setenv(AMRC_DATA_DIR = "/path/to/data")
AMRC_DATA_DIR <- Sys.getenv(
  "AMRC_DATA_DIR",
  unset = "/Users/ajb306/AMR-cartography-results/data"
)

# Build a path inside the data directory, e.g. amrc_data_path("MIC_table_Spneumoniae.csv")
amrc_data_path <- function(...) file.path(AMRC_DATA_DIR, ...)

# ---- Other external working directories -------------------------------------
# Two further author-machine directories referenced by the numbered scripts.
# Both default to the original author paths (so existing local runs are
# unchanged) and are env-overridable exactly like AMRC_DATA_DIR above. Keeping
# them here means no numbered script contains a machine-specific absolute path.

# Python mvLMM project directory: holds the intermediate CSVs written by script
# 31 (…/pythonProject1) and re-read by scripts 30/32.
AMRC_PYPROJECT_DIR <- Sys.getenv(
  "AMRC_PYPROJECT_DIR",
  unset = "/Users/ajb306/PycharmProjects/pythonProject1"
)
amrc_pyproject_path <- function(...) file.path(AMRC_PYPROJECT_DIR, ...)

# Pneumococcus analysis output directory: holds the ranking/overlap support
# tables (Sub_effect_sizes_mv_pneumo.csv, Sig_AA_subs_* etc.) written/read by
# scripts 17, 30 and 32.
AMRC_PNEUMO_ANALYSIS_DIR <- Sys.getenv(
  "AMRC_PNEUMO_ANALYSIS_DIR",
  unset = "/Users/ajb306/Google Drive/PhD - Andrew Balmer/Chapters/AMR cartography/MIC data/Streptococcus pneumoniae analysis"
)
amrc_pneumo_path <- function(...) file.path(AMRC_PNEUMO_ANALYSIS_DIR, ...)

# ---- Canonical map transform constants (audit finding #3) -------------------
# The phenotype/genetic SMACOF maps are rotated and dilated onto an
# interpretable MIC-unit scale. The goodness-of-fit scripts (05-07) DERIVE these
# by regressing map distance on measured (table) distance; downstream scripts
# use the frozen values below. If the maps are ever regenerated, update the
# values here and re-derive in 05-07 — every downstream script reads them from
# this file, so they can no longer drift out of sync.
AMRC_MAP_ROTATION_DEGREES <- 326          # phenotype map rotation, in degrees
AMRC_PHEN_DILATION_SLOPE  <- 0.1842996    # phenotype map: MIC units per map unit = 1 / slope
AMRC_GEN_DILATION_SLOPE   <- 0.01814108   # genetic map dilation slope
