# =============================================================================
# ABT-RISE: Config Loader
#
# Reads clif_config.json at repo root and exposes:
#   site_id   -- from config "site_name"
#   data_dir  -- from "abtrise_input_dir"  (default: ./output_phi/analysis)
#   out_dir   -- from "abtrise_output_dir" (default: ./output_to_share)
#
# Sourced by ABTRISE_run_all.R and ABTRISE_01_setup_c.R.
# Sites: edit clif_config.json at repo root. No edits needed in R code.
# =============================================================================

# --- Auto-install helper -----------------------------------------------------
# The pipeline runs under `Rscript --vanilla`, so no CRAN mirror is set by
# default (getOption("repos") is the unresolved "@CRAN@" placeholder); pass one
# explicitly. Installs only the missing packages -- base packages such as
# `splines` are already present and skipped -- then verifies and fails loudly
# with an actionable message if an install could not succeed.
ensure_packages <- function(pkgs, repos = "https://cloud.r-project.org") {
  missing <- setdiff(pkgs, rownames(installed.packages()))
  if (length(missing) > 0) {
    message("Installing missing R package(s): ", paste(missing, collapse = ", "))
    install.packages(missing, repos = repos)
  }
  still_missing <- setdiff(pkgs, rownames(installed.packages()))
  if (length(still_missing) > 0) {
    stop("Could not install R package(s): ", paste(still_missing, collapse = ", "),
         ".\n  Check internet access / build tools, or install manually with ",
         "install.packages(c(",
         paste0('"', still_missing, '"', collapse = ", "), ")).",
         call. = FALSE)
  }
}

# config.R is sourced before anything else on every path, so bootstrap its own
# dependencies here before the fail-fast library() calls below.
ensure_packages(c("here", "jsonlite"))

suppressPackageStartupMessages({
  library(here)
  library(jsonlite)
})

.abtrise_resolve_path <- function(path, default) {
  p <- if (is.null(path) || !nzchar(path)) default else path
  if (substr(p, 1, 1) %in% c("/", "~")) return(p)
  p <- sub("^\\./", "", p)
  here::here(p)
}

.abtrise_config_path <- here::here("clif_config.json")

if (!file.exists(.abtrise_config_path)) {
  stop(
    "clif_config.json not found at repo root: ", .abtrise_config_path, "\n",
    "Copy clif_config_template.json to clif_config.json and edit it."
  )
}

.abtrise_cfg <- jsonlite::fromJSON(.abtrise_config_path)

if (is.null(.abtrise_cfg$site_name) || !nzchar(.abtrise_cfg$site_name)) {
  stop("clif_config.json: 'site_name' is missing or empty.")
}

site_id  <- .abtrise_cfg$site_name
data_dir <- .abtrise_resolve_path(.abtrise_cfg$abtrise_input_dir,
                                  "./output_phi/analysis")
out_dir  <- .abtrise_resolve_path(.abtrise_cfg$abtrise_output_dir,
                                  "./output_to_share")

if (!dir.exists(data_dir)) {
  stop("abtrise_input_dir does not exist: ", data_dir)
}

.abtrise_required <- c("file1_person_period.parquet",
                       "file2_hospitalization_level.parquet")
.abtrise_missing  <- .abtrise_required[
  !file.exists(file.path(data_dir, .abtrise_required))
]
if (length(.abtrise_missing) > 0) {
  stop("Required input file(s) not found in ", data_dir, ": ",
       paste(.abtrise_missing, collapse = ", "))
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("── ABT-RISE config ──────────────────────────────────────────\n")
cat("  site_id : ", site_id,  "\n", sep = "")
cat("  data_dir: ", data_dir, "\n", sep = "")
cat("  out_dir : ", out_dir,  "\n", sep = "")
cat("─────────────────────────────────────────────────────────────\n\n")
