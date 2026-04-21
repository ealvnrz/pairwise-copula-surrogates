if (!exists("script_dir", inherits = FALSE)) {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
  } else {
    script_dir <- getwd()
  }
}

discover_repo_root <- function(start_dir) {
  candidates <- unique(
    normalizePath(
      c(
        start_dir,
        file.path(start_dir, ".."),
        file.path(start_dir, "..", ".."),
        getwd()
      ),
      winslash = "/",
      mustWork = FALSE
    )
  )

  for (candidate in candidates) {
    if (
      dir.exists(file.path(candidate, "src", "core")) &&
      dir.exists(file.path(candidate, "src", "utils"))
    ) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  stop("Could not locate the repository root from the current script directory.", call. = FALSE)
}

root_dir <- discover_repo_root(script_dir)
src_dir <- file.path(root_dir, "src")
core_dir <- file.path(src_dir, "core")
utils_dir <- file.path(src_dir, "utils")
results_dir <- file.path(root_dir, "results")
result_dir <- results_dir
results_figures_dir <- file.path(results_dir, "figures")
results_summaries_dir <- file.path(results_dir, "summaries")
results_objects_dir <- file.path(results_dir, "objects")

dir.create(results_figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_summaries_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_objects_dir, recursive = TRUE, showWarnings = FALSE)

result_path <- function(name) {
  ext <- tolower(tools::file_ext(name))
  subdir <- switch(
    ext,
    png = results_figures_dir,
    jpg = results_figures_dir,
    jpeg = results_figures_dir,
    csv = results_summaries_dir,
    txt = results_summaries_dir,
    md = results_summaries_dir,
    rds = results_objects_dir,
    results_dir
  )
  file.path(subdir, name)
}

module_paths <- c(
  file.path(core_dir, "clayton_core.R"),
  file.path(core_dir, "spatial_factor_core.R"),
  file.path(core_dir, "surrogates.R"),
  file.path(core_dir, "reduced_wpl.R"),
  file.path(utils_dir, "diagnostics.R"),
  file.path(utils_dir, "simulation.R")
)

for (module_path in module_paths) {
  source(module_path, local = FALSE)
}
