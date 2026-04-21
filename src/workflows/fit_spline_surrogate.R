script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

reference_bundle <- readRDS(result_path("01_exact_local_grid.rds"))

spline_fit <- fit_spline_surrogate(
  local_table = reference_bundle$primary$local_table,
  box = reference_bundle$primary$box,
  knots = list(x = 11L, q = 9L, t = 9L, degree = 3L),
  penalty = list(x = 1e-6, q = 1e-6, t = 1e-6),
  nu = reference_bundle$primary$nu
)

saveRDS(spline_fit, result_path("02_spline_fit.rds"))

