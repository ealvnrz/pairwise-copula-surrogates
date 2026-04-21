script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

reference_bundle <- readRDS(result_path("01_exact_local_grid.rds"))

evaluation_grid <- sample_transformed_points_by_region(
  box = reference_bundle$primary$box,
  n_per_region = 150L,
  seed = 20260410L
)

truncation <- reference_truncation_stability(
  grid_df = evaluation_grid,
  nu = reference_bundle$primary$nu,
  box = reference_bundle$primary$box,
  truncation_levels = c(300L, 400L, 500L, 700L, 900L)
)

write.csv(truncation$truncation_pointwise, result_path("truncation_stability_pointwise.csv"), row.names = FALSE)
write.csv(truncation$by_region, result_path("truncation_stability_summary.csv"), row.names = FALSE)

save_truncation_stability_summary_plot(
  truncation$by_region,
  result_path("fig_truncation_stability_boxplot.png")
)

