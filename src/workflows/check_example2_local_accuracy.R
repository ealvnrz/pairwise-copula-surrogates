script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

factor_bundle <- readRDS(result_path("25_factorcopula_example_box_and_refinement.rds"))
factor_spec <- factor_bundle$factor_spec
box <- readRDS(result_path("factorcopula_example_box.rds"))
spline_fit <- readRDS(result_path("factorcopula_example_spline_fit.rds"))
cheb_bundle <- readRDS(result_path("factorcopula_example_chebyshev_fit.rds"))
highres_spec <- modifyList(factor_spec, list(quadrature_n = 80L))

evaluation_grid <- sample_transformed_points_by_region(
  box = box,
  n_per_region = 200L,
  seed = 20268220L
)

reference_fun <- function(z) {
  reference_transformed_factorcopula(z[1], z[2], z[3], fixed_factor = factor_spec, log = TRUE)
}

spline_errors <- surrogate_local_error_summary_generic(
  fit = spline_fit,
  evaluation_grid = evaluation_grid,
  reference_fun = reference_fun
)
cheb_errors <- surrogate_local_error_summary_generic(
  fit = cheb_bundle$fit,
  evaluation_grid = evaluation_grid,
  reference_fun = reference_fun
)

pointwise <- rbind(
  transform(spline_errors$pointwise, method = "spline", label = "selected"),
  transform(cheb_errors$pointwise, method = "chebyshev", label = "selected")
)

stability_df <- within(
  classify_transformed_region(evaluation_grid, box),
  {
    baseline <- reference_transformed_factorcopula(x, q, t, fixed_factor = factor_spec, log = TRUE)
    highres <- reference_transformed_factorcopula(x, q, t, fixed_factor = highres_spec, log = TRUE)
    abs_diff_highres = abs(baseline - highres)
  }
)

by_region <- do.call(
  rbind,
  lapply(split(pointwise, interaction(pointwise$method, pointwise$region, drop = TRUE)), function(df) {
    data.frame(
      method = df$method[1],
      region = as.character(df$region[1]),
      n = nrow(df),
      mean_value_error = mean(df$value_error),
      median_value_error = median(df$value_error),
      max_value_error = max(df$value_error),
      mean_grad_error = mean(df$grad_error),
      mean_hess_error = mean(df$hess_error),
      stringsAsFactors = FALSE
    )
  })
)
rownames(by_region) <- NULL

overall <- do.call(
  rbind,
  lapply(split(pointwise, pointwise$method), function(df) {
    data.frame(
      method = df$method[1],
      label = "selected",
      mean_value_error = mean(df$value_error),
      mean_grad_error = mean(df$grad_error),
      mean_hess_error = mean(df$hess_error),
      sup_value_error = max(df$value_error),
      sup_grad_error = max(df$grad_error),
      sup_hess_error = max(df$hess_error),
      quadrature_mean_abs_diff_highres = mean(stability_df$abs_diff_highres),
      quadrature_max_abs_diff_highres = max(stability_df$abs_diff_highres),
      stringsAsFactors = FALSE
    )
  })
)
rownames(overall) <- NULL

write.csv(pointwise, result_path("factorcopula_example_local_accuracy_pointwise.csv"), row.names = FALSE)
write.csv(by_region, result_path("factorcopula_example_local_error_by_region.csv"), row.names = FALSE)
write.csv(overall, result_path("factorcopula_example_local_accuracy.csv"), row.names = FALSE)
save_surrogate_error_boxplot(
  pointwise_df = pointwise,
  file = result_path("fig_factorcopula_example_local_error_boxplot.png")
)

