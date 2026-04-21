script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

reference_bundle <- readRDS(result_path("01_exact_local_grid.rds"))
spline_fit <- readRDS(result_path("02_spline_fit.rds"))
cheb_bundle <- readRDS(result_path("03_chebyshev_fit.rds"))

evaluation_grid <- expand.grid(
  x = seq(reference_bundle$primary$box$x[1] + 0.25, reference_bundle$primary$box$x[2] - 0.25, length.out = 7),
  q = seq(max(0.25, reference_bundle$primary$box$q[1] + 0.25), reference_bundle$primary$box$q[2] - 0.25, length.out = 6),
  t = seq(reference_bundle$primary$box$t[1] + 0.10, reference_bundle$primary$box$t[2] - 0.10, length.out = 6),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

reference_external <- reference_vs_geomodels_summary(evaluation_grid, nu = reference_bundle$primary$nu, box = reference_bundle$primary$box)
spline_diag <- surrogate_local_error_summary(spline_fit, nu = reference_bundle$primary$nu, evaluation_grid = evaluation_grid)
cheb_diag <- surrogate_local_error_summary(cheb_bundle$fit, nu = reference_bundle$primary$nu, evaluation_grid = evaluation_grid)

spline_pointwise <- classify_transformed_region(spline_diag$pointwise, reference_bundle$primary$box)
cheb_pointwise <- classify_transformed_region(cheb_diag$pointwise, reference_bundle$primary$box)

local_error_summary <- data.frame(
  method = c("spline", "chebyshev"),
  sup_value_error = c(spline_diag$summary$sup_value_error, cheb_diag$summary$sup_value_error),
  mean_value_error = c(spline_diag$summary$mean_value_error, cheb_diag$summary$mean_value_error),
  sup_grad_error = c(spline_diag$summary$sup_grad_error, cheb_diag$summary$sup_grad_error),
  mean_grad_error = c(spline_diag$summary$mean_grad_error, cheb_diag$summary$mean_grad_error),
  sup_hess_error = c(spline_diag$summary$sup_hess_error, cheb_diag$summary$sup_hess_error),
  mean_hess_error = c(spline_diag$summary$mean_hess_error, cheb_diag$summary$mean_hess_error),
  stringsAsFactors = FALSE
)

local_error_by_region <- do.call(
  rbind,
  lapply(
    names(list(spline = spline_pointwise, chebyshev = cheb_pointwise)),
    function(method_name) {
      df <- switch(method_name, spline = spline_pointwise, chebyshev = cheb_pointwise)
      out <- do.call(
        rbind,
        lapply(split(df, df$region), function(sub_df) {
          data.frame(
            region = as.character(sub_df$region[1]),
            mean_value_error = mean(sub_df$value_error),
            max_value_error = max(sub_df$value_error),
            mean_grad_error = mean(sub_df$grad_error),
            max_grad_error = max(sub_df$grad_error),
            mean_hess_error = mean(sub_df$hess_error),
            max_hess_error = max(sub_df$hess_error),
            stringsAsFactors = FALSE
          )
        })
      )
      out$method <- method_name
      out
    }
  )
)
local_error_by_region <- local_error_by_region[, c("method", "region", "mean_value_error", "max_value_error", "mean_grad_error", "max_grad_error", "mean_hess_error", "max_hess_error")]

error_plot_df <- rbind(
  transform(spline_pointwise, method = "spline"),
  transform(cheb_pointwise, method = "chebyshev")
)
write.csv(reference_external$pointwise, result_path("04_reference_vs_geomodels_pointwise.csv"), row.names = FALSE)
write.csv(reference_external$by_region, result_path("geomodels_comparison_summary.csv"), row.names = FALSE)
write.csv(spline_pointwise, result_path("04_spline_vs_reference_pointwise.csv"), row.names = FALSE)
write.csv(cheb_pointwise, result_path("04_chebyshev_vs_reference_pointwise.csv"), row.names = FALSE)
write.csv(local_error_summary, result_path("04_local_error_summary.csv"), row.names = FALSE)
write.csv(local_error_by_region, result_path("04_local_error_by_region.csv"), row.names = FALSE)

save_surrogate_error_boxplot(
  error_plot_df,
  result_path("fig_surrogate_value_error_boxplot.png")
)

