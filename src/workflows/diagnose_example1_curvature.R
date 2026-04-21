script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

spline_fit <- readRDS(result_path("02_spline_fit.rds"))
cheb_bundle <- readRDS(result_path("03_chebyshev_fit.rds"))
comparison <- readRDS(result_path("05_reduced_wpl_comparison.rds"))

data_reduced <- prepare_reduced_clayton_data(
  comparison$simulation$y,
  comparison$simulation$coords,
  comparison$simulation$X,
  m = 2L
)

reduced_optima <- comparison$primary_optima
reduced_reference <- as.numeric(reduced_optima[reduced_optima$method == "reference", c("delta", "b")])

finite_difference_hessian_steps <- function(fun, x, step, ...) {
  x <- as.numeric(x)
  step <- as.numeric(step)
  if (length(step) == 1L) {
    step <- rep(step, length(x))
  }
  hess <- matrix(0, nrow = length(x), ncol = length(x))
  for (i in seq_along(x)) {
    for (j in seq_along(x)) {
      ei <- rep(0, length(x))
      ej <- rep(0, length(x))
      ei[i] <- step[i]
      ej[j] <- step[j]
      hess[i, j] <- (
        fun(x + ei + ej, ...) -
          fun(x + ei - ej, ...) -
          fun(x - ei + ej, ...) +
          fun(x - ei - ej, ...)
      ) / (4 * step[i] * step[j])
    }
  }
  hess
}

summarize_curvature <- function(method, hess) {
  curvature <- -0.5 * (hess + t(hess))
  eig <- sort(Re(eigen(curvature, symmetric = TRUE, only.values = TRUE)$values))
  data.frame(
    context = "reduced_benchmark",
    method = method,
    lambda_min = min(eig),
    lambda_max = max(eig),
    condition_number = max(abs(eig)) / min(abs(eig)),
    stringsAsFactors = FALSE
  )
}

curvature_summary <- rbind(
  summarize_curvature(
    "reference",
    finite_difference_hessian_steps(
      function(theta) {
        wpl_reference_reduced(theta, data = data_reduced, fixed_beta = comparison$simulation$true$beta, fixed_nu = comparison$simulation$true$nu)
      },
      reduced_reference,
      step = c(0.05, 0.01)
    )
  ),
  summarize_curvature(
    "spline",
    finite_difference_hessian_steps(
      function(theta) {
        wpl_surrogate_reduced(theta, data = data_reduced, surrogate_fit = spline_fit, fixed_beta = comparison$simulation$true$beta, fixed_nu = comparison$simulation$true$nu)
      },
      reduced_reference,
      step = c(0.05, 0.01)
    )
  ),
  summarize_curvature(
    "chebyshev",
    finite_difference_hessian_steps(
      function(theta) {
        wpl_surrogate_reduced(theta, data = data_reduced, surrogate_fit = cheb_bundle$fit, fixed_beta = comparison$simulation$true$beta, fixed_nu = comparison$simulation$true$nu)
      },
      reduced_reference,
      step = c(0.05, 0.01)
    )
  )
)

write.csv(curvature_summary, result_path("criterion_curvature_summary.csv"), row.names = FALSE)

stale_outputs <- c(
  "criterion_geometry_profiles.csv",
  "criterion_geometry_profiles.rds",
  "fig_geometry_profiles.png"
)
for (output_name in stale_outputs) {
  output_path <- result_path(output_name)
  if (file.exists(output_path)) {
    file.remove(output_path)
  }
}

