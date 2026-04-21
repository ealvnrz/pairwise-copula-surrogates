script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

spline_fit <- readRDS(result_path("02_spline_fit.rds"))
cheb_bundle <- readRDS(result_path("03_chebyshev_fit.rds"))

run_opt <- function(objective, start, lower, upper, ...) {
  optimize_reduced_wpl(objective, start = start, lower = lower, upper = upper, ...)
}

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

compute_curvature_proxy <- function(objective, theta, step, ...) {
  hess <- finite_difference_hessian_steps(objective, theta, step = step, ...)
  curvature <- -0.5 * (hess + t(hess))
  eig <- sort(Re(eigen(curvature, symmetric = TRUE, only.values = TRUE)$values))
  positive_definite <- all(is.finite(eig)) && min(eig) > 0
  cov_proxy <- if (positive_definite) solve(curvature) else matrix(NA_real_, nrow(curvature), ncol(curvature))
  list(
    hessian = hess,
    curvature = curvature,
    eigenvalues = eig,
    positive_definite = positive_definite,
    covariance_proxy = cov_proxy
  )
}

sim <- simulate_clayton_beta_example(
  n = 120L,
  beta = 0,
  delta = 6,
  b = 0.35,
  nu = 6,
  seed = 20263300L
)

start <- c(delta = 5.5, b = 0.30)
lower <- c(1.0, 0.05)
upper <- c(20.0, 1.5)
step <- c(0.05, 0.01)

replicate_rows <- lapply(seq_len(10L), function(rep_id) {
  rep_sim <- simulate_clayton_beta_example(
    n = 120L,
    beta = sim$true$beta,
    delta = sim$true$delta,
    b = sim$true$b,
    nu = sim$true$nu,
    seed = 20263300L + rep_id
  )
  rep_data <- prepare_reduced_clayton_data(rep_sim$y, rep_sim$coords, rep_sim$X, m = 2L)

  fits <- list(
    reference = run_opt(
      objective = wpl_reference_reduced,
      start = start,
      lower = lower,
      upper = upper,
      data = rep_data,
      fixed_beta = rep_sim$true$beta,
      fixed_nu = rep_sim$true$nu
    ),
    spline = run_opt(
      objective = wpl_surrogate_reduced,
      start = start,
      lower = lower,
      upper = upper,
      data = rep_data,
      surrogate_fit = spline_fit,
      fixed_beta = rep_sim$true$beta,
      fixed_nu = rep_sim$true$nu
    ),
    chebyshev = run_opt(
      objective = wpl_surrogate_reduced,
      start = start,
      lower = lower,
      upper = upper,
      data = rep_data,
      surrogate_fit = cheb_bundle$fit,
      fixed_beta = rep_sim$true$beta,
      fixed_nu = rep_sim$true$nu
    )
  )

  proxies <- list(
    reference = if (identical(fits$reference$convergence, 0L)) {
      compute_curvature_proxy(
        function(theta, data, fixed_beta, fixed_nu) {
          wpl_reference_reduced(theta, data = data, fixed_beta = fixed_beta, fixed_nu = fixed_nu)
        },
        theta = fits$reference$par,
        step = step,
        data = rep_data,
        fixed_beta = rep_sim$true$beta,
        fixed_nu = rep_sim$true$nu
      )
    } else NULL,
    spline = if (identical(fits$spline$convergence, 0L)) {
      compute_curvature_proxy(
        function(theta, data, surrogate_fit, fixed_beta, fixed_nu) {
          wpl_surrogate_reduced(theta, data = data, surrogate_fit = surrogate_fit, fixed_beta = fixed_beta, fixed_nu = fixed_nu)
        },
        theta = fits$spline$par,
        step = step,
        data = rep_data,
        surrogate_fit = spline_fit,
        fixed_beta = rep_sim$true$beta,
        fixed_nu = rep_sim$true$nu
      )
    } else NULL,
    chebyshev = if (identical(fits$chebyshev$convergence, 0L)) {
      compute_curvature_proxy(
        function(theta, data, surrogate_fit, fixed_beta, fixed_nu) {
          wpl_surrogate_reduced(theta, data = data, surrogate_fit = surrogate_fit, fixed_beta = fixed_beta, fixed_nu = fixed_nu)
        },
        theta = fits$chebyshev$par,
        step = step,
        data = rep_data,
        surrogate_fit = cheb_bundle$fit,
        fixed_beta = rep_sim$true$beta,
        fixed_nu = rep_sim$true$nu
      )
    } else NULL
  )

  ref_cov <- if (!is.null(proxies$reference)) proxies$reference$covariance_proxy else matrix(NA_real_, 2L, 2L)
  ref_se <- if (all(is.finite(ref_cov))) sqrt(diag(ref_cov)) else rep(NA_real_, 2L)

  do.call(
    rbind,
    lapply(c("reference", "spline", "chebyshev"), function(method) {
      fit <- fits[[method]]
      proxy <- proxies[[method]]
      cov_proxy <- if (!is.null(proxy)) proxy$covariance_proxy else matrix(NA_real_, 2L, 2L)
      se_proxy <- if (all(is.finite(cov_proxy))) sqrt(diag(cov_proxy)) else rep(NA_real_, 2L)
      rel_cov_gap <- if (method == "reference" || !all(is.finite(ref_cov)) || !all(is.finite(cov_proxy))) {
        NA_real_
      } else {
        sqrt(sum((cov_proxy - ref_cov)^2)) / sqrt(sum(ref_cov^2))
      }
      data.frame(
        replicate = rep_id,
        method = method,
        converged = identical(fit$convergence, 0L),
        positive_definite = !is.null(proxy) && isTRUE(proxy$positive_definite),
        delta_hat = unname(fit$par[1]),
        b_hat = unname(fit$par[2]),
        lambda_min = if (!is.null(proxy)) min(proxy$eigenvalues) else NA_real_,
        lambda_max = if (!is.null(proxy)) max(proxy$eigenvalues) else NA_real_,
        condition_number = if (!is.null(proxy) && proxy$positive_definite) max(proxy$eigenvalues) / min(proxy$eigenvalues) else NA_real_,
        delta_se_proxy = se_proxy[1],
        b_se_proxy = se_proxy[2],
        relative_covariance_gap_vs_reference = rel_cov_gap,
        abs_delta_se_gap_vs_reference = if (method == "reference") NA_real_ else abs(se_proxy[1] - ref_se[1]),
        abs_b_se_gap_vs_reference = if (method == "reference") NA_real_ else abs(se_proxy[2] - ref_se[2]),
        stringsAsFactors = FALSE
      )
    })
  )
})

replicate_df <- do.call(rbind, replicate_rows)

summary_df <- do.call(
  rbind,
  lapply(split(replicate_df, replicate_df$method), function(sub_df) {
    data.frame(
      method = sub_df$method[1],
      n = nrow(sub_df),
      convergence_rate = mean(sub_df$converged),
      positive_definite_rate = mean(sub_df$positive_definite),
      mean_lambda_min = mean(sub_df$lambda_min, na.rm = TRUE),
      mean_lambda_max = mean(sub_df$lambda_max, na.rm = TRUE),
      mean_condition_number = mean(sub_df$condition_number, na.rm = TRUE),
      mean_delta_se_proxy = mean(sub_df$delta_se_proxy, na.rm = TRUE),
      mean_b_se_proxy = mean(sub_df$b_se_proxy, na.rm = TRUE),
      mean_relative_covariance_gap_vs_reference = mean(sub_df$relative_covariance_gap_vs_reference, na.rm = TRUE),
      mean_abs_delta_se_gap_vs_reference = mean(sub_df$abs_delta_se_gap_vs_reference, na.rm = TRUE),
      mean_abs_b_se_gap_vs_reference = mean(sub_df$abs_b_se_gap_vs_reference, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
rownames(summary_df) <- NULL
summary_df <- summary_df[match(c("reference", "spline", "chebyshev"), summary_df$method), , drop = FALSE]

write.csv(replicate_df, result_path("example1_uncertainty_proxy_replicates.csv"), row.names = FALSE)
write.csv(summary_df, result_path("example1_uncertainty_proxy_summary.csv"), row.names = FALSE)

