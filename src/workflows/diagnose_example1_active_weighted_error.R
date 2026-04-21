script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

spline_fit <- readRDS(result_path("02_spline_fit.rds"))
cheb_bundle <- readRDS(result_path("03_chebyshev_fit.rds"))

stress_corner_q_cut <- 5
stress_corner_t_cut <- 1.875

run_opt_with_time <- function(objective, start, lower, upper, ...) {
  elapsed <- system.time({
    opt <- optimize_reduced_wpl(objective, start = start, lower = lower, upper = upper, ...)
  })[["elapsed"]]
  list(opt = opt, elapsed = elapsed)
}

is_reference_converged <- function(opt_obj) {
  identical(opt_obj$convergence, 0L) &&
    all(is.finite(opt_obj$par)) &&
    is.finite(opt_obj$value)
}

compute_active_pair_local_metrics <- function(component, surrogate_fit) {
  idx <- component$index
  reference_local <- reference_local_clayton(
    component$u[idx[, "i"]],
    component$u[idx[, "j"]],
    component$r,
    nu = component$nu,
    log = TRUE
  )
  surrogate_local <- predict_local_surrogate(
    surrogate_fit,
    component$u[idx[, "i"]],
    component$u[idx[, "j"]],
    component$r,
    nu = component$nu
  )
  abs_error <- abs(surrogate_local - reference_local)
  transformed <- transform_clayton_inputs(
    component$u[idx[, "i"]],
    component$u[idx[, "j"]],
    component$r
  )
  stress_mask <- transformed$q < stress_corner_q_cut & transformed$t > stress_corner_t_cut
  total_weight <- sum(component$weights)

  list(
    weighted_mae = sum(component$weights * abs_error) / total_weight,
    weighted_rmse = sqrt(sum(component$weights * abs_error^2) / total_weight),
    max_abs_err = max(abs_error),
    stress_weight_share = sum(component$weights[stress_mask]) / total_weight,
    total_weight = total_weight,
    n_pairs = length(component$weights)
  )
}

make_metric_row <- function(benchmark, replicate, seed, method, reference_fit, metrics, theta_ref) {
  data.frame(
    benchmark = benchmark,
    replicate = replicate,
    seed = seed,
    method = method,
    reference_converged = is_reference_converged(reference_fit$opt),
    reference_convergence_code = reference_fit$opt$convergence,
    reference_optimization_time = reference_fit$elapsed,
    reference_objective = if (is_reference_converged(reference_fit$opt)) unname(reference_fit$opt$value) else NA_real_,
    weighted_mae = metrics$weighted_mae,
    weighted_rmse = metrics$weighted_rmse,
    max_abs_err = metrics$max_abs_err,
    stress_weight_share = metrics$stress_weight_share,
    total_weight = metrics$total_weight,
    n_pairs = metrics$n_pairs,
    theta_ref_1 = theta_ref[1],
    theta_ref_2 = theta_ref[2],
    theta_ref_3 = if (length(theta_ref) >= 3L) theta_ref[3] else NA_real_,
    theta_ref_4 = if (length(theta_ref) >= 4L) theta_ref[4] else NA_real_,
    stringsAsFactors = FALSE
  )
}

make_missing_row <- function(benchmark, replicate, seed, method, reference_fit) {
  data.frame(
    benchmark = benchmark,
    replicate = replicate,
    seed = seed,
    method = method,
    reference_converged = FALSE,
    reference_convergence_code = reference_fit$opt$convergence,
    reference_optimization_time = reference_fit$elapsed,
    reference_objective = NA_real_,
    weighted_mae = NA_real_,
    weighted_rmse = NA_real_,
    max_abs_err = NA_real_,
    stress_weight_share = NA_real_,
    total_weight = NA_real_,
    n_pairs = NA_integer_,
    theta_ref_1 = NA_real_,
    theta_ref_2 = NA_real_,
    theta_ref_3 = NA_real_,
    theta_ref_4 = NA_real_,
    stringsAsFactors = FALSE
  )
}

quantile95 <- function(x) {
  stats::quantile(x, probs = 0.95, names = FALSE, na.rm = TRUE, type = 7)
}

summarize_weighted_error_results <- function(replicate_df) {
  benchmark_levels <- unique(replicate_df$benchmark)
  method_levels <- c("spline", "chebyshev")
  rows <- vector("list", length(benchmark_levels) * length(method_levels))
  row_id <- 1L

  for (benchmark_name in benchmark_levels) {
    benchmark_df <- replicate_df[replicate_df$benchmark == benchmark_name, , drop = FALSE]
    total_datasets <- length(unique(benchmark_df$replicate))
    used_replicates <- unique(benchmark_df$replicate[benchmark_df$reference_converged])
    n_used <- length(used_replicates)
    reference_convergence_rate <- n_used / total_datasets

    for (method_name in method_levels) {
      method_df <- benchmark_df[
        benchmark_df$method == method_name & benchmark_df$reference_converged,
        ,
        drop = FALSE
      ]
      rows[[row_id]] <- data.frame(
        benchmark = benchmark_name,
        method = method_name,
        n_total = total_datasets,
        n_used = n_used,
        reference_convergence_rate = reference_convergence_rate,
        mean_weighted_mae = mean(method_df$weighted_mae, na.rm = TRUE),
        median_weighted_mae = stats::median(method_df$weighted_mae, na.rm = TRUE),
        p95_weighted_mae = quantile95(method_df$weighted_mae),
        mean_weighted_rmse = mean(method_df$weighted_rmse, na.rm = TRUE),
        mean_max_abs_err = mean(method_df$max_abs_err, na.rm = TRUE),
        mean_stress_weight_share = mean(method_df$stress_weight_share, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
    }
  }

  do.call(rbind, rows)
}

format_sci <- function(x, digits = 2L) {
  ifelse(is.na(x), "--", formatC(x, format = "e", digits = digits))
}

format_share <- function(x, digits = 2L) {
  ifelse(is.na(x), "--", sprintf(paste0("%.", digits, "f\\%%"), 100 * x))
}

evaluate_reduced_benchmark <- function(spline_fit, chebyshev_fit) {
  true_spec <- list(beta = 0, delta = 6, b = 0.35, nu = 6)
  seeds <- 20260410L:20260459L
  start <- c(delta = 5.5, b = 0.30)
  lower <- c(1.0, 0.05)
  upper <- c(20.0, 1.5)

  rows <- lapply(seq_along(seeds), function(rep_id) {
    seed <- seeds[[rep_id]]
    sim <- simulate_clayton_beta_example(
      n = 120L,
      beta = true_spec$beta,
      delta = true_spec$delta,
      b = true_spec$b,
      nu = true_spec$nu,
      seed = seed
    )
    data <- prepare_reduced_clayton_data(sim$y, sim$coords, sim$X, m = 2L)
    reference_fit <- run_opt_with_time(
      objective = wpl_reference_reduced,
      start = start,
      lower = lower,
      upper = upper,
      data = data,
      fixed_beta = sim$true$beta,
      fixed_nu = sim$true$nu
    )

    if (!is_reference_converged(reference_fit$opt)) {
      return(rbind(
        make_missing_row("reduced", rep_id, seed, "spline", reference_fit),
        make_missing_row("reduced", rep_id, seed, "chebyshev", reference_fit)
      ))
    }

    comp <- build_reduced_clayton_components(
      reference_fit$opt$par,
      data = data,
      fixed_beta = sim$true$beta,
      fixed_nu = sim$true$nu
    )
    spline_metrics <- compute_active_pair_local_metrics(comp, spline_fit)
    cheb_metrics <- compute_active_pair_local_metrics(comp, chebyshev_fit)
    rbind(
      make_metric_row("reduced", rep_id, seed, "spline", reference_fit, spline_metrics, reference_fit$opt$par),
      make_metric_row("reduced", rep_id, seed, "chebyshev", reference_fit, cheb_metrics, reference_fit$opt$par)
    )
  })

  do.call(rbind, rows)
}

evaluate_practical_benchmark <- function(spline_fit, chebyshev_fit) {
  true_spec <- list(beta = c(-0.35, 0.9), delta = 6, b = 0.35, nu = 6)
  seeds <- 20261001L:20261050L
  start <- c(beta_intercept = 0, beta_x1 = 0, delta = 6, b = 0.30)
  lower <- c(-5, -5, 1.0, 0.05)
  upper <- c(5, 5, 20.0, 1.5)

  rows <- lapply(seq_along(seeds), function(rep_id) {
    seed <- seeds[[rep_id]]
    sim <- simulate_clayton_beta_example(
      n = 120L,
      beta = true_spec$beta,
      delta = true_spec$delta,
      b = true_spec$b,
      nu = true_spec$nu,
      seed = seed
    )
    data <- prepare_reduced_clayton_data(sim$y, sim$coords, sim$X, m = 2L)
    reference_fit <- run_opt_with_time(
      objective = wpl_reference_practical,
      start = start,
      lower = lower,
      upper = upper,
      data = data,
      fixed_nu = sim$true$nu
    )

    if (!is_reference_converged(reference_fit$opt)) {
      return(rbind(
        make_missing_row("practical", rep_id, seed, "spline", reference_fit),
        make_missing_row("practical", rep_id, seed, "chebyshev", reference_fit)
      ))
    }

    comp <- build_practical_clayton_components(
      reference_fit$opt$par,
      data = data,
      fixed_nu = sim$true$nu
    )
    spline_metrics <- compute_active_pair_local_metrics(comp, spline_fit)
    cheb_metrics <- compute_active_pair_local_metrics(comp, chebyshev_fit)
    rbind(
      make_metric_row("practical", rep_id, seed, "spline", reference_fit, spline_metrics, reference_fit$opt$par),
      make_metric_row("practical", rep_id, seed, "chebyshev", reference_fit, cheb_metrics, reference_fit$opt$par)
    )
  })

  do.call(rbind, rows)
}

reduced_results <- evaluate_reduced_benchmark(spline_fit, cheb_bundle$fit)
practical_results <- evaluate_practical_benchmark(spline_fit, cheb_bundle$fit)

replicate_results <- rbind(reduced_results, practical_results)
summary_results <- summarize_weighted_error_results(replicate_results)

write.csv(
  replicate_results,
  result_path("example1_active_pair_weighted_error_replicates.csv"),
  row.names = FALSE
)
write.csv(
  summary_results,
  result_path("example1_active_pair_weighted_error_summary.csv"),
  row.names = FALSE
)

