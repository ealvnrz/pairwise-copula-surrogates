script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

run_opt_with_time <- function(objective, start, lower, upper, ...) {
  elapsed <- system.time({
    opt <- optimize_reduced_wpl(objective, start = start, lower = lower, upper = upper, ...)
  })[["elapsed"]]
  list(opt = opt, elapsed = elapsed)
}

refit_selected_spline <- function(local_table, box, dims, degree, penalty, nu = 0) {
  fit_spline_surrogate(
    local_table = local_table,
    box = box,
    knots = utils::modifyList(as.list(dims), list(degree = degree)),
    penalty = penalty,
    nu = nu
  )
}

collect_penalty_rows <- function(example, penalty_label, local_mean_abs_error, local_sup_abs_error, par, exact_criterion, elapsed) {
  data.frame(
    example = example,
    penalty = penalty_label,
    local_mean_abs_error = local_mean_abs_error,
    local_sup_abs_error = local_sup_abs_error,
    par1 = unname(par[1]),
    par2 = unname(par[2]),
    exact_criterion_at_optimum = exact_criterion,
    optimization_time = elapsed,
    stringsAsFactors = FALSE
  )
}

example1_reference_bundle <- readRDS(result_path("01_exact_local_grid.rds"))
example1_base_fit <- readRDS(result_path("02_spline_fit.rds"))
example1_eval <- unique(utils::read.csv(result_path("04_spline_vs_reference_pointwise.csv"))[, c("x", "q", "t", "reference")])

example1_penalty_zero <- refit_selected_spline(
  local_table = example1_reference_bundle$primary$local_table,
  box = example1_reference_bundle$primary$box,
  dims = example1_base_fit$dims,
  degree = example1_base_fit$basis[[1]]$degree,
  penalty = list(x = 0, q = 0, t = 0),
  nu = example1_base_fit$nu
)
example1_penalty_tiny <- refit_selected_spline(
  local_table = example1_reference_bundle$primary$local_table,
  box = example1_reference_bundle$primary$box,
  dims = example1_base_fit$dims,
  degree = example1_base_fit$basis[[1]]$degree,
  penalty = list(x = 1e-6, q = 1e-6, t = 1e-6),
  nu = example1_base_fit$nu
)

example1_pred_zero <- evaluate_transformed_surrogate(example1_penalty_zero, example1_eval$x, example1_eval$q, example1_eval$t)
example1_pred_tiny <- evaluate_transformed_surrogate(example1_penalty_tiny, example1_eval$x, example1_eval$q, example1_eval$t)

example1_sim <- simulate_clayton_beta_example(
  n = 120L,
  beta = c(-0.35, 0.9),
  delta = 6,
  b = 0.35,
  nu = 6,
  seed = 20260409L
)
example1_data <- prepare_reduced_clayton_data(example1_sim$y, example1_sim$coords, example1_sim$X, m = 2L)
example1_start <- c(beta_intercept = 0, beta_x1 = 0, delta = 6, b = 0.30)
example1_lower <- c(-5, -5, 1.0, 0.05)
example1_upper <- c(5, 5, 20.0, 1.5)

example1_opt_zero <- run_opt_with_time(
  objective = wpl_surrogate_practical,
  start = example1_start,
  lower = example1_lower,
  upper = example1_upper,
  data = example1_data,
  surrogate_fit = example1_penalty_zero,
  fixed_nu = example1_sim$true$nu
)
example1_opt_tiny <- run_opt_with_time(
  objective = wpl_surrogate_practical,
  start = example1_start,
  lower = example1_lower,
  upper = example1_upper,
  data = example1_data,
  surrogate_fit = example1_penalty_tiny,
  fixed_nu = example1_sim$true$nu
)

example1_exact_zero <- wpl_reference_practical(example1_opt_zero$opt$par, example1_data, fixed_nu = example1_sim$true$nu)
example1_exact_tiny <- wpl_reference_practical(example1_opt_tiny$opt$par, example1_data, fixed_nu = example1_sim$true$nu)

factor_bundle <- readRDS(result_path("25_factorcopula_example_box_and_refinement.rds"))
example2_base_fit <- readRDS(result_path("factorcopula_example_spline_fit.rds"))
example2_eval_raw <- utils::read.csv(result_path("factorcopula_example_local_accuracy_pointwise.csv"))
example2_eval <- unique(example2_eval_raw[, c("x", "q", "t", "reference")])
example2_table <- build_factorcopula_local_table(
  box = factor_bundle$box,
  grid = list(x = 27L, q = 27L, t = 27L),
  fixed_factor = factor_bundle$factor_spec
)

example2_penalty_zero <- refit_selected_spline(
  local_table = example2_table,
  box = factor_bundle$box,
  dims = example2_base_fit$dims,
  degree = example2_base_fit$basis[[1]]$degree,
  penalty = list(x = 0, q = 0, t = 0),
  nu = example2_base_fit$nu
)
example2_penalty_tiny <- refit_selected_spline(
  local_table = example2_table,
  box = factor_bundle$box,
  dims = example2_base_fit$dims,
  degree = example2_base_fit$basis[[1]]$degree,
  penalty = list(x = 1e-6, q = 1e-6, t = 1e-6),
  nu = example2_base_fit$nu
)

example2_pred_zero <- evaluate_transformed_surrogate(example2_penalty_zero, example2_eval$x, example2_eval$q, example2_eval$t)
example2_pred_tiny <- evaluate_transformed_surrogate(example2_penalty_tiny, example2_eval$x, example2_eval$q, example2_eval$t)

example2_sim <- simulate_spatial_factorcopula_example(
  n = 200L,
  R = 60L,
  alpha = 1.2,
  b = 0.30,
  factor_spec = factor_bundle$factor_spec,
  seed = 20268301L
)
example2_data <- prepare_pairwise_factor_data(
  example2_sim$Y,
  example2_sim$coords,
  m = 6L,
  fixed_factor = factor_bundle$factor_spec
)
example2_start <- c(alpha = 1.2, b = 0.30)
example2_lower <- c(alpha = 0.25, b = 0.05)
example2_upper <- c(alpha = 2.0, b = 1.5)

run_factor_with_cache <- function(fit) {
  elapsed <- system.time({
    cache <- prepare_factorcopula_surrogate_cache(example2_data, fit)
    opt <- optimize_reduced_wpl(
      objective = wpl_surrogate_factorcopula,
      start = example2_start,
      lower = example2_lower,
      upper = example2_upper,
      data = example2_data,
      surrogate_fit = cache,
      fixed_factor = factor_bundle$factor_spec
    )
  })[["elapsed"]]
  list(opt = opt, elapsed = elapsed)
}

example2_opt_zero <- run_factor_with_cache(example2_penalty_zero)
example2_opt_tiny <- run_factor_with_cache(example2_penalty_tiny)

example2_exact_zero <- wpl_reference_factorcopula(example2_opt_zero$opt$par, example2_data, fixed_factor = factor_bundle$factor_spec)
example2_exact_tiny <- wpl_reference_factorcopula(example2_opt_tiny$opt$par, example2_data, fixed_factor = factor_bundle$factor_spec)

by_penalty <- do.call(
  rbind,
  list(
    collect_penalty_rows(
      "example1_clayton",
      "zero",
      mean(abs(example1_pred_zero - example1_eval$reference)),
      max(abs(example1_pred_zero - example1_eval$reference)),
      example1_opt_zero$opt$par,
      example1_exact_zero,
      example1_opt_zero$elapsed
    ),
    collect_penalty_rows(
      "example1_clayton",
      "tiny",
      mean(abs(example1_pred_tiny - example1_eval$reference)),
      max(abs(example1_pred_tiny - example1_eval$reference)),
      example1_opt_tiny$opt$par,
      example1_exact_tiny,
      example1_opt_tiny$elapsed
    ),
    collect_penalty_rows(
      "example2_factorcopula",
      "zero",
      mean(abs(example2_pred_zero - example2_eval$reference)),
      max(abs(example2_pred_zero - example2_eval$reference)),
      example2_opt_zero$opt$par,
      example2_exact_zero,
      example2_opt_zero$elapsed
    ),
    collect_penalty_rows(
      "example2_factorcopula",
      "tiny",
      mean(abs(example2_pred_tiny - example2_eval$reference)),
      max(abs(example2_pred_tiny - example2_eval$reference)),
      example2_opt_tiny$opt$par,
      example2_exact_tiny,
      example2_opt_tiny$elapsed
    )
  )
)

summary_rows <- rbind(
  data.frame(
    example = "example1_clayton",
    mean_abs_local_prediction_diff = mean(abs(example1_pred_zero - example1_pred_tiny)),
    sup_abs_local_prediction_diff = max(abs(example1_pred_zero - example1_pred_tiny)),
    abs_par1_diff = abs(unname(example1_opt_zero$opt$par[1] - example1_opt_tiny$opt$par[1])),
    abs_par2_diff = abs(unname(example1_opt_zero$opt$par[2] - example1_opt_tiny$opt$par[2])),
    abs_par3_diff = abs(unname(example1_opt_zero$opt$par[3] - example1_opt_tiny$opt$par[3])),
    abs_par4_diff = abs(unname(example1_opt_zero$opt$par[4] - example1_opt_tiny$opt$par[4])),
    abs_exact_criterion_diff = abs(example1_exact_zero - example1_exact_tiny),
    rel_exact_criterion_diff_pct = 100 * abs(example1_exact_zero - example1_exact_tiny) / max(abs(example1_exact_tiny), 1e-8),
    stringsAsFactors = FALSE
  ),
  data.frame(
    example = "example2_factorcopula",
    mean_abs_local_prediction_diff = mean(abs(example2_pred_zero - example2_pred_tiny)),
    sup_abs_local_prediction_diff = max(abs(example2_pred_zero - example2_pred_tiny)),
    abs_par1_diff = abs(unname(example2_opt_zero$opt$par[1] - example2_opt_tiny$opt$par[1])),
    abs_par2_diff = abs(unname(example2_opt_zero$opt$par[2] - example2_opt_tiny$opt$par[2])),
    abs_par3_diff = NA_real_,
    abs_par4_diff = NA_real_,
    abs_exact_criterion_diff = abs(example2_exact_zero - example2_exact_tiny),
    rel_exact_criterion_diff_pct = 100 * abs(example2_exact_zero - example2_exact_tiny) / max(abs(example2_exact_tiny), 1e-8),
    stringsAsFactors = FALSE
  )
)

utils::write.csv(by_penalty, result_path("spline_penalty_sensitivity_by_penalty.csv"), row.names = FALSE)
utils::write.csv(summary_rows, result_path("spline_penalty_sensitivity_summary.csv"), row.names = FALSE)

