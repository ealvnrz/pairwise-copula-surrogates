script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

alpha_true <- 1.2
b_true <- 0.30
R_repl <- 60L
m_neighbors <- 6L
n_obs <- 200L
n_replicates <- 50L

factor_bundle <- readRDS(result_path("25_factorcopula_example_box_and_refinement.rds"))
factor_spec <- factor_bundle$factor_spec
spline_fit <- readRDS(result_path("factorcopula_example_spline_fit.rds"))
cheb_bundle <- readRDS(result_path("factorcopula_example_chebyshev_fit.rds"))

run_exact_with_time <- function(data, start, lower, upper) {
  run_opt_with_timing(
    objective = wpl_reference_factorcopula,
    start = start,
    lower = lower,
    upper = upper,
    data = data,
    fixed_factor = factor_spec
  )
}

run_surrogate_with_time <- function(data, fit, start, lower, upper) {
  cache <- prepare_factorcopula_surrogate_cache(data, fit)
  run_opt_with_timing(
    objective = wpl_surrogate_factorcopula,
    start = start,
    lower = lower,
    upper = upper,
    data = data,
    surrogate_fit = cache,
    fixed_factor = factor_spec
  )
}

collect_method_row <- function(method, par, elapsed, converged, data_obj) {
  reference_theta <- data_obj$reference_fit$opt$par
  criterion_at_theta <- wpl_reference_factorcopula(par, data_obj$data, fixed_factor = factor_spec)
  criterion_at_ref <- wpl_reference_factorcopula(reference_theta, data_obj$data, fixed_factor = factor_spec)
  data.frame(
    method = method,
    alpha = par[1],
    b = par[2],
    alpha_gap_vs_reference = par[1] - reference_theta[1],
    b_gap_vs_reference = par[2] - reference_theta[2],
    alpha_error_vs_truth = par[1] - alpha_true,
    b_error_vs_truth = par[2] - b_true,
    criterion_gap_exact = criterion_at_ref - criterion_at_theta,
    criterion_gap_exact_pct = 100 * abs(criterion_at_ref - criterion_at_theta) / pmax(abs(criterion_at_ref), 1e-8),
    optimization_time = elapsed,
    converged = converged,
    fncount = extract_optim_count(data_obj[[paste0(method, "_fit")]]$opt, "function"),
    grcount = extract_optim_count(data_obj[[paste0(method, "_fit")]]$opt, "gradient"),
    stringsAsFactors = FALSE
  )
}

fit_dataset <- function(seed) {
  sim <- simulate_spatial_factorcopula_example(
    n = n_obs,
    R = R_repl,
    alpha = alpha_true,
    b = b_true,
    factor_spec = factor_spec,
    seed = seed
  )
  data <- prepare_pairwise_factor_data(sim$Y, sim$coords, m = m_neighbors, fixed_factor = factor_spec)

  start <- c(alpha = alpha_true, b = b_true)
  lower <- c(alpha = 0.25, b = 0.05)
  upper <- c(alpha = 2.0, b = 1.5)

  reference_fit <- run_exact_with_time(data, start = start, lower = lower, upper = upper)
  spline_fit_opt <- run_surrogate_with_time(data, fit = spline_fit, start = start, lower = lower, upper = upper)
  cheb_fit_opt <- run_surrogate_with_time(data, fit = cheb_bundle$fit, start = start, lower = lower, upper = upper)

  list(
    simulation = sim,
    data = data,
    reference_fit = reference_fit,
    spline_fit = spline_fit_opt,
    chebyshev_fit = cheb_fit_opt
  )
}

illustrative <- fit_dataset(20268301L)
illustrative_optima <- rbind(
  collect_method_row("reference", illustrative$reference_fit$opt$par, illustrative$reference_fit$elapsed, identical(illustrative$reference_fit$opt$convergence, 0L), illustrative),
  collect_method_row("spline", illustrative$spline_fit$opt$par, illustrative$spline_fit$elapsed, identical(illustrative$spline_fit$opt$convergence, 0L), illustrative),
  collect_method_row("chebyshev", illustrative$chebyshev_fit$opt$par, illustrative$chebyshev_fit$elapsed, identical(illustrative$chebyshev_fit$opt$convergence, 0L), illustrative)
)

replicate_objs <- lapply(seq_len(n_replicates), function(rep_id) {
  fit_dataset(20268400L + rep_id)
})

replicate_rows <- do.call(
  rbind,
  lapply(seq_along(replicate_objs), function(rep_id) {
    obj <- replicate_objs[[rep_id]]
    rbind(
      transform(collect_method_row("reference", obj$reference_fit$opt$par, obj$reference_fit$elapsed, identical(obj$reference_fit$opt$convergence, 0L), obj), replicate = rep_id),
      transform(collect_method_row("spline", obj$spline_fit$opt$par, obj$spline_fit$elapsed, identical(obj$spline_fit$opt$convergence, 0L), obj), replicate = rep_id),
      transform(collect_method_row("chebyshev", obj$chebyshev_fit$opt$par, obj$chebyshev_fit$elapsed, identical(obj$chebyshev_fit$opt$convergence, 0L), obj), replicate = rep_id)
    )
  })
)

practical_summary <- do.call(
  rbind,
  lapply(split(replicate_rows, replicate_rows$method), function(df) {
    data.frame(
      method = df$method[1],
      n = nrow(df),
      convergence_rate = mean(df$converged),
      mean_abs_alpha_gap_vs_reference = if (df$method[1] == "reference") NA_real_ else mean(abs(df$alpha_gap_vs_reference)),
      mean_abs_b_gap_vs_reference = if (df$method[1] == "reference") NA_real_ else mean(abs(df$b_gap_vs_reference)),
      mean_abs_alpha_error_vs_truth = mean(abs(df$alpha_error_vs_truth)),
      mean_abs_b_error_vs_truth = mean(abs(df$b_error_vs_truth)),
      mean_criterion_gap_exact = if (df$method[1] == "reference") 0 else mean(abs(df$criterion_gap_exact)),
      mean_criterion_gap_exact_pct = if (df$method[1] == "reference") 0 else mean(abs(df$criterion_gap_exact_pct)),
      mean_time = mean(df$optimization_time),
      sd_time = stats::sd(df$optimization_time),
      mean_fncount = mean(df$fncount, na.rm = TRUE),
      mean_grcount = mean(df$grcount, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
rownames(practical_summary) <- NULL
practical_summary <- practical_summary[match(c("reference", "spline", "chebyshev"), practical_summary$method), , drop = FALSE]

write.csv(illustrative_optima, result_path("factorcopula_example_illustrative_optima.csv"), row.names = FALSE)
write.csv(replicate_rows, result_path("factorcopula_example_replicates.csv"), row.names = FALSE)
write.csv(practical_summary, result_path("factorcopula_example_practical_summary.csv"), row.names = FALSE)

