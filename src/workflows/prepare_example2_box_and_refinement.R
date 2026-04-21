script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

alpha_true <- 1.2
b_true <- 0.30
R_repl <- 60L
m_neighbors <- 6L
n_pilot_box <- 100L
n_pilot_obs <- 100L
n_pilot_refinement <- 8L
factor_spec <- default_factorcopula_spec(
  pareto_shape = 2.5,
  pareto_scale = 1.0,
  quadrature_n = 40L,
  lookup_grid_size = 4000L,
  lower_tail_target = 1e-6,
  upper_tail_target = 1e-6
)

run_opt_with_time <- function(objective, start, lower, upper, ...) {
  elapsed <- system.time({
    opt <- optimize_reduced_wpl(objective, start = start, lower = lower, upper = upper, ...)
  })[["elapsed"]]
  list(opt = opt, elapsed = elapsed)
}

run_surrogate_with_time <- function(data, fit, start, lower, upper) {
  elapsed <- system.time({
    cache <- prepare_factorcopula_surrogate_cache(data, fit)
    opt <- optimize_reduced_wpl(
      objective = wpl_surrogate_factorcopula,
      start = start,
      lower = lower,
      upper = upper,
      data = data,
      surrogate_fit = cache,
      fixed_factor = factor_spec
    )
  })[["elapsed"]]
  list(opt = opt, elapsed = elapsed)
}

derive_factorcopula_box <- function() {
  x_vals <- vector("list", n_pilot_box)
  q_vals <- vector("list", n_pilot_box)
  t_vals <- vector("list", n_pilot_box)

  for (rep_id in seq_len(n_pilot_box)) {
    sim <- simulate_spatial_factorcopula_example(
      n = n_pilot_obs,
      R = R_repl,
      alpha = alpha_true,
      b = b_true,
      factor_spec = factor_spec,
      seed = 20268000L + rep_id
    )
    design <- build_unique_nn_pair_design(sim$coords, m = m_neighbors)
    idx <- design$pair_index
    rho <- spatial_factor_correlation(design$distances, alpha = alpha_true, b = b_true)

    transformed <- transform_factorcopula_inputs(
      u = as.vector(sim$u_exact[, idx[, "i"], drop = FALSE]),
      v = as.vector(sim$u_exact[, idx[, "j"], drop = FALSE]),
      rho = rep(rho, each = R_repl),
      fixed_factor = factor_spec
    )
    x_vals[[rep_id]] <- transformed$x
    q_vals[[rep_id]] <- transformed$q
    t_vals[[rep_id]] <- transformed$t
  }

  inflate_limits <- function(values, lower_floor = -Inf) {
    probs <- stats::quantile(unlist(values, use.names = FALSE), probs = c(0.005, 0.995), names = FALSE, type = 8)
    span <- diff(probs)
    pad <- 0.05 * max(span, 1e-8)
    c(max(lower_floor, probs[1] - pad), probs[2] + pad)
  }

  list(
    x = inflate_limits(x_vals),
    q = inflate_limits(q_vals, lower_floor = 0),
    t = inflate_limits(t_vals, lower_floor = 0)
  )
}

fit_candidate_surrogate <- function(method, level, box) {
  if (identical(method, "spline")) {
    grid_size <- switch(as.character(level), `15` = 15L, `21` = 21L, `27` = 27L)
    knot_df <- switch(as.character(level), `15` = 6L, `21` = 8L, `27` = 10L)
    local_table <- build_factorcopula_local_table(
      box = box,
      grid = list(x = grid_size, q = grid_size, t = grid_size),
      fixed_factor = factor_spec
    )
    fit_time <- system.time({
      fit <- fit_spline_surrogate(
        local_table = local_table,
        box = box,
        knots = list(x = knot_df, q = knot_df, t = knot_df, degree = 3L),
        penalty = list(x = 1e-6, q = 1e-6, t = 1e-6),
        nu = 0
      )
    })[["elapsed"]]
    return(list(
      method = method,
      label = sprintf("spline_%d^3", grid_size),
      grid_label = sprintf("%d^3", grid_size),
      fit = fit,
      fit_seconds = fit_time
    ))
  }

  degree <- switch(as.character(level), `10` = 10L, `14` = 14L, `18` = 18L)
  cheb_table <- expand.grid(
    x = chebyshev_nodes(degree, box$x),
    q = chebyshev_nodes(degree, box$q),
    t = chebyshev_nodes(degree, box$t),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  cheb_table$value <- reference_transformed_factorcopula(
    cheb_table$x,
    cheb_table$q,
    cheb_table$t,
    fixed_factor = factor_spec,
    log = TRUE
  )
  fit_time <- system.time({
    fit <- fit_chebyshev_surrogate(
      local_table = cheb_table,
      box = box,
      degree = list(x = degree, q = degree, t = degree),
      nu = 0
    )
  })[["elapsed"]]
  list(
    method = method,
    label = sprintf("chebyshev_%d^3", degree),
    grid_label = sprintf("%d^3", degree),
    fit = fit,
    fit_seconds = fit_time
  )
}

evaluate_candidate <- function(candidate, evaluation_grid, pilot_datasets, reference_opts, start, lower, upper) {
  reference_fun <- function(z) {
    reference_transformed_factorcopula(z[1], z[2], z[3], fixed_factor = factor_spec, log = TRUE)
  }
  local_errors <- surrogate_local_error_summary_generic(
    fit = candidate$fit,
    evaluation_grid = evaluation_grid,
    reference_fun = reference_fun
  )

  replicate_rows <- lapply(seq_along(pilot_datasets), function(rep_id) {
    data_obj <- pilot_datasets[[rep_id]]
    ref_opt <- reference_opts[[rep_id]]
    sur_opt <- run_surrogate_with_time(
      data = data_obj$data,
      fit = candidate$fit,
      start = start,
      lower = lower,
      upper = upper
    )
    ref_par <- ref_opt$opt$par
    sur_par <- sur_opt$opt$par
    data.frame(
      replicate = rep_id,
      method = candidate$method,
      label = candidate$label,
      alpha_gap_vs_reference = sur_par[1] - ref_par[1],
      b_gap_vs_reference = sur_par[2] - ref_par[2],
      alpha_error_vs_truth = sur_par[1] - alpha_true,
      b_error_vs_truth = sur_par[2] - b_true,
      criterion_gap_exact = wpl_reference_factorcopula(ref_par, data_obj$data, fixed_factor = factor_spec) -
        wpl_reference_factorcopula(sur_par, data_obj$data, fixed_factor = factor_spec),
      optimization_time = sur_opt$elapsed,
      converged = identical(sur_opt$opt$convergence, 0L),
      stringsAsFactors = FALSE
    )
  })
  replicate_df <- do.call(rbind, replicate_rows)

  candidate$summary <- data.frame(
    method = candidate$method,
    label = candidate$label,
    grid_label = candidate$grid_label,
    fit_seconds = candidate$fit_seconds,
    mean_local_value_error = local_errors$summary$mean_value_error,
    mean_local_grad_error = local_errors$summary$mean_grad_error,
    mean_local_hess_error = local_errors$summary$mean_hess_error,
    mean_abs_alpha_gap = mean(abs(replicate_df$alpha_gap_vs_reference)),
    mean_abs_b_gap = mean(abs(replicate_df$b_gap_vs_reference)),
    mean_parameter_gap = mean(c(mean(abs(replicate_df$alpha_gap_vs_reference)), mean(abs(replicate_df$b_gap_vs_reference)))),
    mean_criterion_gap_exact = mean(abs(replicate_df$criterion_gap_exact)),
    mean_optimization_time = mean(replicate_df$optimization_time),
    convergence_rate = mean(replicate_df$converged),
    total_cost_proxy = candidate$fit_seconds + mean(replicate_df$optimization_time),
    stringsAsFactors = FALSE
  )
  candidate$replicate_df <- replicate_df
  candidate
}

select_candidate <- function(summary_df) {
  best_local <- min(summary_df$mean_local_value_error)
  best_gap <- min(summary_df$mean_parameter_gap)
  eligible <- summary_df[
    summary_df$mean_local_value_error <= 1.10 * best_local &
      summary_df$mean_parameter_gap <= 1.10 * best_gap,
    ,
    drop = FALSE
  ]
  if (!nrow(eligible)) {
    eligible <- summary_df
  }
  eligible[order(eligible$total_cost_proxy, eligible$fit_seconds), , drop = FALSE][1, , drop = FALSE]
}

box <- derive_factorcopula_box()
box_df <- data.frame(
  coordinate = c("x", "q", "t"),
  lower = c(box$x[1], box$q[1], box$t[1]),
  upper = c(box$x[2], box$q[2], box$t[2]),
  stringsAsFactors = FALSE
)
write.csv(box_df, result_path("factorcopula_example_box.csv"), row.names = FALSE)
saveRDS(box, result_path("factorcopula_example_box.rds"))

evaluation_grid <- sample_transformed_points_by_region(
  box = box,
  n_per_region = 160L,
  seed = 20268090L
)

pilot_sims <- lapply(seq_len(n_pilot_refinement), function(rep_id) {
  sim <- simulate_spatial_factorcopula_example(
    n = 200L,
    R = R_repl,
    alpha = alpha_true,
    b = b_true,
    factor_spec = factor_spec,
    seed = 20268100L + rep_id
  )
  list(
    sim = sim,
    data = prepare_pairwise_factor_data(sim$Y, sim$coords, m = m_neighbors, fixed_factor = factor_spec)
  )
})

start <- c(alpha = alpha_true, b = b_true)
lower <- c(alpha = 0.25, b = 0.05)
upper <- c(alpha = 2.0, b = 1.5)

reference_opts <- lapply(pilot_sims, function(obj) {
  run_opt_with_time(
    objective = wpl_reference_factorcopula,
    start = start,
    lower = lower,
    upper = upper,
    data = obj$data,
    fixed_factor = factor_spec
  )
})

candidate_specs <- c(
  lapply(c(15L, 21L, 27L), function(level) fit_candidate_surrogate("spline", level, box)),
  lapply(c(10L, 14L, 18L), function(level) fit_candidate_surrogate("chebyshev", level, box))
)

evaluated_candidates <- lapply(
  candidate_specs,
  evaluate_candidate,
  evaluation_grid = evaluation_grid,
  pilot_datasets = pilot_sims,
  reference_opts = reference_opts,
  start = start,
  lower = lower,
  upper = upper
)

refinement_summary <- do.call(rbind, lapply(evaluated_candidates, `[[`, "summary"))
selected_rows <- do.call(
  rbind,
  lapply(split(refinement_summary, refinement_summary$method), select_candidate)
)
refinement_summary$selected <- refinement_summary$label %in% selected_rows$label
write.csv(refinement_summary, result_path("factorcopula_example_refinement_summary.csv"), row.names = FALSE)

selected_summary <- refinement_summary[refinement_summary$selected, , drop = FALSE]
selected_candidates <- lapply(split(selected_summary, selected_summary$method), function(row_df) {
  label <- row_df$label[1]
  evaluated_candidates[[match(label, vapply(evaluated_candidates, `[[`, character(1), "label"))]]
})

saveRDS(
  list(
    box = box,
    refinement_summary = refinement_summary,
    selected = selected_rows,
    factor_spec = factor_spec
  ),
  result_path("25_factorcopula_example_box_and_refinement.rds")
)
saveRDS(selected_candidates[["spline"]]$fit, result_path("factorcopula_example_spline_fit.rds"))
saveRDS(list(fit = selected_candidates[["chebyshev"]]$fit), result_path("factorcopula_example_chebyshev_fit.rds"))

