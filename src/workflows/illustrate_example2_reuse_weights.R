script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

alpha_true <- 1.2
b_true <- 0.30
R_repl <- 60L
m_neighbors_primary <- 6L
m_neighbors_fallback <- 4L
n_obs <- 200L
n_replicates <- 5L
seed_block <- 20268501L + seq_len(n_replicates) - 1L

factor_bundle <- readRDS(result_path("25_factorcopula_example_box_and_refinement.rds"))
factor_spec <- factor_bundle$factor_spec
spline_fit <- readRDS(result_path("factorcopula_example_spline_fit.rds"))
cheb_bundle <- readRDS(result_path("factorcopula_example_chebyshev_fit.rds"))

start <- c(alpha = alpha_true, b = b_true)
lower <- c(alpha = 0.25, b = 0.05)
upper <- c(alpha = 2.0, b = 1.5)

normalized_inverse_distance_weights <- function(design) {
  raw <- 1 / pmax(as.numeric(design$distances), .Machine$double.eps)
  raw / mean(raw)
}

run_exact_with_time <- function(data, weights, wrapper_label) {
  elapsed <- system.time({
    opt <- optimize_reduced_wpl(
      objective = wpl_reference_factorcopula,
      start = start,
      lower = lower,
      upper = upper,
      data = data,
      fixed_factor = factor_spec,
      weights = weights
    )
  })[["elapsed"]]
  list(opt = opt, elapsed = elapsed, wrapper = wrapper_label)
}

run_surrogate_with_time <- function(data, fit, weights, wrapper_label) {
  cache <- prepare_factorcopula_surrogate_cache(data, fit)
  elapsed <- system.time({
    opt <- optimize_reduced_wpl(
      objective = wpl_surrogate_factorcopula,
      start = start,
      lower = lower,
      upper = upper,
      data = data,
      surrogate_fit = cache,
      fixed_factor = factor_spec,
      weights = weights
    )
  })[["elapsed"]]
  list(opt = opt, elapsed = elapsed, cache = cache, wrapper = wrapper_label)
}

collect_method_row <- function(method, par, elapsed, converged, data_obj, weights) {
  reference_theta <- data_obj$reference_fit$opt$par
  criterion_at_theta <- wpl_reference_factorcopula(par, data_obj$data, fixed_factor = factor_spec, weights = weights)
  criterion_at_ref <- wpl_reference_factorcopula(reference_theta, data_obj$data, fixed_factor = factor_spec, weights = weights)
  data.frame(
    method = method,
    alpha = par[1],
    b = par[2],
    alpha_gap_vs_reference = par[1] - reference_theta[1],
    b_gap_vs_reference = par[2] - reference_theta[2],
    criterion_gap_exact_pct = 100 * abs(criterion_at_ref - criterion_at_theta) / pmax(abs(criterion_at_ref), 1e-8),
    optimization_time = elapsed,
    converged = converged,
    stringsAsFactors = FALSE
  )
}

fit_dataset <- function(seed, m_neighbors, use_inverse_distance_weights) {
  sim <- simulate_spatial_factorcopula_example(
    n = n_obs,
    R = R_repl,
    alpha = alpha_true,
    b = b_true,
    factor_spec = factor_spec,
    seed = seed
  )
  data <- prepare_pairwise_factor_data(sim$Y, sim$coords, m = m_neighbors, fixed_factor = factor_spec)
  weights <- if (isTRUE(use_inverse_distance_weights)) {
    normalized_inverse_distance_weights(data$design)
  } else {
    rep(1, nrow(data$design$pair_index))
  }

  list(
    simulation = sim,
    data = data,
    weights = weights,
    weight_ratio = max(weights) / min(weights),
    reference_fit = run_exact_with_time(
      data = data,
      weights = weights,
      wrapper_label = if (isTRUE(use_inverse_distance_weights)) "normalized_inverse_distance" else "unit_weights_m4"
    ),
    spline_fit = run_surrogate_with_time(
      data = data,
      fit = spline_fit,
      weights = weights,
      wrapper_label = if (isTRUE(use_inverse_distance_weights)) "normalized_inverse_distance" else "unit_weights_m4"
    ),
    chebyshev_fit = run_surrogate_with_time(
      data = data,
      fit = cheb_bundle$fit,
      weights = weights,
      wrapper_label = if (isTRUE(use_inverse_distance_weights)) "normalized_inverse_distance" else "unit_weights_m4"
    )
  )
}

summarize_rows <- function(replicate_rows) {
  out <- do.call(
    rbind,
    lapply(split(replicate_rows, replicate_rows$method), function(df) {
      data.frame(
        method = df$method[1],
        n = nrow(df),
        convergence_rate = mean(df$converged),
        mean_abs_alpha_gap_vs_reference = if (df$method[1] == "reference") NA_real_ else mean(abs(df$alpha_gap_vs_reference)),
        mean_abs_b_gap_vs_reference = if (df$method[1] == "reference") NA_real_ else mean(abs(df$b_gap_vs_reference)),
        mean_criterion_gap_exact_pct = if (df$method[1] == "reference") 0 else mean(abs(df$criterion_gap_exact_pct)),
        mean_time = mean(df$optimization_time),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out[match(c("reference", "spline", "chebyshev"), out$method), , drop = FALSE]
}

run_configuration <- function(m_neighbors, use_inverse_distance_weights) {
  replicate_objs <- lapply(seed_block, function(seed) {
    fit_dataset(
      seed = seed,
      m_neighbors = m_neighbors,
      use_inverse_distance_weights = use_inverse_distance_weights
    )
  })

  replicate_rows <- do.call(
    rbind,
    lapply(seq_along(replicate_objs), function(rep_id) {
      obj <- replicate_objs[[rep_id]]
      weights <- obj$weights
      rbind(
        transform(
          collect_method_row(
            "reference",
            obj$reference_fit$opt$par,
            obj$reference_fit$elapsed,
            identical(obj$reference_fit$opt$convergence, 0L),
            obj,
            weights
          ),
          replicate = rep_id,
          wrapper = obj$reference_fit$wrapper,
          m_neighbors = m_neighbors,
          weight_ratio = obj$weight_ratio
        ),
        transform(
          collect_method_row(
            "spline",
            obj$spline_fit$opt$par,
            obj$spline_fit$elapsed,
            identical(obj$spline_fit$opt$convergence, 0L),
            obj,
            weights
          ),
          replicate = rep_id,
          wrapper = obj$spline_fit$wrapper,
          m_neighbors = m_neighbors,
          weight_ratio = obj$weight_ratio
        ),
        transform(
          collect_method_row(
            "chebyshev",
            obj$chebyshev_fit$opt$par,
            obj$chebyshev_fit$elapsed,
            identical(obj$chebyshev_fit$opt$convergence, 0L),
            obj,
            weights
          ),
          replicate = rep_id,
          wrapper = obj$chebyshev_fit$wrapper,
          m_neighbors = m_neighbors,
          weight_ratio = obj$weight_ratio
        )
      )
    })
  )

  list(
    replicate_rows = replicate_rows,
    summary = summarize_rows(replicate_rows),
    replicate_objects = replicate_objs
  )
}

primary_run <- run_configuration(
  m_neighbors = m_neighbors_primary,
  use_inverse_distance_weights = TRUE
)

primary_summary <- primary_run$summary
primary_weight_ratios <- vapply(primary_run$replicate_objects, function(obj) obj$weight_ratio, numeric(1))
primary_is_flat <- all(primary_weight_ratios <= 1.10)
primary_is_unstable <- any(primary_summary$convergence_rate < 0.8)
fallback_needed <- isTRUE(primary_is_flat || primary_is_unstable)

final_run <- primary_run
fallback_reason <- "none"

if (fallback_needed) {
  final_run <- run_configuration(
    m_neighbors = m_neighbors_fallback,
    use_inverse_distance_weights = FALSE
  )
  fallback_reason <- paste(
    c(
      if (primary_is_flat) "inverse-distance weights effectively flat on all datasets" else NULL,
      if (primary_is_unstable) "primary weighted benchmark unstable under convergence threshold" else NULL
    ),
    collapse = "; "
  )
}

summary_df <- final_run$summary
replicate_df <- final_run$replicate_rows

write.csv(
  replicate_df,
  result_path("factorcopula_reuse_weights_replicates.csv"),
  row.names = FALSE
)
write.csv(
  summary_df,
  result_path("factorcopula_reuse_weights_summary.csv"),
  row.names = FALSE
)

