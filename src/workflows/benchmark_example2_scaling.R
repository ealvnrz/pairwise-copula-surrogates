script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

alpha_true <- 1.2
b_true <- 0.30
R_repl <- 60L
m_neighbors <- 6L
n_grid <- c(50L, 100L, 200L, 400L)
n_replicates <- 10L

factor_bundle <- readRDS(result_path("25_factorcopula_example_box_and_refinement.rds"))
factor_spec <- factor_bundle$factor_spec
spline_fit <- readRDS(result_path("factorcopula_example_spline_fit.rds"))
cheb_bundle <- readRDS(result_path("factorcopula_example_chebyshev_fit.rds"))
practical_summary <- read.csv(result_path("factorcopula_example_practical_summary.csv"), stringsAsFactors = FALSE)
replicate_rows <- read.csv(result_path("factorcopula_example_replicates.csv"), stringsAsFactors = FALSE)

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

rows <- lapply(n_grid, function(n_val) {
  do.call(
    rbind,
    lapply(seq_len(n_replicates), function(rep_id) {
      sim <- simulate_spatial_factorcopula_example(
        n = n_val,
        R = R_repl,
        alpha = alpha_true,
        b = b_true,
        factor_spec = factor_spec,
        seed = 20268500L + 100L * match(n_val, n_grid) + rep_id
      )
      data <- prepare_pairwise_factor_data(sim$Y, sim$coords, m = m_neighbors, fixed_factor = factor_spec)
      start <- c(alpha = alpha_true, b = b_true)
      lower <- c(alpha = 0.25, b = 0.05)
      upper <- c(alpha = 2.0, b = 1.5)

      reference_fit <- run_exact_with_time(data, start = start, lower = lower, upper = upper)
      spline_fit_opt <- run_surrogate_with_time(data, fit = spline_fit, start = start, lower = lower, upper = upper)
      chebyshev_fit_opt <- run_surrogate_with_time(data, fit = cheb_bundle$fit, start = start, lower = lower, upper = upper)

      data.frame(
        n = n_val,
        replicate = rep_id,
        method = c("reference", "spline", "chebyshev"),
        mean_time = c(reference_fit$elapsed, spline_fit_opt$elapsed, chebyshev_fit_opt$elapsed),
        fncount = c(reference_fit$fncount, spline_fit_opt$fncount, chebyshev_fit_opt$fncount),
        grcount = c(reference_fit$grcount, spline_fit_opt$grcount, chebyshev_fit_opt$grcount),
        convergence_code = c(reference_fit$convergence, spline_fit_opt$convergence, chebyshev_fit_opt$convergence),
        converged = c(
          identical(reference_fit$opt$convergence, 0L),
          identical(spline_fit_opt$opt$convergence, 0L),
          identical(chebyshev_fit_opt$opt$convergence, 0L)
        ),
        stringsAsFactors = FALSE
      )
    })
  )
})

scaling_raw <- do.call(rbind, rows)
scaling_summary <- do.call(
  rbind,
  lapply(split(scaling_raw, list(scaling_raw$n, scaling_raw$method), drop = TRUE), function(sub_df) {
    converged_df <- sub_df[sub_df$converged, , drop = FALSE]
    data.frame(
      n = sub_df$n[1],
      method = sub_df$method[1],
      n_replicates = nrow(sub_df),
      n_converged = sum(sub_df$converged),
      mean_time = mean(sub_df$mean_time),
      sd_time = stats::sd(sub_df$mean_time),
      mean_time_converged = if (nrow(converged_df) > 0L) mean(converged_df$mean_time) else NA_real_,
      sd_time_converged = if (nrow(converged_df) > 1L) stats::sd(converged_df$mean_time) else NA_real_,
      convergence_rate = mean(sub_df$converged),
      mean_fncount_converged = if (nrow(converged_df) > 0L) mean(converged_df$fncount, na.rm = TRUE) else NA_real_,
      mean_grcount_converged = if (nrow(converged_df) > 0L) mean(converged_df$grcount, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
)
rownames(scaling_summary) <- NULL
scaling_summary <- runtime_scaling_summary(split(scaling_summary, seq_len(nrow(scaling_summary))))

matched_rows <- do.call(
  rbind,
  lapply(n_grid, function(n_val) {
    sub_df <- scaling_raw[scaling_raw$n == n_val, , drop = FALSE]
    ref_df <- sub_df[sub_df$method == "reference", c("replicate", "mean_time", "converged")]
    names(ref_df) <- c("replicate", "reference_time", "reference_converged")
    do.call(
      rbind,
      lapply(c("spline", "chebyshev"), function(method_name) {
        method_df <- sub_df[sub_df$method == method_name, c("replicate", "mean_time", "converged")]
        names(method_df) <- c("replicate", "method_time", "method_converged")
        matched <- merge(method_df, ref_df, by = "replicate", all = FALSE, sort = FALSE)
        matched <- matched[matched$reference_converged & matched$method_converged, , drop = FALSE]
        data.frame(
          n = n_val,
          method = method_name,
          matched_n = nrow(matched),
          matched_reference_mean_time = if (nrow(matched) > 0L) mean(matched$reference_time) else NA_real_,
          matched_mean_time = if (nrow(matched) > 0L) mean(matched$method_time) else NA_real_,
          matched_speedup_vs_reference = if (nrow(matched) > 0L) mean(matched$reference_time) / mean(matched$method_time) else NA_real_,
          stringsAsFactors = FALSE
        )
      })
    )
  })
)

scaling_summary <- merge(
  scaling_summary,
  matched_rows,
  by = c("n", "method"),
  all.x = TRUE,
  sort = FALSE
)
scaling_summary$matched_n[scaling_summary$method == "reference"] <- scaling_summary$n_converged[scaling_summary$method == "reference"]
scaling_summary$matched_reference_mean_time[scaling_summary$method == "reference"] <- scaling_summary$mean_time_converged[scaling_summary$method == "reference"]
scaling_summary$matched_mean_time[scaling_summary$method == "reference"] <- scaling_summary$mean_time_converged[scaling_summary$method == "reference"]
scaling_summary$matched_speedup_vs_reference[scaling_summary$method == "reference"] <- 1
scaling_summary <- scaling_summary[order(scaling_summary$n, match(scaling_summary$method, c("reference", "spline", "chebyshev"))), , drop = FALSE]

write.csv(scaling_raw, result_path("factorcopula_example_runtime_scaling_replicates.csv"), row.names = FALSE)
write.csv(scaling_summary, result_path("factorcopula_example_runtime_scaling.csv"), row.names = FALSE)

save_runtime_scaling_plot(
  scaling_summary,
  result_path("fig_factorcopula_example_runtime_scaling.png"),
  raw_df = scaling_raw
)

