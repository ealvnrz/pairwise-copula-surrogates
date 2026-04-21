script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

spline_fit <- readRDS(result_path("02_spline_fit.rds"))
cheb_bundle <- readRDS(result_path("03_chebyshev_fit.rds"))

n_grid <- c(50L, 100L, 200L, 400L, 800L)
beta_true <- c(-0.35, 0.9)
n_replicates <- 3L
rows <- lapply(n_grid, function(n_val) {
  replicate_rows <- lapply(seq_len(n_replicates), function(rep_id) {
    sim <- simulate_clayton_beta_example(
      n = n_val,
      beta = beta_true,
      delta = 6,
      b = 0.35,
      nu = 6,
      seed = 20262000 + 100L * match(n_val, n_grid) + rep_id
    )
    data <- prepare_reduced_clayton_data(sim$y, sim$coords, sim$X, m = 2L)
    # A common warm start isolates computational scaling from optimizer transients.
    start <- c(beta_intercept = beta_true[1], beta_x1 = beta_true[2], delta = 6, b = 0.35)
    lower <- c(-5, -5, 1.0, 0.05)
    upper <- c(5, 5, 20.0, 1.5)

    ref_fit <- run_opt_with_timing(
      objective = wpl_reference_practical,
      start = start,
      lower = lower,
      upper = upper,
      data = data,
      fixed_nu = sim$true$nu
    )
    spline_fit_opt <- run_opt_with_timing(
      objective = wpl_surrogate_practical,
      start = start,
      lower = lower,
      upper = upper,
      data = data,
      surrogate_fit = spline_fit,
      fixed_nu = sim$true$nu
    )
    cheb_fit_opt <- run_opt_with_timing(
      objective = wpl_surrogate_practical,
      start = start,
      lower = lower,
      upper = upper,
      data = data,
      surrogate_fit = cheb_bundle$fit,
      fixed_nu = sim$true$nu
    )

    data.frame(
      n = n_val,
      replicate = rep_id,
      method = c("reference", "spline", "chebyshev"),
      mean_time = c(ref_fit$elapsed, spline_fit_opt$elapsed, cheb_fit_opt$elapsed),
      fncount = c(ref_fit$fncount, spline_fit_opt$fncount, cheb_fit_opt$fncount),
      grcount = c(ref_fit$grcount, spline_fit_opt$grcount, cheb_fit_opt$grcount),
      convergence_code = c(ref_fit$convergence, spline_fit_opt$convergence, cheb_fit_opt$convergence),
      converged = c(
        identical(ref_fit$opt$convergence, 0L),
        identical(spline_fit_opt$opt$convergence, 0L),
        identical(cheb_fit_opt$opt$convergence, 0L)
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, replicate_rows)
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

write.csv(scaling_raw, result_path("runtime_scaling_replicates.csv"), row.names = FALSE)
write.csv(scaling_summary, result_path("runtime_scaling_summary.csv"), row.names = FALSE)

save_runtime_scaling_plot(
  scaling_summary,
  result_path("fig_runtime_scaling.png"),
  raw_df = scaling_raw
)

