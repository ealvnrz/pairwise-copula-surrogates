script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

spline_fit <- readRDS(result_path("02_spline_fit.rds"))
cheb_bundle <- readRDS(result_path("03_chebyshev_fit.rds"))

sim <- simulate_clayton_beta_example(
  n = 120L,
  beta = c(-0.35, 0.9),
  delta = 6,
  b = 0.35,
  nu = 6,
  seed = 20260409
)
data <- prepare_reduced_clayton_data(sim$y, sim$coords, sim$X, m = 2L)

start <- c(beta_intercept = 0, beta_x1 = 0, delta = 6, b = 0.30)
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

illustrative_optima <- summarize_practical_optima(
  reference_opt = ref_fit$opt,
  spline_opt = spline_fit_opt$opt,
  chebyshev_opt = cheb_fit_opt$opt,
  data = data,
  fixed_nu = sim$true$nu
)
illustrative_optima$optimization_time <- c(ref_fit$elapsed, spline_fit_opt$elapsed, cheb_fit_opt$elapsed)
illustrative_optima$fncount <- c(ref_fit$fncount, spline_fit_opt$fncount, cheb_fit_opt$fncount)
illustrative_optima$grcount <- c(ref_fit$grcount, spline_fit_opt$grcount, cheb_fit_opt$grcount)

replicates <- lapply(seq_len(50L), function(rep_id) {
  rep_sim <- simulate_clayton_beta_example(
    n = 120L,
    beta = sim$true$beta,
    delta = sim$true$delta,
    b = sim$true$b,
    nu = sim$true$nu,
    seed = 20261000 + rep_id
  )
  rep_data <- prepare_reduced_clayton_data(rep_sim$y, rep_sim$coords, rep_sim$X, m = 2L)
  list(
    simulation = rep_sim,
    data = rep_data,
    reference_opt = run_opt_with_timing(
      objective = wpl_reference_practical,
      start = start,
      lower = lower,
      upper = upper,
      data = rep_data,
      fixed_nu = rep_sim$true$nu
    ),
    spline_opt = run_opt_with_timing(
      objective = wpl_surrogate_practical,
      start = start,
      lower = lower,
      upper = upper,
      data = rep_data,
      surrogate_fit = spline_fit,
      fixed_nu = rep_sim$true$nu
    ),
    chebyshev_opt = run_opt_with_timing(
      objective = wpl_surrogate_practical,
      start = start,
      lower = lower,
      upper = upper,
      data = rep_data,
      surrogate_fit = cheb_bundle$fit,
      fixed_nu = rep_sim$true$nu
    )
  )
})

replicate_gap_df <- summarize_practical_replicate_gaps(
  lapply(replicates, function(obj) {
    list(
      data = obj$data,
      reference_opt = obj$reference_opt$opt,
      spline_opt = obj$spline_opt$opt,
      chebyshev_opt = obj$chebyshev_opt$opt
    )
  }),
  truth = sim$true,
  parameter_names = practical_parameter_names(data)
)
timing_df <- do.call(
  rbind,
  lapply(seq_along(replicates), function(k) {
    data.frame(
      replicate = k,
      method = c("reference", "spline", "chebyshev"),
      optimization_time = c(
        replicates[[k]]$reference_opt$elapsed,
        replicates[[k]]$spline_opt$elapsed,
        replicates[[k]]$chebyshev_opt$elapsed
      ),
      fncount = c(
        replicates[[k]]$reference_opt$fncount,
        replicates[[k]]$spline_opt$fncount,
        replicates[[k]]$chebyshev_opt$fncount
      ),
      grcount = c(
        replicates[[k]]$reference_opt$grcount,
        replicates[[k]]$spline_opt$grcount,
        replicates[[k]]$chebyshev_opt$grcount
      ),
      stringsAsFactors = FALSE
    )
  })
)
summary_df <- summarize_gap_table(
  replicate_gap_df,
  gap_columns = c(
    "beta_intercept_gap_vs_reference",
    "beta_x1_gap_vs_reference",
    "delta_gap_vs_reference",
    "b_gap_vs_reference"
  ),
  error_columns = c(
    "beta_intercept_error_vs_truth",
    "beta_x1_error_vs_truth",
    "delta_error_vs_truth",
    "b_error_vs_truth"
  )
)
criterion_summary <- do.call(
  rbind,
  lapply(split(replicate_gap_df, replicate_gap_df$method), function(sub_df) {
    data.frame(
      method = sub_df$method[1],
      mean_abs_exact_criterion_gap = mean(abs(sub_df$criterion_gap_exact), na.rm = TRUE),
      mean_abs_exact_criterion_gap_pct = mean(abs(sub_df$criterion_gap_exact_pct), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
rownames(criterion_summary) <- NULL
summary_df <- merge(summary_df, criterion_summary, by = "method", all.x = TRUE, sort = FALSE)
summary_df <- merge(
  summary_df,
  stats::aggregate(optimization_time ~ method, data = timing_df, FUN = mean),
  by = "method",
  all.x = TRUE,
  sort = FALSE
)
summary_df <- merge(
  summary_df,
  stats::aggregate(fncount ~ method, data = timing_df, FUN = mean),
  by = "method",
  all.x = TRUE,
  sort = FALSE
)
summary_df <- merge(
  summary_df,
  stats::aggregate(grcount ~ method, data = timing_df, FUN = mean),
  by = "method",
  all.x = TRUE,
  sort = FALSE
)
summary_df <- summary_df[match(c("reference", "spline", "chebyshev"), summary_df$method), , drop = FALSE]

write.csv(illustrative_optima, result_path("practical_beta_free_optima.csv"), row.names = FALSE)
write.csv(replicate_gap_df, result_path("practical_beta_free_replicates.csv"), row.names = FALSE)
write.csv(summary_df, result_path("practical_beta_free_summary.csv"), row.names = FALSE)

