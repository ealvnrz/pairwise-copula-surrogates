script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

spline_fit <- readRDS(result_path("02_spline_fit.rds"))
cheb_bundle <- readRDS(result_path("03_chebyshev_fit.rds"))

sim <- simulate_clayton_beta_example(
  n = 120L,
  beta = 0,
  delta = 6,
  b = 0.35,
  nu = 6,
  seed = 20260409
)
data <- prepare_reduced_clayton_data(sim$y, sim$coords, sim$X, m = 2L)

delta_grid <- seq(4.5, 7.8, length.out = 51)
b_grid <- seq(0.20, 0.50, length.out = 51)
surface <- expand.grid(delta = delta_grid, b = b_grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
parameter_grid <- split(surface[, c("delta", "b")], seq_len(nrow(surface)))

surface$reference <- vapply(parameter_grid, function(par) {
  wpl_reference_reduced(unlist(par, use.names = FALSE), data = data, fixed_beta = sim$true$beta, fixed_nu = sim$true$nu)
}, numeric(1))
surface$spline <- vapply(parameter_grid, function(par) {
  wpl_surrogate_reduced(unlist(par, use.names = FALSE), data = data, surrogate_fit = spline_fit, fixed_beta = sim$true$beta, fixed_nu = sim$true$nu)
}, numeric(1))
surface$chebyshev <- vapply(parameter_grid, function(par) {
  wpl_surrogate_reduced(unlist(par, use.names = FALSE), data = data, surrogate_fit = cheb_bundle$fit, fixed_beta = sim$true$beta, fixed_nu = sim$true$nu)
}, numeric(1))

reference_opt <- optimize_reduced_wpl(
  objective = wpl_reference_reduced,
  start = c(delta = 5.5, b = 0.30),
  lower = c(1.0, 0.05),
  upper = c(20.0, 1.5),
  data = data,
  fixed_beta = sim$true$beta,
  fixed_nu = sim$true$nu
)
spline_opt <- optimize_reduced_wpl(
  objective = wpl_surrogate_reduced,
  start = c(delta = 5.5, b = 0.30),
  lower = c(1.0, 0.05),
  upper = c(20.0, 1.5),
  data = data,
  surrogate_fit = spline_fit,
  fixed_beta = sim$true$beta,
  fixed_nu = sim$true$nu
)
cheb_opt <- optimize_reduced_wpl(
  objective = wpl_surrogate_reduced,
  start = c(delta = 5.5, b = 0.30),
  lower = c(1.0, 0.05),
  upper = c(20.0, 1.5),
  data = data,
  surrogate_fit = cheb_bundle$fit,
  fixed_beta = sim$true$beta,
  fixed_nu = sim$true$nu
)

primary_optima <- summarize_reduced_optima(
  reference_opt = reference_opt,
  spline_opt = spline_opt,
  chebyshev_opt = cheb_opt,
  data = data,
  fixed_beta = sim$true$beta,
  fixed_nu = sim$true$nu
)
save_reduced_contour_plot(
  surface,
  primary_optima,
  result_path("fig_reduced_wpl_contours.png")
)

primary_replicates <- lapply(seq_len(50L), function(rep_id) {
  rep_sim <- simulate_clayton_beta_example(
    n = 120L,
    beta = sim$true$beta,
    delta = sim$true$delta,
    b = sim$true$b,
    nu = sim$true$nu,
    seed = 20260409 + rep_id
  )
  rep_data <- prepare_reduced_clayton_data(rep_sim$y, rep_sim$coords, rep_sim$X, m = 2L)
  list(
    simulation = rep_sim,
    data = rep_data,
    reference_opt = optimize_reduced_wpl(
      objective = wpl_reference_reduced,
      start = c(delta = 5.5, b = 0.30),
      lower = c(1.0, 0.05),
      upper = c(20.0, 1.5),
      data = rep_data,
      fixed_beta = rep_sim$true$beta,
      fixed_nu = rep_sim$true$nu
    ),
    spline_opt = optimize_reduced_wpl(
      objective = wpl_surrogate_reduced,
      start = c(delta = 5.5, b = 0.30),
      lower = c(1.0, 0.05),
      upper = c(20.0, 1.5),
      data = rep_data,
      surrogate_fit = spline_fit,
      fixed_beta = rep_sim$true$beta,
      fixed_nu = rep_sim$true$nu
    ),
    chebyshev_opt = optimize_reduced_wpl(
      objective = wpl_surrogate_reduced,
      start = c(delta = 5.5, b = 0.30),
      lower = c(1.0, 0.05),
      upper = c(20.0, 1.5),
      data = rep_data,
      surrogate_fit = cheb_bundle$fit,
      fixed_beta = rep_sim$true$beta,
      fixed_nu = rep_sim$true$nu
    )
  )
})

replicate_gap_df <- summarize_replicate_gaps(primary_replicates, truth = sim$true)
replicate_gap_summary <- summarize_gap_table(
  replicate_gap_df,
  gap_columns = c("delta_gap_vs_reference", "b_gap_vs_reference"),
  error_columns = c("delta_error_vs_truth", "b_error_vs_truth")
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
replicate_gap_summary <- merge(replicate_gap_summary, criterion_summary, by = "method", all.x = TRUE, sort = FALSE)
replicate_gap_summary <- replicate_gap_summary[match(c("reference", "spline", "chebyshev"), replicate_gap_summary$method), , drop = FALSE]

write.csv(primary_optima, result_path("parameter_gap_summary.csv"), row.names = FALSE)
write.csv(replicate_gap_df, result_path("parameter_gap_replicates.csv"), row.names = FALSE)
write.csv(replicate_gap_summary, result_path("parameter_gap_replicate_summary.csv"), row.names = FALSE)

saveRDS(
  list(
    simulation = sim,
    illustrative_seed = 20260409L,
    surface = surface,
    primary_optima = primary_optima,
    replicate_gaps = replicate_gap_df,
    replicate_gap_summary = replicate_gap_summary
  ),
  result_path("05_reduced_wpl_comparison.rds")
)

