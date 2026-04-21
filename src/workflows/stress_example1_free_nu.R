script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

spline_fit <- readRDS(result_path("02_spline_fit_stress.rds"))
cheb_bundle <- readRDS(result_path("03_chebyshev_fit_stress.rds"))

run_opt_with_time <- function(objective, start, lower, upper, ...) {
  elapsed <- system.time({
    opt <- optimize_reduced_wpl(objective, start = start, lower = lower, upper = upper, ...)
  })[["elapsed"]]
  list(opt = opt, elapsed = elapsed)
}

summarize_stress_optima <- function(reference_opt, spline_opt, chebyshev_opt, data, fixed_beta) {
  optima <- rbind(
    data.frame(
      method = "reference",
      delta = unname(reference_opt$par[1]),
      b = unname(reference_opt$par[2]),
      nu = unname(reference_opt$par[3]),
      converged = identical(reference_opt$convergence, 0L),
      objective_min = unname(reference_opt$value),
      stringsAsFactors = FALSE
    ),
    data.frame(
      method = "spline",
      delta = unname(spline_opt$par[1]),
      b = unname(spline_opt$par[2]),
      nu = unname(spline_opt$par[3]),
      converged = identical(spline_opt$convergence, 0L),
      objective_min = unname(spline_opt$value),
      stringsAsFactors = FALSE
    ),
    data.frame(
      method = "chebyshev",
      delta = unname(chebyshev_opt$par[1]),
      b = unname(chebyshev_opt$par[2]),
      nu = unname(chebyshev_opt$par[3]),
      converged = identical(chebyshev_opt$convergence, 0L),
      objective_min = unname(chebyshev_opt$value),
      stringsAsFactors = FALSE
    )
  )

  reference_values <- data.frame(
    method = optima$method,
    reference_criterion = c(
      wpl_reference_stress(reference_opt$par, data, fixed_beta = fixed_beta),
      wpl_reference_stress(spline_opt$par, data, fixed_beta = fixed_beta),
      wpl_reference_stress(chebyshev_opt$par, data, fixed_beta = fixed_beta)
    ),
    stringsAsFactors = FALSE
  )

  gap_df <- merge(optima, reference_values, by = "method", sort = FALSE)
  ref_row <- gap_df[gap_df$method == "reference", ]
  gap_df$delta_gap_vs_reference <- gap_df$delta - ref_row$delta
  gap_df$b_gap_vs_reference <- gap_df$b - ref_row$b
  gap_df$nu_gap_vs_reference <- gap_df$nu - ref_row$nu
  gap_df$criterion_gap_exact <- ref_row$reference_criterion - gap_df$reference_criterion
  gap_df$criterion_gap_exact_pct <- 100 * gap_df$criterion_gap_exact / pmax(abs(ref_row$reference_criterion), 1e-12)
  gap_df
}

summarize_stress_replicates <- function(replicates, truth) {
  rows <- lapply(seq_along(replicates), function(k) {
    rep_obj <- replicates[[k]]
    ref_par <- rep_obj$reference_opt$opt$par
    ref_converged <- identical(rep_obj$reference_opt$opt$convergence, 0L)
    ref_exact <- if (ref_converged) {
      wpl_reference_stress(ref_par, rep_obj$data, fixed_beta = truth$beta)
    } else {
      NA_real_
    }
    methods <- list(
      reference = rep_obj$reference_opt,
      spline = rep_obj$spline_opt,
      chebyshev = rep_obj$chebyshev_opt
    )

    do.call(
      rbind,
      lapply(names(methods), function(method) {
        opt_obj <- methods[[method]]
        opt <- opt_obj$opt
        exact_value <- wpl_reference_stress(opt$par, rep_obj$data, fixed_beta = truth$beta)
        data.frame(
          replicate = k,
          method = method,
          reference_converged = ref_converged,
          converged = identical(opt$convergence, 0L),
          delta = unname(opt$par[1]),
          b = unname(opt$par[2]),
          nu = unname(opt$par[3]),
          delta_gap_vs_reference = if (ref_converged) unname(opt$par[1] - ref_par[1]) else NA_real_,
          b_gap_vs_reference = if (ref_converged) unname(opt$par[2] - ref_par[2]) else NA_real_,
          nu_gap_vs_reference = if (ref_converged) unname(opt$par[3] - ref_par[3]) else NA_real_,
          delta_error_vs_truth = unname(opt$par[1] - truth$delta),
          b_error_vs_truth = unname(opt$par[2] - truth$b),
          nu_error_vs_truth = unname(opt$par[3] - truth$nu),
          criterion_gap_exact = if (ref_converged) unname(ref_exact - exact_value) else NA_real_,
          criterion_gap_exact_pct = if (ref_converged) unname(100 * (ref_exact - exact_value) / pmax(abs(ref_exact), 1e-12)) else NA_real_,
          abs_criterion_gap_exact = if (ref_converged) abs(unname(ref_exact - exact_value)) else NA_real_,
          abs_criterion_gap_exact_pct = if (ref_converged) abs(unname(100 * (ref_exact - exact_value) / pmax(abs(ref_exact), 1e-12))) else NA_real_,
          optimization_time = opt_obj$elapsed,
          stringsAsFactors = FALSE
        )
      })
    )
  })
  do.call(rbind, rows)
}

summarize_stress_table <- function(df) {
  split_df <- split(df, df$method)
  out <- do.call(
    rbind,
    lapply(split_df, function(sub_df) {
      valid_ref <- sub_df$reference_converged
      data.frame(
        method = sub_df$method[1],
        n_total = nrow(sub_df),
        n_reference_converged = sum(valid_ref),
        convergence_rate = mean(sub_df$converged),
        mean_abs_delta_gap_vs_reference = mean(abs(sub_df$delta_gap_vs_reference), na.rm = TRUE),
        sd_delta_gap_vs_reference = stats::sd(sub_df$delta_gap_vs_reference, na.rm = TRUE),
        mean_abs_b_gap_vs_reference = mean(abs(sub_df$b_gap_vs_reference), na.rm = TRUE),
        sd_b_gap_vs_reference = stats::sd(sub_df$b_gap_vs_reference, na.rm = TRUE),
        mean_abs_nu_gap_vs_reference = mean(abs(sub_df$nu_gap_vs_reference), na.rm = TRUE),
        sd_nu_gap_vs_reference = stats::sd(sub_df$nu_gap_vs_reference, na.rm = TRUE),
        mean_abs_delta_error_vs_truth = mean(abs(sub_df$delta_error_vs_truth)),
        sd_delta_error_vs_truth = stats::sd(sub_df$delta_error_vs_truth),
        mean_abs_b_error_vs_truth = mean(abs(sub_df$b_error_vs_truth)),
        sd_b_error_vs_truth = stats::sd(sub_df$b_error_vs_truth),
        mean_abs_nu_error_vs_truth = mean(abs(sub_df$nu_error_vs_truth)),
        sd_nu_error_vs_truth = stats::sd(sub_df$nu_error_vs_truth),
        mean_abs_exact_criterion_gap = mean(sub_df$abs_criterion_gap_exact, na.rm = TRUE),
        mean_abs_exact_criterion_gap_pct = mean(sub_df$abs_criterion_gap_exact_pct, na.rm = TRUE),
        mean_time = mean(sub_df$optimization_time),
        sd_time = stats::sd(sub_df$optimization_time),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out[match(c("reference", "spline", "chebyshev"), out$method), , drop = FALSE]
}

sim <- simulate_clayton_beta_example(
  n = 120L,
  beta = 0,
  delta = 6,
  b = 0.35,
  nu = 5,
  seed = 20263200L
)

start <- c(delta = 5.5, b = 0.30, nu = 5.0)
lower <- c(1.0, 0.05, 3.0)
upper <- c(20.0, 1.5, 5.0)

illustrative_pool <- lapply(0:20, function(offset) {
  illustrative_sim <- simulate_clayton_beta_example(
    n = 120L,
    beta = sim$true$beta,
    delta = sim$true$delta,
    b = sim$true$b,
    nu = sim$true$nu,
    seed = 20263200L + offset
  )
  illustrative_data <- prepare_reduced_clayton_data(illustrative_sim$y, illustrative_sim$coords, illustrative_sim$X, m = 2L)
  illustrative_reference <- run_opt_with_time(
    objective = wpl_reference_stress,
    start = start,
    lower = lower,
    upper = upper,
    data = illustrative_data,
    fixed_beta = illustrative_sim$true$beta
  )
  list(
    seed = 20263200L + offset,
    simulation = illustrative_sim,
    data = illustrative_data,
    reference_opt = illustrative_reference
  )
})
illustrative_choice <- illustrative_pool[[which(vapply(illustrative_pool, function(obj) identical(obj$reference_opt$opt$convergence, 0L), logical(1)))[1]]]
sim <- illustrative_choice$simulation
data <- illustrative_choice$data
reference_opt <- illustrative_choice$reference_opt
spline_opt <- run_opt_with_time(
  objective = wpl_surrogate_stress,
  start = start,
  lower = lower,
  upper = upper,
  data = data,
  surrogate_fit = spline_fit,
  fixed_beta = sim$true$beta
)
chebyshev_opt <- run_opt_with_time(
  objective = wpl_surrogate_stress,
  start = start,
  lower = lower,
  upper = upper,
  data = data,
  surrogate_fit = cheb_bundle$fit,
  fixed_beta = sim$true$beta
)

illustrative_optima <- summarize_stress_optima(
  reference_opt = reference_opt$opt,
  spline_opt = spline_opt$opt,
  chebyshev_opt = chebyshev_opt$opt,
  data = data,
  fixed_beta = sim$true$beta
)
illustrative_optima$optimization_time <- c(reference_opt$elapsed, spline_opt$elapsed, chebyshev_opt$elapsed)

replicates <- lapply(seq_len(20L), function(rep_id) {
  rep_sim <- simulate_clayton_beta_example(
    n = 120L,
    beta = sim$true$beta,
    delta = sim$true$delta,
    b = sim$true$b,
    nu = sim$true$nu,
    seed = 20263200L + rep_id
  )
  rep_data <- prepare_reduced_clayton_data(rep_sim$y, rep_sim$coords, rep_sim$X, m = 2L)
  list(
    seed = 20263200L + rep_id,
    data = rep_data,
    reference_opt = run_opt_with_time(
      objective = wpl_reference_stress,
      start = start,
      lower = lower,
      upper = upper,
      data = rep_data,
      fixed_beta = rep_sim$true$beta
    ),
    spline_opt = run_opt_with_time(
      objective = wpl_surrogate_stress,
      start = start,
      lower = lower,
      upper = upper,
      data = rep_data,
      surrogate_fit = spline_fit,
      fixed_beta = rep_sim$true$beta
    ),
    chebyshev_opt = run_opt_with_time(
      objective = wpl_surrogate_stress,
      start = start,
      lower = lower,
      upper = upper,
      data = rep_data,
      surrogate_fit = cheb_bundle$fit,
      fixed_beta = rep_sim$true$beta
    )
  )
})

replicate_df <- summarize_stress_replicates(replicates, truth = sim$true)
summary_df <- summarize_stress_table(replicate_df)

write.csv(illustrative_optima, result_path("example1_stress_nu_free_optima.csv"), row.names = FALSE)
write.csv(replicate_df, result_path("example1_stress_nu_free_replicates.csv"), row.names = FALSE)
write.csv(summary_df, result_path("example1_stress_nu_free_summary.csv"), row.names = FALSE)

