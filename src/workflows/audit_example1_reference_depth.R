script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

spline_fit <- readRDS(result_path("02_spline_fit.rds"))
cheb_bundle <- readRDS(result_path("03_chebyshev_fit.rds"))

reference_caps <- c(400L, 700L, 900L)
subset_ids <- seq_len(10L)

fmt_num <- function(x, digits = 4L) {
  ifelse(is.na(x), "---", formatC(x, format = "f", digits = digits))
}

fmt_metric <- function(x, digits = 4L) {
  ifelse(
    is.na(x),
    "---",
    ifelse(abs(x) < 1e-3, formatC(x, format = "e", digits = 2), formatC(x, format = "f", digits = digits))
  )
}

quantile95 <- function(x) {
  stats::quantile(x, probs = 0.95, names = FALSE, na.rm = TRUE, type = 7)
}

bind_rows_fill <- function(dfs) {
  dfs <- dfs[!vapply(dfs, is.null, logical(1))]
  if (length(dfs) == 0L) {
    return(data.frame())
  }
  all_names <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  dfs <- lapply(dfs, function(df) {
    missing <- setdiff(all_names, names(df))
    for (name in missing) {
      df[[name]] <- NA
    }
    df[all_names]
  })
  do.call(rbind, dfs)
}

benchmark_specs <- list(
  reduced = list(
    truth = list(beta = 0, delta = 6, b = 0.35, nu = 6),
    seed_fun = function(rep_id) 20260409L + rep_id,
    simulate = function(seed) {
      simulate_clayton_beta_example(
        n = 120L,
        beta = 0,
        delta = 6,
        b = 0.35,
        nu = 6,
        seed = seed
      )
    },
    prepare = function(sim) prepare_reduced_clayton_data(sim$y, sim$coords, sim$X, m = 2L),
    reference_objective = wpl_reference_reduced,
    spline_objective = wpl_surrogate_reduced,
    chebyshev_objective = wpl_surrogate_reduced,
    start = c(delta = 5.5, b = 0.30),
    lower = c(1.0, 0.05),
    upper = c(20.0, 1.5),
    reference_args = function(sim, data) list(data = data, fixed_beta = sim$true$beta, fixed_nu = sim$true$nu),
    surrogate_args = function(sim, data, fit) list(data = data, surrogate_fit = fit, fixed_beta = sim$true$beta, fixed_nu = sim$true$nu),
    parameter_names = c("delta", "b"),
    beta_parameter_names = character(0)
  ),
  practical = list(
    truth = list(beta = c(-0.35, 0.9), delta = 6, b = 0.35, nu = 6),
    seed_fun = function(rep_id) 20261000L + rep_id,
    simulate = function(seed) {
      simulate_clayton_beta_example(
        n = 120L,
        beta = c(-0.35, 0.9),
        delta = 6,
        b = 0.35,
        nu = 6,
        seed = seed
      )
    },
    prepare = function(sim) prepare_reduced_clayton_data(sim$y, sim$coords, sim$X, m = 2L),
    reference_objective = wpl_reference_practical,
    spline_objective = wpl_surrogate_practical,
    chebyshev_objective = wpl_surrogate_practical,
    start = c(beta_intercept = 0, beta_x1 = 0, delta = 6, b = 0.30),
    lower = c(-5, -5, 1.0, 0.05),
    upper = c(5, 5, 20.0, 1.5),
    reference_args = function(sim, data) list(data = data, fixed_nu = sim$true$nu),
    surrogate_args = function(sim, data, fit) list(data = data, surrogate_fit = fit, fixed_nu = sim$true$nu),
    parameter_names = c("beta_intercept", "beta_x1", "delta", "b"),
    beta_parameter_names = c("beta_intercept", "beta_x1")
  )
)

collect_reference_row <- function(benchmark, replicate, seed, cap, run_obj, baseline_par = NULL, objective, objective_args) {
  par <- unname(run_obj$opt$par)
  out <- data.frame(
    benchmark = benchmark,
    replicate = replicate,
    seed = seed,
    reference_cap = cap,
    method = "reference",
    converged = identical(run_obj$opt$convergence, 0L),
    optimization_time = run_obj$elapsed,
    fncount = run_obj$fncount,
    grcount = run_obj$grcount,
    criterion_value = if (identical(run_obj$opt$convergence, 0L)) {
      do.call(objective, c(list(theta = par), objective_args, list(max_m = cap, max_n = cap)))
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
  for (j in seq_along(par)) {
    out[[paste0("theta_", j)]] <- par[j]
    out[[paste0("abs_shift_vs_cap400_", j)]] <- if (is.null(baseline_par)) 0 else abs(par[j] - baseline_par[j])
  }
  out
}

collect_surrogate_row <- function(benchmark, replicate, seed, cap, method, run_obj, reference_row, objective, objective_args) {
  par <- unname(run_obj$opt$par)
  ref_par <- unlist(reference_row[paste0("theta_", seq_along(par))], use.names = FALSE)
  converged <- identical(run_obj$opt$convergence, 0L)
  reference_converged <- isTRUE(reference_row$converged)
  criterion_at_theta <- if (converged && reference_converged) {
    do.call(objective, c(list(theta = par), objective_args, list(max_m = cap, max_n = cap)))
  } else {
    NA_real_
  }
  criterion_at_ref <- if (reference_converged) reference_row$criterion_value else NA_real_

  out <- data.frame(
    benchmark = benchmark,
    replicate = replicate,
    seed = seed,
    reference_cap = cap,
    method = method,
    converged = converged,
    reference_converged = reference_converged,
    optimization_time = run_obj$elapsed,
    fncount = run_obj$fncount,
    grcount = run_obj$grcount,
    criterion_value = criterion_at_theta,
    criterion_gap_exact = if (reference_converged && converged) criterion_at_ref - criterion_at_theta else NA_real_,
    criterion_gap_exact_pct = if (reference_converged && converged) {
      100 * abs(criterion_at_ref - criterion_at_theta) / pmax(abs(criterion_at_ref), 1e-8)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
  for (j in seq_along(par)) {
    out[[paste0("theta_", j)]] <- par[j]
    out[[paste0("abs_gap_vs_reference_", j)]] <- if (reference_converged && converged) abs(par[j] - ref_par[j]) else NA_real_
  }
  out
}

run_depth_check <- function(benchmark_name, spec) {
  rows <- vector("list", length(subset_ids))
  for (idx in seq_along(subset_ids)) {
    rep_id <- subset_ids[[idx]]
    seed <- spec$seed_fun(rep_id)
    sim <- spec$simulate(seed)
    data <- spec$prepare(sim)

    spline_run <- do.call(
      run_opt_with_timing,
      c(
        list(
          objective = spec$spline_objective,
          start = spec$start,
          lower = spec$lower,
          upper = spec$upper
        ),
        spec$surrogate_args(sim, data, spline_fit)
      )
    )
    cheb_run <- do.call(
      run_opt_with_timing,
      c(
        list(
          objective = spec$chebyshev_objective,
          start = spec$start,
          lower = spec$lower,
          upper = spec$upper
        ),
        spec$surrogate_args(sim, data, cheb_bundle$fit)
      )
    )

    reference_runs <- lapply(reference_caps, function(cap) {
      do.call(
        run_opt_with_timing,
        c(
          list(
            objective = spec$reference_objective,
            start = spec$start,
            lower = spec$lower,
            upper = spec$upper
          ),
          spec$reference_args(sim, data),
          list(max_m = cap, max_n = cap)
        )
      )
    })
    names(reference_runs) <- as.character(reference_caps)
    baseline_par <- unname(reference_runs[[as.character(reference_caps[1])]]$opt$par)
    objective_args <- spec$reference_args(sim, data)

    reference_rows <- do.call(
      rbind,
      lapply(reference_caps, function(cap) {
        collect_reference_row(
          benchmark = benchmark_name,
          replicate = rep_id,
          seed = seed,
          cap = cap,
          run_obj = reference_runs[[as.character(cap)]],
          baseline_par = if (identical(cap, reference_caps[1])) NULL else baseline_par,
          objective = spec$reference_objective,
          objective_args = objective_args
        )
      })
    )

    cap400_row <- reference_rows[reference_rows$reference_cap == reference_caps[1], , drop = FALSE]
    if (isTRUE(cap400_row$converged)) {
      cap400_values <- vapply(reference_caps, function(cap) {
        do.call(
          spec$reference_objective,
          c(
            list(theta = baseline_par),
            objective_args,
            list(max_m = cap, max_n = cap)
          )
        )
      }, numeric(1))
      reference_rows$criterion_improvement_over_cap400_pct <- ifelse(
        reference_rows$reference_cap == reference_caps[1] | !reference_rows$converged,
        0,
        100 * abs(reference_rows$criterion_value - cap400_values[match(reference_rows$reference_cap, reference_caps)]) /
          pmax(abs(reference_rows$criterion_value), 1e-8)
      )
    } else {
      reference_rows$criterion_improvement_over_cap400_pct <- NA_real_
    }

    surrogate_rows <- do.call(
      rbind,
      lapply(reference_caps, function(cap) {
        reference_row <- reference_rows[reference_rows$reference_cap == cap, , drop = FALSE]
        rbind(
          collect_surrogate_row(
            benchmark = benchmark_name,
            replicate = rep_id,
            seed = seed,
            cap = cap,
            method = "spline",
            run_obj = spline_run,
            reference_row = reference_row,
            objective = spec$reference_objective,
            objective_args = objective_args
          ),
          collect_surrogate_row(
            benchmark = benchmark_name,
            replicate = rep_id,
            seed = seed,
            cap = cap,
            method = "chebyshev",
            run_obj = cheb_run,
            reference_row = reference_row,
            objective = spec$reference_objective,
            objective_args = objective_args
          )
        )
      })
    )

    rows[[idx]] <- bind_rows_fill(list(reference_rows, surrogate_rows))
  }

  bind_rows_fill(rows)
}

depth_replicates <- bind_rows_fill(
  lapply(names(benchmark_specs), function(name) run_depth_check(name, benchmark_specs[[name]]))
)

summarize_reference_depth <- function(df, spec_map) {
  reference_df <- df[df$method == "reference", , drop = FALSE]
  rows <- do.call(
    rbind,
    lapply(split(reference_df, list(reference_df$benchmark, reference_df$reference_cap), drop = TRUE), function(sub_df) {
      benchmark_name <- sub_df$benchmark[1]
      param_cols <- paste0("abs_shift_vs_cap400_", seq_along(spec_map[[benchmark_name]]$parameter_names))
      per_dataset_max_shift <- apply(sub_df[param_cols], 1, function(row) max(row, na.rm = TRUE))
      if (all(!is.finite(per_dataset_max_shift))) {
        per_dataset_max_shift <- rep(NA_real_, nrow(sub_df))
      }
      data.frame(
        benchmark = benchmark_name,
        reference_cap = sub_df$reference_cap[1],
        n = nrow(sub_df),
        convergence_rate = mean(sub_df$converged),
        mean_max_abs_shift_vs_cap400 = mean(per_dataset_max_shift, na.rm = TRUE),
        p95_max_abs_shift_vs_cap400 = quantile95(per_dataset_max_shift),
        max_abs_shift_vs_cap400 = max(per_dataset_max_shift, na.rm = TRUE),
        mean_exact_improvement_over_cap400_pct = mean(sub_df$criterion_improvement_over_cap400_pct, na.rm = TRUE),
        p95_exact_improvement_over_cap400_pct = quantile95(sub_df$criterion_improvement_over_cap400_pct),
        max_exact_improvement_over_cap400_pct = max(sub_df$criterion_improvement_over_cap400_pct, na.rm = TRUE),
        mean_time = mean(sub_df$optimization_time, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(rows) <- NULL
  rows[order(match(rows$benchmark, c("reduced", "practical")), rows$reference_cap), , drop = FALSE]
}

summarize_surrogate_depth <- function(df, spec_map) {
  surrogate_df <- df[df$method != "reference", , drop = FALSE]
  rows <- do.call(
    rbind,
    lapply(split(surrogate_df, list(surrogate_df$benchmark, surrogate_df$reference_cap, surrogate_df$method), drop = TRUE), function(sub_df) {
      benchmark_name <- sub_df$benchmark[1]
      data.frame(
        benchmark = benchmark_name,
        reference_cap = sub_df$reference_cap[1],
        method = sub_df$method[1],
        n = nrow(sub_df),
        matched_convergence_rate = mean(sub_df$converged & sub_df$reference_converged),
        mean_abs_beta_gap_vs_reference = if (benchmark_name == "practical") {
          mean(c(sub_df$abs_gap_vs_reference_1, sub_df$abs_gap_vs_reference_2), na.rm = TRUE)
        } else {
          NA_real_
        },
        mean_abs_delta_gap_vs_reference = if (benchmark_name == "reduced") {
          mean(sub_df$abs_gap_vs_reference_1, na.rm = TRUE)
        } else {
          mean(sub_df$abs_gap_vs_reference_3, na.rm = TRUE)
        },
        mean_abs_b_gap_vs_reference = if (benchmark_name == "reduced") {
          mean(sub_df$abs_gap_vs_reference_2, na.rm = TRUE)
        } else {
          mean(sub_df$abs_gap_vs_reference_4, na.rm = TRUE)
        },
        mean_exact_loss_pct = mean(sub_df$criterion_gap_exact_pct, na.rm = TRUE),
        mean_time = mean(sub_df$optimization_time, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(rows) <- NULL
  rows[order(match(rows$benchmark, c("reduced", "practical")), rows$reference_cap, match(rows$method, c("spline", "chebyshev"))), , drop = FALSE]
}

reference_summary <- summarize_reference_depth(depth_replicates, benchmark_specs)
surrogate_summary <- summarize_surrogate_depth(depth_replicates, benchmark_specs)

write.csv(depth_replicates, result_path("example1_reference_depth_subset_replicates.csv"), row.names = FALSE)
write.csv(reference_summary, result_path("example1_reference_depth_reference_summary.csv"), row.names = FALSE)
write.csv(surrogate_summary, result_path("example1_reference_depth_surrogate_summary.csv"), row.names = FALSE)

