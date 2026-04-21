classify_transformed_region <- function(grid_df, box) {
  q_low_cut <- box$q[1] + 0.20 * diff(box$q)
  q_stable_cut <- box$q[1] + 0.35 * diff(box$q)
  t_high_cut <- box$t[1] + 0.75 * diff(box$t)
  t_stable_cut <- box$t[1] + 0.50 * diff(box$t)

  region <- ifelse(
    grid_df$t >= t_high_cut & grid_df$q <= q_low_cut,
    "high_t_low_q",
    ifelse(
      grid_df$t >= t_high_cut,
      "high_dependence",
      ifelse(
        grid_df$q <= q_low_cut,
        "low_separation",
        ifelse(
          grid_df$t <= t_stable_cut & grid_df$q >= q_stable_cut,
          "interior_stable",
          "transition"
        )
      )
    )
  )

  transform(
    grid_df,
    region = factor(
      region,
      levels = c("interior_stable", "transition", "low_separation", "high_dependence", "high_t_low_q")
    ),
    q_low_cut = q_low_cut,
    t_high_cut = t_high_cut
  )
}

sample_transformed_points_by_region <- function(box, n_per_region = 150L, seed = 20260410L, max_draws = 250000L, batch_size = 4000L) {
  target_regions <- c("interior_stable", "transition", "low_separation", "high_dependence", "high_t_low_q")
  collected <- setNames(vector("list", length(target_regions)), target_regions)
  counts <- setNames(integer(length(target_regions)), target_regions)
  set.seed(seed)
  draws <- 0L

  while (any(counts < n_per_region) && draws < max_draws) {
    batch_n <- min(batch_size, max_draws - draws)
    batch <- data.frame(
      x = stats::runif(batch_n, min = box$x[1], max = box$x[2]),
      q = stats::runif(batch_n, min = box$q[1], max = box$q[2]),
      t = stats::runif(batch_n, min = box$t[1], max = box$t[2]),
      stringsAsFactors = FALSE
    )
    classified <- classify_transformed_region(batch, box)
    draws <- draws + batch_n

    for (region_name in target_regions) {
      need <- n_per_region - counts[[region_name]]
      if (need <= 0L) {
        next
      }
      region_rows <- classified[classified$region == region_name, c("x", "q", "t", "region"), drop = FALSE]
      if (!nrow(region_rows)) {
        next
      }
      take <- min(need, nrow(region_rows))
      collected[[region_name]] <- rbind(collected[[region_name]], region_rows[seq_len(take), , drop = FALSE])
      counts[[region_name]] <- nrow(collected[[region_name]])
    }
  }

  if (any(counts < n_per_region)) {
    stop("Could not collect enough transformed evaluation points in every region.", call. = FALSE)
  }

  sampled <- do.call(rbind, collected)
  rownames(sampled) <- NULL
  sampled$region <- factor(sampled$region, levels = target_regions)
  sampled
}

summarize_region_metric <- function(df, value_col, tolerances = c(1e-8, 1e-6, 1e-4)) {
  out <- do.call(
    rbind,
    lapply(split(df, df$region), function(sub_df) {
      values <- sub_df[[value_col]]
      data.frame(
        region = as.character(sub_df$region[1]),
        n = nrow(sub_df),
        mean = mean(values),
        median = median(values),
        max = max(values),
        share_gt_1e8 = mean(values > tolerances[1]),
        share_gt_1e6 = mean(values > tolerances[2]),
        share_gt_1e4 = mean(values > tolerances[3]),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out
}

format_scientific_floor <- function(x, floor = 1e-15, digits = 2) {
  ifelse(
    x <= floor,
    sprintf("<%.0e", floor),
    formatC(x, format = "e", digits = digits)
  )
}

surrogate_local_error_summary_generic <- function(fit, evaluation_grid, reference_fun, surrogate_fun = NULL, reference_step = 1e-4) {
  if (is.null(surrogate_fun)) {
    surrogate_fun <- function(z, fit) {
      evaluate_transformed_surrogate(fit, z[1], z[2], z[3])
    }
  }

  reference_values <- vapply(
    seq_len(nrow(evaluation_grid)),
    function(i) reference_fun(as.numeric(evaluation_grid[i, c("x", "q", "t")])),
    numeric(1)
  )
  surrogate_values <- vapply(
    seq_len(nrow(evaluation_grid)),
    function(i) surrogate_fun(as.numeric(evaluation_grid[i, c("x", "q", "t")]), fit),
    numeric(1)
  )

  grad_error <- numeric(nrow(evaluation_grid))
  hess_error <- numeric(nrow(evaluation_grid))
  value_error <- abs(surrogate_values - reference_values)

  for (i in seq_len(nrow(evaluation_grid))) {
    z0 <- as.numeric(evaluation_grid[i, c("x", "q", "t")])
    grad_reference <- finite_difference_gradient(reference_fun, z0, step = reference_step)
    grad_sur <- finite_difference_gradient(surrogate_fun, z0, step = reference_step, fit = fit)
    hess_reference <- finite_difference_hessian(reference_fun, z0, step = reference_step)
    hess_sur <- finite_difference_hessian(surrogate_fun, z0, step = reference_step, fit = fit)
    grad_error[i] <- sqrt(sum((grad_sur - grad_reference)^2))
    hess_error[i] <- sqrt(sum((hess_sur - hess_reference)^2))
  }

  list(
    pointwise = cbind(
      evaluation_grid,
      reference = reference_values,
      surrogate = surrogate_values,
      value_error = value_error,
      grad_error = grad_error,
      hess_error = hess_error
    ),
    summary = list(
      sup_value_error = max(value_error),
      sup_grad_error = max(grad_error),
      sup_hess_error = max(hess_error),
      mean_value_error = mean(value_error),
      mean_grad_error = mean(grad_error),
      mean_hess_error = mean(hess_error)
    )
  )
}

surrogate_local_error_summary <- function(fit, nu, evaluation_grid, reference_step = 1e-4, reference_tol = 1e-12, max_m = 700L, max_n = 700L) {
  surrogate_local_error_summary_generic(
    fit = fit,
    evaluation_grid = evaluation_grid,
    reference_fun = function(z) {
      reference_transformed_clayton(z[1], z[2], z[3], nu = nu, log = TRUE, tol = reference_tol, max_m = max_m, max_n = max_n)
    },
    reference_step = reference_step
  )
}

reference_vs_geomodels_summary <- function(grid_df, nu, box, tolerances = c(1e-6, 1e-4, 1e-2)) {
  uvrt <- inverse_transform_clayton_inputs(grid_df$x, grid_df$q, grid_df$t, branch = "positive")
  external_reference <- geomodels_external_reference_clayton(uvrt$u, uvrt$v, uvrt$r, nu = nu, log = TRUE)
  series_reference <- reference_transformed_clayton(grid_df$x, grid_df$q, grid_df$t, nu = nu, log = TRUE)
  pointwise <- transform(
    classify_transformed_region(grid_df, box),
    series_reference = series_reference,
    geomodels_reference = external_reference,
    abs_diff = abs(series_reference - external_reference),
    rel_diff = abs(series_reference - external_reference) / pmax(1, abs(series_reference))
  )

  by_region <- do.call(
    rbind,
    lapply(split(pointwise, pointwise$region), function(df) {
      data.frame(
        region = as.character(df$region[1]),
        n = nrow(df),
        mean_abs_diff = mean(df$abs_diff),
        median_abs_diff = median(df$abs_diff),
        max_abs_diff = max(df$abs_diff),
        share_gt_1e6 = mean(df$abs_diff > tolerances[1]),
        share_gt_1e4 = mean(df$abs_diff > tolerances[2]),
        share_gt_1e2 = mean(df$abs_diff > tolerances[3]),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(by_region) <- NULL

  worst_points <- utils::head(pointwise[order(pointwise$abs_diff, decreasing = TRUE), ], 10L)
  rule_text <- sprintf(
    "Use the primary benchmark outside the stress corner t > %.3f and q < %.2f; that corner is retained as a numerical audit region.",
    unique(pointwise$t_high_cut)[1],
    unique(pointwise$q_low_cut)[1]
  )

  list(
    pointwise = pointwise,
    summary = list(
      mean_abs_diff = mean(pointwise$abs_diff),
      median_abs_diff = median(pointwise$abs_diff),
      max_abs_diff = max(pointwise$abs_diff),
      mean_rel_diff = mean(pointwise$rel_diff),
      share_gt_1e6 = mean(pointwise$abs_diff > tolerances[1]),
      share_gt_1e4 = mean(pointwise$abs_diff > tolerances[2]),
      share_gt_1e2 = mean(pointwise$abs_diff > tolerances[3]),
      operational_rule = rule_text
    ),
    by_region = by_region,
    worst_points = worst_points,
    recommended_box = list(
      x = box$x,
      q = c(box$q[1] + 0.20 * diff(box$q), box$q[2]),
      t = c(box$t[1], box$t[1] + 0.75 * diff(box$t))
    )
  )
}

reference_truncation_stability <- function(grid_df, nu, box, truncation_levels = c(300L, 500L, 700L, 900L), tol = 1e-12) {
  pointwise <- classify_transformed_region(grid_df, box)
  pointwise$reference_value <- reference_transformed_clayton(
    pointwise$x,
    pointwise$q,
    pointwise$t,
    nu = nu,
    log = TRUE,
    tol = tol,
    max_m = max(truncation_levels),
    max_n = max(truncation_levels)
  )

  level_values <- lapply(truncation_levels, function(level) {
    reference_transformed_clayton(
      pointwise$x,
      pointwise$q,
      pointwise$t,
      nu = nu,
      log = TRUE,
      tol = tol,
      max_m = level,
      max_n = level
    )
  })
  names(level_values) <- paste0("m", truncation_levels)

  long_df <- do.call(
    rbind,
    lapply(seq_along(truncation_levels), function(k) {
      level <- truncation_levels[k]
      values <- level_values[[k]]
      data.frame(
        x = pointwise$x,
        q = pointwise$q,
        t = pointwise$t,
        region = pointwise$region,
        truncation_level = level,
        abs_diff_vs_max = abs(values - pointwise$reference_value),
        stringsAsFactors = FALSE
      )
    })
  )

  by_region <- do.call(
    rbind,
    lapply(split(long_df, list(long_df$region, long_df$truncation_level), drop = TRUE), function(sub_df) {
      data.frame(
        region = as.character(sub_df$region[1]),
        truncation_level = sub_df$truncation_level[1],
        n = nrow(sub_df),
        mean_abs_diff = mean(sub_df$abs_diff_vs_max),
        median_abs_diff = median(sub_df$abs_diff_vs_max),
        max_abs_diff = max(sub_df$abs_diff_vs_max),
        share_gt_1e6 = mean(sub_df$abs_diff_vs_max > 1e-6),
        share_gt_1e4 = mean(sub_df$abs_diff_vs_max > 1e-4),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(by_region) <- NULL

  list(
    pointwise = pointwise,
    truncation_pointwise = long_df,
    by_region = by_region,
    max_level = max(truncation_levels),
    operational_rule = sprintf(
      "Primary benchmark region excludes the stress corner t > %.3f and q < %.2f, where truncation differences concentrate.",
      unique(pointwise$t_high_cut)[1],
      unique(pointwise$q_low_cut)[1]
    )
  )
}

extract_optimum_row <- function(opt, method) {
  data.frame(
    method = method,
    delta = unname(opt$par[1]),
    b = unname(opt$par[2]),
    converged = identical(opt$convergence, 0L),
    objective_min = unname(opt$value),
    stringsAsFactors = FALSE
  )
}

summarize_reduced_optima <- function(reference_opt, spline_opt, chebyshev_opt, data, fixed_beta, fixed_nu) {
  optima <- rbind(
    extract_optimum_row(reference_opt, "reference"),
    extract_optimum_row(spline_opt, "spline"),
    extract_optimum_row(chebyshev_opt, "chebyshev")
  )

  reference_values <- data.frame(
    method = optima$method,
    reference_criterion = c(
      wpl_reference_reduced(reference_opt$par, data, fixed_beta = fixed_beta, fixed_nu = fixed_nu),
      wpl_reference_reduced(spline_opt$par, data, fixed_beta = fixed_beta, fixed_nu = fixed_nu),
      wpl_reference_reduced(chebyshev_opt$par, data, fixed_beta = fixed_beta, fixed_nu = fixed_nu)
    ),
    stringsAsFactors = FALSE
  )

  gap_df <- merge(optima, reference_values, by = "method", sort = FALSE)
  ref_row <- gap_df[gap_df$method == "reference", ]
  gap_df$delta_gap_vs_reference <- gap_df$delta - ref_row$delta
  gap_df$b_gap_vs_reference <- gap_df$b - ref_row$b
  gap_df
}

summarize_replicate_gaps <- function(replicates, truth) {
  rows <- lapply(seq_along(replicates), function(k) {
    rep_obj <- replicates[[k]]
    ref_par <- rep_obj$reference_opt$par
    ref_criterion <- if (!is.null(rep_obj$data)) {
      wpl_reference_reduced(ref_par, rep_obj$data, fixed_beta = truth$beta, fixed_nu = truth$nu)
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
        opt <- methods[[method]]
        data.frame(
          replicate = k,
          method = method,
          delta = unname(opt$par[1]),
          b = unname(opt$par[2]),
          converged = identical(opt$convergence, 0L),
          delta_gap_vs_reference = unname(opt$par[1] - ref_par[1]),
          b_gap_vs_reference = unname(opt$par[2] - ref_par[2]),
          delta_error_vs_truth = unname(opt$par[1] - truth$delta),
          b_error_vs_truth = unname(opt$par[2] - truth$b),
          criterion_gap_exact = if (!is.na(ref_criterion)) ref_criterion - wpl_reference_reduced(opt$par, rep_obj$data, fixed_beta = truth$beta, fixed_nu = truth$nu) else NA_real_,
          criterion_gap_exact_pct = if (!is.na(ref_criterion)) {
            100 * abs(ref_criterion - wpl_reference_reduced(opt$par, rep_obj$data, fixed_beta = truth$beta, fixed_nu = truth$nu)) /
              pmax(abs(ref_criterion), 1e-8)
          } else {
            NA_real_
          },
          stringsAsFactors = FALSE
        )
      })
    )
  })
  do.call(rbind, rows)
}

extract_named_optimum_row <- function(opt, method, parameter_names) {
  data.frame(
    method = method,
    stats::setNames(as.list(unname(opt$par)), parameter_names),
    converged = identical(opt$convergence, 0L),
    objective_min = unname(opt$value),
    stringsAsFactors = FALSE
  )
}

summarize_practical_optima <- function(reference_opt, spline_opt, chebyshev_opt, data, fixed_nu) {
  parameter_names <- practical_parameter_names(data)
  optima <- rbind(
    extract_named_optimum_row(reference_opt, "reference", parameter_names),
    extract_named_optimum_row(spline_opt, "spline", parameter_names),
    extract_named_optimum_row(chebyshev_opt, "chebyshev", parameter_names)
  )

  reference_values <- data.frame(
    method = optima$method,
    reference_criterion = c(
      wpl_reference_practical(reference_opt$par, data, fixed_nu = fixed_nu),
      wpl_reference_practical(spline_opt$par, data, fixed_nu = fixed_nu),
      wpl_reference_practical(chebyshev_opt$par, data, fixed_nu = fixed_nu)
    ),
    stringsAsFactors = FALSE
  )
  gap_df <- merge(optima, reference_values, by = "method", sort = FALSE)
  ref_row <- gap_df[gap_df$method == "reference", ]
  for (name in parameter_names) {
    gap_df[[paste0(name, "_gap_vs_reference")]] <- gap_df[[name]] - ref_row[[name]]
  }
  gap_df
}

summarize_practical_replicate_gaps <- function(replicates, truth, parameter_names) {
  rows <- lapply(seq_along(replicates), function(k) {
    rep_obj <- replicates[[k]]
    ref_par <- rep_obj$reference_opt$par
    ref_criterion <- if (!is.null(rep_obj$data)) {
      wpl_reference_practical(ref_par, rep_obj$data, fixed_nu = truth$nu)
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
        opt <- methods[[method]]
        row <- list(
          replicate = k,
          method = method,
          converged = identical(opt$convergence, 0L)
        )
        for (j in seq_along(parameter_names)) {
          row[[parameter_names[j]]] <- unname(opt$par[j])
          row[[paste0(parameter_names[j], "_gap_vs_reference")]] <- unname(opt$par[j] - ref_par[j])
        }
        row[["delta_error_vs_truth"]] <- unname(opt$par[length(parameter_names) - 1L] - truth$delta)
        row[["b_error_vs_truth"]] <- unname(opt$par[length(parameter_names)] - truth$b)
        beta_truth <- truth$beta
        for (j in seq_along(beta_truth)) {
          row[[paste0(parameter_names[j], "_error_vs_truth")]] <- unname(opt$par[j] - beta_truth[j])
        }
        row[["criterion_gap_exact"]] <- if (!is.na(ref_criterion)) {
          ref_criterion - wpl_reference_practical(opt$par, rep_obj$data, fixed_nu = truth$nu)
        } else {
          NA_real_
        }
        row[["criterion_gap_exact_pct"]] <- if (!is.na(ref_criterion)) {
          100 * abs(ref_criterion - wpl_reference_practical(opt$par, rep_obj$data, fixed_nu = truth$nu)) /
            pmax(abs(ref_criterion), 1e-8)
        } else {
          NA_real_
        }
        as.data.frame(row, stringsAsFactors = FALSE)
      })
    )
  })
  do.call(rbind, rows)
}

benchmark_one <- function(expr, times = 20L, inner = 1L) {
  expr_call <- substitute(expr)
  eval_env <- parent.frame()
  elapsed <- replicate(times, {
    system.time({
      for (k in seq_len(inner)) {
        eval(expr_call, envir = eval_env)
      }
    })[["elapsed"]] / inner
  })
  c(mean = mean(elapsed), median = median(elapsed), min = min(elapsed), max = max(elapsed))
}

extract_optim_count <- function(opt, kind) {
  counts <- opt$counts
  if (is.null(counts) || length(counts) == 0L) {
    return(NA_real_)
  }
  if (!is.null(names(counts)) && kind %in% names(counts)) {
    return(unname(counts[[kind]]))
  }
  if (identical(kind, "function") && length(counts) >= 1L) {
    return(unname(counts[[1]]))
  }
  if (identical(kind, "gradient") && length(counts) >= 2L) {
    return(unname(counts[[2]]))
  }
  NA_real_
}

run_opt_with_timing <- function(objective, start, lower, upper, ...) {
  elapsed <- system.time({
    opt <- optimize_reduced_wpl(objective, start = start, lower = lower, upper = upper, ...)
  })[["elapsed"]]

  list(
    opt = opt,
    elapsed = elapsed,
    fncount = extract_optim_count(opt, "function"),
    grcount = extract_optim_count(opt, "gradient"),
    convergence = opt$convergence,
    message = if (!is.null(opt$message)) opt$message else NA_character_
  )
}

summarize_timed_runs <- function(run_list) {
  elapsed <- vapply(run_list, function(run) run$elapsed, numeric(1))
  fncount <- vapply(run_list, function(run) run$fncount, numeric(1))
  grcount <- vapply(run_list, function(run) run$grcount, numeric(1))
  converged <- vapply(run_list, function(run) identical(run$opt$convergence, 0L), logical(1))

  c(
    mean = mean(elapsed),
    median = median(elapsed),
    min = min(elapsed),
    max = max(elapsed),
    convergence_rate = mean(converged),
    mean_fncount = mean(fncount, na.rm = TRUE),
    mean_grcount = mean(grcount, na.rm = TRUE)
  )
}

runtime_benchmark_summary <- function(data, truth, spline_fit, chebyshev_fit, optimization_times = 5L) {
  idx <- data$design$pair_index
  mu <- compute_beta_mean(data$X, truth$beta)
  alpha <- mu * truth$delta
  gamma <- (1 - mu) * truth$delta
  u <- pbeta(data$y, shape1 = alpha, shape2 = gamma)
  r <- wendland_r4(data$design$distances, truth$b)

  local <- rbind(
    reference = benchmark_one(reference_local_clayton(u[idx[, "i"]], u[idx[, "j"]], r, nu = truth$nu, log = TRUE), inner = 200L),
    spline = benchmark_one(predict_spline_surrogate(spline_fit, u[idx[, "i"]], u[idx[, "j"]], r), inner = 200L),
    chebyshev = benchmark_one(predict_chebyshev_surrogate(chebyshev_fit, u[idx[, "i"]], u[idx[, "j"]], r), inner = 200L)
  )

  criterion <- rbind(
    reference = benchmark_one(wpl_reference_reduced(c(truth$delta, truth$b), data, fixed_beta = truth$beta, fixed_nu = truth$nu), inner = 100L),
    spline = benchmark_one(wpl_surrogate_reduced(c(truth$delta, truth$b), data, surrogate_fit = spline_fit, fixed_beta = truth$beta, fixed_nu = truth$nu), inner = 100L),
    chebyshev = benchmark_one(wpl_surrogate_reduced(c(truth$delta, truth$b), data, surrogate_fit = chebyshev_fit, fixed_beta = truth$beta, fixed_nu = truth$nu), inner = 100L)
  )

  optimization_runs <- list(
    reference = replicate(
      optimization_times,
      run_opt_with_timing(
        wpl_reference_reduced,
        c(5.5, 0.30),
        c(1.0, 0.05),
        c(20.0, 1.5),
        data = data,
        fixed_beta = truth$beta,
        fixed_nu = truth$nu
      ),
      simplify = FALSE
    ),
    spline = replicate(
      optimization_times,
      run_opt_with_timing(
        wpl_surrogate_reduced,
        c(5.5, 0.30),
        c(1.0, 0.05),
        c(20.0, 1.5),
        data = data,
        surrogate_fit = spline_fit,
        fixed_beta = truth$beta,
        fixed_nu = truth$nu
      ),
      simplify = FALSE
    ),
    chebyshev = replicate(
      optimization_times,
      run_opt_with_timing(
        wpl_surrogate_reduced,
        c(5.5, 0.30),
        c(1.0, 0.05),
        c(20.0, 1.5),
        data = data,
        surrogate_fit = chebyshev_fit,
        fixed_beta = truth$beta,
        fixed_nu = truth$nu
      ),
      simplify = FALSE
    )
  )
  optimization <- do.call(
    rbind,
    lapply(optimization_runs, summarize_timed_runs)
  )

  list(local = local, criterion = criterion, optimization = optimization)
}

runtime_scaling_summary <- function(rows) {
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  df$speedup_vs_reference <- NA_real_
  for (n_val in unique(df$n)) {
    ref_mean <- df$mean_time[df$n == n_val & df$method == "reference"][1]
    df$speedup_vs_reference[df$n == n_val] <- ref_mean / df$mean_time[df$n == n_val]
  }
  df
}

summarize_gap_table <- function(df, gap_columns, error_columns = NULL) {
  split_df <- split(df, df$method)
  out <- do.call(
    rbind,
    lapply(split_df, function(sub_df) {
      stats_list <- list(method = sub_df$method[1], n = nrow(sub_df), convergence_rate = mean(sub_df$converged))
      for (col_name in gap_columns) {
        stats_list[[paste0(col_name, "_mean_abs")]] <- mean(abs(sub_df[[col_name]]))
        stats_list[[paste0(col_name, "_sd")]] <- stats::sd(sub_df[[col_name]])
      }
      if (!is.null(error_columns)) {
        for (col_name in error_columns) {
          stats_list[[paste0(col_name, "_mean_abs")]] <- mean(abs(sub_df[[col_name]]))
          stats_list[[paste0(col_name, "_sd")]] <- stats::sd(sub_df[[col_name]])
        }
      }
      as.data.frame(stats_list, stringsAsFactors = FALSE)
    })
  )
  rownames(out) <- NULL
  out
}

save_truncation_stability_summary_plot <- function(summary_df, file) {
  summary_df <- summary_df[summary_df$truncation_level != max(summary_df$truncation_level), , drop = FALSE]
  region_levels <- c("interior_stable", "transition", "low_separation", "high_dependence", "high_t_low_q")
  region_labels <- c("Interior stable", "Transition", "Low separation", "High dependence", "High t, low q")
  region_cols <- c(
    interior_stable = "#636363",
    transition = "#9ecae1",
    low_separation = "#6baed6",
    high_dependence = "#2171b5",
    high_t_low_q = "#cb181d"
  )
  region_pch <- c(
    interior_stable = 15,
    transition = 16,
    low_separation = 17,
    high_dependence = 18,
    high_t_low_q = 19
  )
  trunc_levels <- sort(unique(summary_df$truncation_level))
  y_floor <- 1e-15

  grDevices::png(file, width = 1400, height = 700, res = 140)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 0.68))

  for (metric in c("mean_abs_diff", "max_abs_diff")) {
    graphics::par(mar = c(5, 5, 2.6, 1.2), xpd = NA)
    metric_values <- pmax(summary_df[[metric]], y_floor)
    ylim <- c(min(metric_values, na.rm = TRUE), max(metric_values, na.rm = TRUE) * 1.10)
    graphics::plot(
      NA,
      xlim = range(trunc_levels),
      ylim = ylim,
      log = "y",
      xlab = "Truncation level",
      ylab = if (metric == "mean_abs_diff") "Mean absolute difference" else "Max absolute difference",
      cex.lab = 1.15,
      cex.axis = 1.05
    )
    graphics::axis(1, at = trunc_levels, labels = trunc_levels)
    graphics::mtext(
      if (metric == "mean_abs_diff") "Mean discrepancy by region" else "Worst-case discrepancy by region",
      side = 3,
      line = 0.6,
      adj = 0,
      font = 2,
      cex = 1.0
    )
    for (region_name in region_levels) {
      sub_df <- summary_df[summary_df$region == region_name, , drop = FALSE]
      sub_df <- sub_df[order(sub_df$truncation_level), , drop = FALSE]
      graphics::lines(
        sub_df$truncation_level,
        pmax(sub_df[[metric]], y_floor),
        type = "b",
        lwd = 2,
        pch = region_pch[[region_name]],
        col = region_cols[[region_name]]
      )
    }
  }

  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::legend(
    "center",
    legend = region_labels,
    col = region_cols[region_levels],
    pch = region_pch[region_levels],
    lwd = 2,
    cex = 1.0,
    bty = "n",
    y.intersp = 1.2
  )
}

save_surrogate_error_boxplot <- function(pointwise_df, file) {
  region_levels <- c("interior_stable", "transition", "low_separation", "high_dependence", "high_t_low_q")
  region_labels <- c("Interior", "Transition", "Low\nseparation", "High\ndependence", "High t,\nlow q")
  method_levels <- c("spline", "chebyshev")
  method_labels <- c(spline = "Spline", chebyshev = "Chebyshev")
  method_cols <- c(spline = "#D95F02", chebyshev = "#1B9E77")

  pointwise_df$region <- factor(pointwise_df$region, levels = region_levels)
  pointwise_df$method <- factor(pointwise_df$method, levels = method_levels)

  centers <- seq_along(region_levels)
  offset <- 0.18
  positions <- c(rbind(centers - offset, centers + offset))
  split_values <- unlist(
    lapply(region_levels, function(region_name) {
      lapply(method_levels, function(method_name) {
        pointwise_df$value_error[pointwise_df$region == region_name & pointwise_df$method == method_name]
      })
    }),
    recursive = FALSE
  )
  split_cols <- rep(unname(method_cols[method_levels]), times = length(region_levels))

  grDevices::png(file, width = 1450, height = 760, res = 140)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(7, 5, 2.5, 6), xpd = NA)
  graphics::boxplot(
    split_values,
    at = positions,
    xaxt = "n",
    log = "y",
    col = split_cols,
    border = "grey20",
    boxwex = 0.28,
    outline = FALSE,
    ylab = "Absolute local value error",
    xlab = "Transformed region",
    cex.lab = 1.15,
    cex.axis = 1.05
  )
  graphics::axis(1, at = centers, labels = region_labels, tick = FALSE, cex.axis = 1.0)
  graphics::abline(v = seq(1.5, length(region_levels) - 0.5, by = 1), col = "grey90", lty = 3)
  graphics::legend(
    "topleft",
    legend = unname(method_labels[method_levels]),
    fill = unname(method_cols[method_levels]),
    border = "grey20",
    bty = "n",
    cex = 0.95
  )
}

save_reduced_contour_plot <- function(surface_df, optima_df, file) {
  delta_vals <- sort(unique(surface_df$delta))
  b_vals <- sort(unique(surface_df$b))
  z_mat <- tapply(surface_df$reference, list(surface_df$b, surface_df$delta), identity)
  z_centered <- z_mat - max(z_mat, na.rm = TRUE)

  point_pch <- c(reference = 21, spline = 24, chebyshev = 22)
  point_bg <- c(
    reference = grDevices::adjustcolor("#111111", alpha.f = 0.90),
    spline = grDevices::adjustcolor("#8C8C8C", alpha.f = 0.85),
    chebyshev = grDevices::adjustcolor("#FFFFFF", alpha.f = 0.95)
  )

  plot_matrix <- t(z_centered)
  z_limits <- range(plot_matrix, finite = TRUE)
  contour_levels <- pretty(z_limits, n = 9)
  palette_fun <- grDevices::colorRampPalette(c("#F7FBFF", "#9ECAE1", "#3182BD", "#08519C"))

  grDevices::png(file, width = 1500, height = 860, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(1, 0.34))
  graphics::par(mar = c(5.2, 5.2, 2.8, 1.2), xpd = NA)

  graphics::image(
    x = delta_vals,
    y = b_vals,
    z = plot_matrix,
    col = palette_fun(72),
    zlim = z_limits,
    xlab = expression(delta),
    ylab = "b",
    xlim = c(min(delta_vals), max(delta_vals)),
    ylim = c(min(b_vals), max(b_vals)),
    useRaster = TRUE,
    cex.lab = 1.2,
    cex.axis = 1.08
  )
  graphics::contour(
    x = delta_vals,
    y = b_vals,
    z = plot_matrix,
    levels = contour_levels,
    add = TRUE,
    drawlabels = FALSE,
    col = "grey25",
    lwd = 1.3
  )

  point_order <- c("chebyshev", "spline", "reference")
  optima_df <- optima_df[match(point_order, optima_df$method), , drop = FALSE]
  for (k in seq_len(nrow(optima_df))) {
    method_name <- optima_df$method[k]
    graphics::points(
      optima_df$delta[k],
      optima_df$b[k],
      pch = point_pch[[method_name]],
      bg = point_bg[[method_name]],
      col = "black",
      lwd = 1.2,
      cex = 1.35
    )
  }
  graphics::mtext(
    "Centered reduced criterion surface with the three reported optima",
    side = 3,
    line = 0.8,
    adj = 0,
    font = 2,
    cex = 1.0
  )
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::legend(
    "center",
    legend = c("Reference", "Spline", "Chebyshev"),
    pch = c(point_pch[["reference"]], point_pch[["spline"]], point_pch[["chebyshev"]]),
    pt.bg = c(point_bg[["reference"]], point_bg[["spline"]], point_bg[["chebyshev"]]),
    pt.cex = 1.45,
    cex = 1.00,
    lwd = 1.1,
    bty = "n",
    y.intersp = 1.35
  )
}

save_runtime_scaling_plot <- function(scaling_df, file, raw_df = NULL, reference_threshold = 0.5) {
  methods <- unique(scaling_df$method)
  cols <- c(reference = "#111111", spline = "#d95f02", chebyshev = "#1b9e77")
  pch_map <- c(reference = 16, spline = 17, chebyshev = 15)
  stress_fill <- grDevices::adjustcolor("grey70", alpha.f = 0.25)

  draw_error_bars <- function(x, lower, upper, width, col) {
    graphics::segments(x, lower, x, upper, col = col, lwd = 1.1)
    graphics::segments(x - width, lower, x + width, lower, col = col, lwd = 1.1)
    graphics::segments(x - width, upper, x + width, upper, col = col, lwd = 1.1)
  }

  shade_stress_regions <- function(x_vals, y_max, stress_mask) {
    if (!any(stress_mask)) {
      return(invisible(NULL))
    }
    unique_x <- sort(unique(x_vals))
    if (length(unique_x) == 1L) {
      dx <- 0.5
    } else {
      dx <- min(diff(unique_x)) / 2
    }
    for (x0 in unique_x[stress_mask]) {
      graphics::rect(
        xleft = x0 - dx,
        xright = x0 + dx,
        ybottom = 0,
        ytop = y_max,
        col = stress_fill,
        border = NA
      )
    }
  }

  time_iqr <- NULL
  speed_iqr <- NULL
  if (!is.null(raw_df)) {
    time_iqr <- do.call(
      rbind,
      lapply(split(raw_df, list(raw_df$n, raw_df$method), drop = TRUE), function(sub_df) {
        data.frame(
          n = sub_df$n[1],
          method = sub_df$method[1],
          q25_time = as.numeric(stats::quantile(sub_df$mean_time, 0.25, names = FALSE, type = 7)),
          q75_time = as.numeric(stats::quantile(sub_df$mean_time, 0.75, names = FALSE, type = 7)),
          stringsAsFactors = FALSE
        )
      })
    )
    rownames(time_iqr) <- NULL

    ref_df <- raw_df[raw_df$method == "reference", c("n", "replicate", "mean_time")]
    names(ref_df)[3] <- "reference_time"
    ref_conv_df <- raw_df[raw_df$method == "reference", c("n", "replicate", "converged")]
    names(ref_conv_df)[3] <- "reference_converged"
    method_df <- raw_df[raw_df$method != "reference", c("n", "replicate", "method", "mean_time", "converged")]
    names(method_df)[5] <- "method_converged"
    speed_raw <- merge(method_df, ref_df, by = c("n", "replicate"), all.x = TRUE, all.y = FALSE)
    speed_raw <- merge(speed_raw, ref_conv_df, by = c("n", "replicate"), all.x = TRUE, all.y = FALSE)
    speed_raw <- speed_raw[speed_raw$reference_converged & speed_raw$method_converged, , drop = FALSE]
    speed_raw$speedup <- speed_raw$reference_time / speed_raw$mean_time
    speed_iqr <- do.call(
      rbind,
      lapply(split(speed_raw, list(speed_raw$n, speed_raw$method), drop = TRUE), function(sub_df) {
        data.frame(
          n = sub_df$n[1],
          method = sub_df$method[1],
          q25_speedup = as.numeric(stats::quantile(sub_df$speedup, 0.25, names = FALSE, type = 7)),
          q75_speedup = as.numeric(stats::quantile(sub_df$speedup, 0.75, names = FALSE, type = 7)),
          stringsAsFactors = FALSE
        )
      })
    )
    rownames(speed_iqr) <- NULL
  }

  grDevices::png(file, width = 1500, height = 680, res = 140)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 0.35))
  cap_width <- 0.04 * diff(range(scaling_df$n))
  reference_df <- scaling_df[scaling_df$method == "reference", c("n", "convergence_rate"), drop = FALSE]
  stress_mask <- reference_df$convergence_rate < reference_threshold

  left_upper <- max(c(
    scaling_df$mean_time,
    if (!is.null(time_iqr)) time_iqr$q75_time else NA_real_
  ), na.rm = TRUE) * 1.08
  graphics::par(mar = c(5, 5, 2, 1))
  graphics::plot(
    NA,
    xlim = range(scaling_df$n),
    ylim = c(0, left_upper),
    xlab = "n",
    ylab = "Mean optimization time (s)",
    cex.lab = 1.15,
    cex.axis = 1.05
  )
  shade_stress_regions(reference_df$n, left_upper, stress_mask)
  graphics::mtext("A. Mean optimization time", side = 3, line = 0.6, adj = 0, font = 2, cex = 1.0)
  for (method in methods) {
    sub_df <- scaling_df[scaling_df$method == method, ]
    if (!is.null(time_iqr)) {
      sub_iqr <- time_iqr[time_iqr$method == method, , drop = FALSE]
      sub_iqr <- sub_iqr[match(sub_df$n, sub_iqr$n), , drop = FALSE]
      draw_error_bars(sub_df$n, sub_iqr$q25_time, sub_iqr$q75_time, width = cap_width, col = cols[[method]])
    }
    graphics::lines(sub_df$n, sub_df$mean_time, type = "b", lwd = 2.2, pch = pch_map[[method]], col = cols[[method]])
  }

  right_upper <- max(c(
    scaling_df$matched_speedup_vs_reference,
    if (!is.null(speed_iqr)) speed_iqr$q75_speedup else NA_real_
  ), na.rm = TRUE) * 1.08
  graphics::par(mar = c(5, 5, 2, 1))
  graphics::plot(
    NA,
    xlim = range(scaling_df$n),
    ylim = c(0, right_upper),
    xlab = "n",
    ylab = "Matched-convergence speedup",
    cex.lab = 1.15,
    cex.axis = 1.05
  )
  shade_stress_regions(reference_df$n, right_upper, stress_mask)
  graphics::mtext("B. Matched-convergence speedup", side = 3, line = 0.6, adj = 0, font = 2, cex = 1.0)
  for (method in setdiff(methods, "reference")) {
    sub_df <- scaling_df[scaling_df$method == method, ]
    valid_speedup <- reference_df$convergence_rate[match(sub_df$n, reference_df$n)] >= reference_threshold &
      !is.na(sub_df$matched_speedup_vs_reference)
    if (!is.null(speed_iqr)) {
      sub_iqr <- speed_iqr[speed_iqr$method == method, , drop = FALSE]
      sub_iqr <- sub_iqr[match(sub_df$n, sub_iqr$n), , drop = FALSE]
      if (any(valid_speedup)) {
        draw_error_bars(
          sub_df$n[valid_speedup],
          sub_iqr$q25_speedup[valid_speedup],
          sub_iqr$q75_speedup[valid_speedup],
          width = cap_width,
          col = cols[[method]]
        )
      }
    }
    if (any(valid_speedup)) {
      graphics::lines(
        sub_df$n[valid_speedup],
        sub_df$matched_speedup_vs_reference[valid_speedup],
        type = "b",
        lwd = 2.2,
        pch = pch_map[[method]],
        col = cols[[method]]
      )
    }
  }
  graphics::abline(h = 1, lty = 2, col = "grey40")
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::legend(
    "center",
    legend = c("Reference", "Spline", "Chebyshev"),
    col = unname(cols[c("reference", "spline", "chebyshev")]),
    pch = unname(pch_map[c("reference", "spline", "chebyshev")]),
    lwd = 2.2,
    pt.cex = 1.2,
    bty = "n",
    cex = 1.0
  )
}
