script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

box <- list(
  x = c(-4, 4),
  q = c(0, 25),
  t = c(0, 2.5)
)
nu <- 6
deepest_level <- 900L
evaluation_grid <- sample_transformed_points_by_region(
  box = box,
  n_per_region = 120L,
  seed = 20260412L
)
primary_grid <- evaluation_grid[evaluation_grid$region != "high_t_low_q", , drop = FALSE]

deepest_values <- reference_transformed_clayton(
  primary_grid$x,
  primary_grid$q,
  primary_grid$t,
  nu = nu,
  log = TRUE,
  max_m = deepest_level,
  max_n = deepest_level
)

reference_levels <- c(300L, 400L, 500L, 700L, deepest_level)
reference_refinement <- do.call(
  rbind,
  lapply(reference_levels, function(level) {
    elapsed <- system.time({
      values <- reference_transformed_clayton(
        primary_grid$x,
        primary_grid$q,
        primary_grid$t,
        nu = nu,
        log = TRUE,
        max_m = level,
        max_n = level
      )
    })[["elapsed"]]
    diff <- abs(values - deepest_values)
    data.frame(
      truncation_level = level,
      n = nrow(primary_grid),
      mean_abs_diff = mean(diff),
      median_abs_diff = median(diff),
      max_abs_diff = max(diff),
      share_gt_1e6 = mean(diff > 1e-6),
      elapsed_seconds = elapsed,
      region_scope = "primary_benchmark",
      stringsAsFactors = FALSE
    )
  })
)

sim <- simulate_clayton_beta_example(
  n = 120L,
  beta = 0,
  delta = 6,
  b = 0.35,
  nu = nu,
  seed = 20260409
)
data <- prepare_reduced_clayton_data(sim$y, sim$coords, sim$X, m = 2L)
reference_opt <- optimize_reduced_wpl(
  objective = wpl_reference_reduced,
  start = c(delta = 5.5, b = 0.30),
  lower = c(1.0, 0.05),
  upper = c(20.0, 1.5),
  data = data,
  fixed_beta = sim$true$beta,
  fixed_nu = sim$true$nu
)
reference_value_at_reference_opt <- wpl_reference_reduced(
  reference_opt$par,
  data = data,
  fixed_beta = sim$true$beta,
  fixed_nu = sim$true$nu
)

config_specs <- list(
  coarse = list(
    reference_cap = 400L,
    table_grid = list(x = 13L, q = 11L, t = 11L),
    spline_knots = list(x = 9L, q = 7L, t = 7L, degree = 3L),
    spline_penalty = list(x = 1e-5, q = 1e-5, t = 1e-5),
    cheb_degree = list(x = 6L, q = 6L, t = 6L)
  ),
  baseline = list(
    reference_cap = 700L,
    table_grid = list(x = 17L, q = 13L, t = 13L),
    spline_knots = list(x = 11L, q = 9L, t = 9L, degree = 3L),
    spline_penalty = list(x = 1e-6, q = 1e-6, t = 1e-6),
    cheb_degree = list(x = 8L, q = 8L, t = 8L)
  ),
  fine = list(
    reference_cap = 900L,
    table_grid = list(x = 21L, q = 15L, t = 15L),
    spline_knots = list(x = 13L, q = 11L, t = 11L, degree = 3L),
    spline_penalty = list(x = 1e-7, q = 1e-7, t = 1e-7),
    cheb_degree = list(x = 10L, q = 10L, t = 10L)
  )
)

summarize_config <- function(config_name, method_name, cap, table_points, table_time, fit_time, fit_obj, opt_obj) {
  local_values <- evaluate_transformed_surrogate(
    fit_obj,
    primary_grid$x,
    primary_grid$q,
    primary_grid$t
  )
  local_error <- abs(local_values - deepest_values)
  reference_value_at_current_optimum <- wpl_reference_reduced(
    opt_obj$opt$par,
    data = data,
    fixed_beta = sim$true$beta,
    fixed_nu = sim$true$nu
  )
  data.frame(
    config = config_name,
    method = method_name,
    reference_cap = cap,
    table_points = table_points,
    table_build_seconds = table_time,
    fit_seconds = fit_time,
    total_offline_seconds = table_time + fit_time,
    mean_value_error = mean(local_error),
    median_value_error = median(local_error),
    max_value_error = max(local_error),
    optimization_time_seconds = opt_obj$elapsed,
    converged = identical(opt_obj$opt$convergence, 0L),
    delta_gap_vs_reference = unname(opt_obj$opt$par[1] - reference_opt$par[1]),
    b_gap_vs_reference = unname(opt_obj$opt$par[2] - reference_opt$par[2]),
    reference_criterion_drop = reference_value_at_reference_opt - reference_value_at_current_optimum,
    stringsAsFactors = FALSE
  )
}

surrogate_refinement <- do.call(
  rbind,
  lapply(names(config_specs), function(config_name) {
    spec <- config_specs[[config_name]]

    spline_table_time <- system.time({
      local_table <- build_clayton_local_table(
        box = box,
        grid = spec$table_grid,
        nu = nu,
        backend = "series_reference",
        tol = 1e-12,
        max_m = spec$reference_cap,
        max_n = spec$reference_cap
      )
    })[["elapsed"]]
    spline_fit_time <- system.time({
      spline_fit <- fit_spline_surrogate(
        local_table = local_table,
        box = box,
        knots = spec$spline_knots,
        penalty = spec$spline_penalty,
        nu = nu
      )
    })[["elapsed"]]
    spline_opt <- system.time({
      opt <- optimize_reduced_wpl(
        objective = wpl_surrogate_reduced,
        start = c(delta = 5.5, b = 0.30),
        lower = c(1.0, 0.05),
        upper = c(20.0, 1.5),
        data = data,
        surrogate_fit = spline_fit,
        fixed_beta = sim$true$beta,
        fixed_nu = sim$true$nu
      )
    })[["elapsed"]]
    spline_row <- summarize_config(
      config_name = config_name,
      method_name = "spline",
      cap = spec$reference_cap,
      table_points = nrow(local_table),
      table_time = spline_table_time,
      fit_time = spline_fit_time,
      fit_obj = spline_fit,
      opt_obj = list(opt = opt, elapsed = spline_opt)
    )

    cheb_table_time <- system.time({
      cheb_table <- expand.grid(
        x = chebyshev_nodes(spec$cheb_degree$x, box$x),
        q = chebyshev_nodes(spec$cheb_degree$q, box$q),
        t = chebyshev_nodes(spec$cheb_degree$t, box$t),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      cheb_table$value <- reference_transformed_clayton(
        cheb_table$x,
        cheb_table$q,
        cheb_table$t,
        nu = nu,
        log = TRUE,
        max_m = spec$reference_cap,
        max_n = spec$reference_cap
      )
      cheb_table$nu <- nu
    })[["elapsed"]]
    cheb_fit_time <- system.time({
      cheb_fit <- fit_chebyshev_surrogate(
        local_table = cheb_table,
        box = box,
        degree = spec$cheb_degree,
        nu = nu
      )
    })[["elapsed"]]
    cheb_opt <- system.time({
      opt <- optimize_reduced_wpl(
        objective = wpl_surrogate_reduced,
        start = c(delta = 5.5, b = 0.30),
        lower = c(1.0, 0.05),
        upper = c(20.0, 1.5),
        data = data,
        surrogate_fit = cheb_fit,
        fixed_beta = sim$true$beta,
        fixed_nu = sim$true$nu
      )
    })[["elapsed"]]
    cheb_row <- summarize_config(
      config_name = config_name,
      method_name = "chebyshev",
      cap = spec$reference_cap,
      table_points = nrow(cheb_table),
      table_time = cheb_table_time,
      fit_time = cheb_fit_time,
      fit_obj = cheb_fit,
      opt_obj = list(opt = opt, elapsed = cheb_opt)
    )

    rbind(spline_row, cheb_row)
  })
)
surrogate_refinement$config <- factor(surrogate_refinement$config, levels = names(config_specs))
surrogate_refinement$method <- factor(surrogate_refinement$method, levels = c("spline", "chebyshev"))
surrogate_refinement <- surrogate_refinement[order(surrogate_refinement$config, surrogate_refinement$method), , drop = FALSE]
surrogate_refinement$config <- as.character(surrogate_refinement$config)
surrogate_refinement$method <- as.character(surrogate_refinement$method)

save_refinement_plot <- function(reference_df, surrogate_df, file) {
  config_levels <- c("coarse", "baseline", "fine")
  config_x <- seq_along(config_levels)
  method_cols <- c(spline = "#d95f02", chebyshev = "#1b9e77")
  method_pch <- c(spline = 17, chebyshev = 15)
  y_floor <- 1e-15

  grDevices::png(file, width = 1200, height = 600, res = 120)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))

  graphics::plot(
    reference_df$truncation_level,
    pmax(reference_df$mean_abs_diff, 1e-15),
    type = "b",
    log = "y",
    ylim = range(c(reference_df$mean_abs_diff, reference_df$max_abs_diff, 1e-15), finite = TRUE),
    xlab = "Truncation level",
    ylab = "Absolute difference versus level 900",
    pch = 16,
    lwd = 2,
    col = "#3182bd",
    main = "Reference refinement on the primary region"
  )
  graphics::lines(
    reference_df$truncation_level,
    pmax(reference_df$max_abs_diff, 1e-15),
    type = "b",
    pch = 17,
    lwd = 2,
    col = "#08519c"
  )
  graphics::legend(
    "topright",
    legend = c("Mean absolute difference", "Max absolute difference"),
    col = c("#3182bd", "#08519c"),
    pch = c(16, 17),
    lwd = 2,
    bty = "n"
  )

  graphics::plot(
    NA,
    xlim = c(0.8, 3.2),
    ylim = range(pmax(surrogate_df$mean_value_error, y_floor), finite = TRUE),
    log = "y",
    xlab = "Surrogate resolution",
    ylab = "Mean local value error",
    xaxt = "n",
    main = "Surrogate refinement on the primary region"
  )
  graphics::axis(1, at = config_x, labels = c("Coarse", "Baseline", "Fine"))
  for (method_name in names(method_cols)) {
    sub_df <- surrogate_df[surrogate_df$method == method_name, , drop = FALSE]
      graphics::lines(
        config_x,
        pmax(sub_df$mean_value_error, y_floor),
        type = "b",
        pch = method_pch[[method_name]],
        lwd = 2,
      col = method_cols[[method_name]]
    )
  }
  graphics::legend(
    "topright",
    legend = c("Spline", "Chebyshev"),
    col = method_cols[c("spline", "chebyshev")],
    pch = method_pch[c("spline", "chebyshev")],
    lwd = 2,
    bty = "n"
  )
}

write.csv(reference_refinement, result_path("reference_refinement_summary.csv"), row.names = FALSE)
write.csv(surrogate_refinement, result_path("surrogate_refinement_summary.csv"), row.names = FALSE)

save_refinement_plot(
  reference_refinement,
  surrogate_refinement,
  result_path("fig_refinement_diagnostics.png")
)

