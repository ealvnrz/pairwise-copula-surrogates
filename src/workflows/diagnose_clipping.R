script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

spline_fit <- readRDS(result_path("02_spline_fit.rds"))
factor_bundle <- readRDS(result_path("25_factorcopula_example_box_and_refinement.rds"))

fmt_num <- function(x, digits = 4L) {
  ifelse(is.na(x), "---", formatC(x, format = "f", digits = digits))
}

weighted_clip_summary <- function(x, q, t, weights, box) {
  outside_x <- x < box$x[1] | x > box$x[2]
  outside_q <- q < box$q[1] | q > box$q[2]
  outside_t <- t < box$t[1] | t > box$t[2]
  outside_any <- outside_x | outside_q | outside_t
  total_weight <- sum(weights)

  data.frame(
    clip_share_any = sum(weights[outside_any]) / total_weight,
    clip_share_x = sum(weights[outside_x]) / total_weight,
    clip_share_q = sum(weights[outside_q]) / total_weight,
    clip_share_t = sum(weights[outside_t]) / total_weight,
    n_points = length(weights),
    total_weight = total_weight,
    stringsAsFactors = FALSE
  )
}

quantile95 <- function(x) stats::quantile(x, probs = 0.95, names = FALSE, na.rm = TRUE, type = 7)

collect_example1_reduced <- function() {
  ref_rows <- read.csv(result_path("parameter_gap_replicates.csv"), stringsAsFactors = FALSE)
  ref_rows <- ref_rows[ref_rows$method == "reference", , drop = FALSE]
  do.call(
    rbind,
    lapply(seq_len(nrow(ref_rows)), function(k) {
      row <- ref_rows[k, ]
      if (!isTRUE(row$converged)) {
        return(NULL)
      }
      seed <- 20260409L + row$replicate
      sim <- simulate_clayton_beta_example(n = 120L, beta = 0, delta = 6, b = 0.35, nu = 6, seed = seed)
      data <- prepare_reduced_clayton_data(sim$y, sim$coords, sim$X, m = 2L)
      comp <- build_reduced_clayton_components(c(row$delta, row$b), data, fixed_beta = sim$true$beta, fixed_nu = sim$true$nu)
      idx <- comp$index
      transformed <- transform_clayton_inputs(comp$u[idx[, "i"]], comp$u[idx[, "j"]], comp$r)
      cbind(
        benchmark = "Example 1 reduced",
        replicate = row$replicate,
        weighted_clip_summary(transformed$x, transformed$q, transformed$t, comp$weights, spline_fit$box),
        stringsAsFactors = FALSE
      )
    })
  )
}

collect_example1_practical <- function() {
  ref_rows <- read.csv(result_path("practical_beta_free_replicates.csv"), stringsAsFactors = FALSE)
  ref_rows <- ref_rows[ref_rows$method == "reference", , drop = FALSE]
  do.call(
    rbind,
    lapply(seq_len(nrow(ref_rows)), function(k) {
      row <- ref_rows[k, ]
      if (!isTRUE(row$converged)) {
        return(NULL)
      }
      seed <- 20261000L + row$replicate
      sim <- simulate_clayton_beta_example(n = 120L, beta = c(-0.35, 0.9), delta = 6, b = 0.35, nu = 6, seed = seed)
      data <- prepare_reduced_clayton_data(sim$y, sim$coords, sim$X, m = 2L)
      theta <- c(row$beta_intercept, row$beta_x1, row$delta, row$b)
      comp <- build_practical_clayton_components(theta, data, fixed_nu = sim$true$nu)
      idx <- comp$index
      transformed <- transform_clayton_inputs(comp$u[idx[, "i"]], comp$u[idx[, "j"]], comp$r)
      cbind(
        benchmark = "Example 1 practical",
        replicate = row$replicate,
        weighted_clip_summary(transformed$x, transformed$q, transformed$t, comp$weights, spline_fit$box),
        stringsAsFactors = FALSE
      )
    })
  )
}

collect_example2_practical <- function() {
  ref_rows <- read.csv(result_path("factorcopula_example_replicates.csv"), stringsAsFactors = FALSE)
  ref_rows <- ref_rows[ref_rows$method == "reference", , drop = FALSE]
  factor_spec <- factor_bundle$factor_spec
  factor_box <- factor_bundle$box
  do.call(
    rbind,
    lapply(seq_len(nrow(ref_rows)), function(k) {
      row <- ref_rows[k, ]
      if (!isTRUE(row$converged)) {
        return(NULL)
      }
      seed <- 20268400L + row$replicate
      sim <- simulate_spatial_factorcopula_example(n = 200L, R = 60L, alpha = 1.2, b = 0.30, factor_spec = factor_spec, seed = seed)
      data <- prepare_pairwise_factor_data(sim$Y, sim$coords, m = 6L, fixed_factor = factor_spec)
      comp <- build_factorcopula_components(c(row$alpha, row$b), data, fixed_factor = factor_spec)
      cbind(
        benchmark = "Example 2 practical",
        replicate = row$replicate,
        weighted_clip_summary(comp$x_vec, comp$q_vec, comp$t_vec, comp$weights_vec, factor_box),
        stringsAsFactors = FALSE
      )
    })
  )
}

clip_replicates <- rbind(
  collect_example1_reduced(),
  collect_example1_practical(),
  collect_example2_practical()
)

clip_summary <- do.call(
  rbind,
  lapply(split(clip_replicates, clip_replicates$benchmark), function(sub_df) {
    data.frame(
      benchmark = sub_df$benchmark[1],
      n = nrow(sub_df),
      mean_clip_share_any = mean(sub_df$clip_share_any, na.rm = TRUE),
      p95_clip_share_any = quantile95(sub_df$clip_share_any),
      max_clip_share_any = max(sub_df$clip_share_any, na.rm = TRUE),
      mean_clip_share_x = mean(sub_df$clip_share_x, na.rm = TRUE),
      mean_clip_share_q = mean(sub_df$clip_share_q, na.rm = TRUE),
      mean_clip_share_t = mean(sub_df$clip_share_t, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
rownames(clip_summary) <- NULL

write.csv(clip_replicates, result_path("clipping_diagnostics_replicates.csv"), row.names = FALSE)
write.csv(clip_summary, result_path("clipping_diagnostics_summary.csv"), row.names = FALSE)

