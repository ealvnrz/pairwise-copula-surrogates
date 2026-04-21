script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

spline_fit <- readRDS(result_path("02_spline_fit.rds"))
cheb_bundle <- readRDS(result_path("03_chebyshev_fit.rds"))

benchmark_rows <- read.csv(result_path("parameter_gap_replicates.csv"), stringsAsFactors = FALSE)
subset_ids <- seq_len(10L)
grad_step <- 0.02
hess_step <- 0.04

format_num <- function(x, digits = 4L) {
  ifelse(is.na(x), "---", formatC(x, format = "f", digits = digits))
}

quantile95 <- function(x) {
  stats::quantile(x, probs = 0.95, names = FALSE, na.rm = TRUE, type = 7)
}

compute_discrepancy <- function(rep_id, method_name) {
  ref_row <- benchmark_rows[benchmark_rows$replicate == rep_id & benchmark_rows$method == "reference", , drop = FALSE]
  method_row <- benchmark_rows[benchmark_rows$replicate == rep_id & benchmark_rows$method == method_name, , drop = FALSE]

  if (!nrow(ref_row) || !nrow(method_row)) {
    return(NULL)
  }

  seed <- 20260409L + rep_id
  sim <- simulate_clayton_beta_example(
    n = 120L,
    beta = 0,
    delta = 6,
    b = 0.35,
    nu = 6,
    seed = seed
  )
  data <- prepare_reduced_clayton_data(sim$y, sim$coords, sim$X, m = 2L)
  total_weight <- sum(data$design$weights)

  objective_ref <- function(theta) {
    wpl_reference_reduced(theta, data = data, fixed_beta = sim$true$beta, fixed_nu = sim$true$nu) / total_weight
  }
  objective_sur <- function(theta) {
    wpl_surrogate_reduced(theta, data = data, surrogate_fit = if (identical(method_name, "spline")) spline_fit else cheb_bundle$fit, fixed_beta = sim$true$beta, fixed_nu = sim$true$nu) / total_weight
  }

  theta_ref <- c(ref_row$delta, ref_row$b)
  theta_sur <- c(method_row$delta, method_row$b)
  reference_converged <- isTRUE(ref_row$converged)
  method_converged <- isTRUE(method_row$converged)

  if (!reference_converged || !method_converged) {
    return(data.frame(
      replicate = rep_id,
      method = method_name,
      reference_converged = reference_converged,
      method_converged = method_converged,
      score_gap_at_reference = NA_real_,
      hessian_gap_at_reference = NA_real_,
      relative_hessian_gap_at_reference = NA_real_,
      score_gap_at_surrogate = NA_real_,
      hessian_gap_at_surrogate = NA_real_,
      relative_hessian_gap_at_surrogate = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  grad_ref_ref <- finite_difference_gradient(objective_ref, theta_ref, step = grad_step)
  grad_sur_ref <- finite_difference_gradient(objective_sur, theta_ref, step = grad_step)
  hess_ref_ref <- finite_difference_hessian(objective_ref, theta_ref, step = hess_step)
  hess_sur_ref <- finite_difference_hessian(objective_sur, theta_ref, step = hess_step)

  grad_ref_sur <- finite_difference_gradient(objective_ref, theta_sur, step = grad_step)
  grad_sur_sur <- finite_difference_gradient(objective_sur, theta_sur, step = grad_step)
  hess_ref_sur <- finite_difference_hessian(objective_ref, theta_sur, step = hess_step)
  hess_sur_sur <- finite_difference_hessian(objective_sur, theta_sur, step = hess_step)

  data.frame(
    replicate = rep_id,
    method = method_name,
    reference_converged = reference_converged,
    method_converged = method_converged,
    score_gap_at_reference = sqrt(sum((grad_sur_ref - grad_ref_ref)^2)),
    hessian_gap_at_reference = sqrt(sum((hess_sur_ref - hess_ref_ref)^2)),
    relative_hessian_gap_at_reference = sqrt(sum((hess_sur_ref - hess_ref_ref)^2)) / pmax(sqrt(sum(hess_ref_ref^2)), 1e-12),
    score_gap_at_surrogate = sqrt(sum((grad_sur_sur - grad_ref_sur)^2)),
    hessian_gap_at_surrogate = sqrt(sum((hess_sur_sur - hess_ref_sur)^2)),
    relative_hessian_gap_at_surrogate = sqrt(sum((hess_sur_sur - hess_ref_sur)^2)) / pmax(sqrt(sum(hess_ref_sur^2)), 1e-12),
    stringsAsFactors = FALSE
  )
}

replicate_df <- do.call(
  rbind,
  unlist(
    lapply(subset_ids, function(rep_id) {
      lapply(c("spline", "chebyshev"), function(method_name) {
        compute_discrepancy(rep_id, method_name)
      })
    }),
    recursive = FALSE
  )
)

summary_df <- do.call(
  rbind,
  lapply(split(replicate_df, replicate_df$method), function(sub_df) {
    used_df <- sub_df[sub_df$reference_converged & sub_df$method_converged, , drop = FALSE]
    data.frame(
      method = sub_df$method[1],
      n_used = nrow(used_df),
      mean_score_gap_at_reference = mean(used_df$score_gap_at_reference, na.rm = TRUE),
      p95_score_gap_at_reference = quantile95(used_df$score_gap_at_reference),
      mean_relative_hessian_gap_at_reference = mean(used_df$relative_hessian_gap_at_reference, na.rm = TRUE),
      p95_relative_hessian_gap_at_reference = quantile95(used_df$relative_hessian_gap_at_reference),
      mean_score_gap_at_surrogate = mean(used_df$score_gap_at_surrogate, na.rm = TRUE),
      mean_relative_hessian_gap_at_surrogate = mean(used_df$relative_hessian_gap_at_surrogate, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
rownames(summary_df) <- NULL
summary_df <- summary_df[match(c("spline", "chebyshev"), summary_df$method), , drop = FALSE]

write.csv(replicate_df, result_path("example1_criterion_derivative_replicates.csv"), row.names = FALSE)
write.csv(summary_df, result_path("example1_criterion_derivative_summary.csv"), row.names = FALSE)

