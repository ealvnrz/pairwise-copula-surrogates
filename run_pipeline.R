args <- commandArgs(trailingOnly = TRUE)
scope <- if (length(args) > 0L) args[[1]] else "all"

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
workflow_dir <- file.path(repo_root, "src", "workflows")
rscript_bin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

workflow_groups <- list(
  core = c(
    "build_exact_local_grid.R",
    "fit_spline_surrogate.R",
    "fit_chebyshev_surrogate.R",
    "check_local_accuracy.R",
    "benchmark_example1_reduced.R",
    "benchmark_example1_microbenchmark.R",
    "audit_example1_truncation_stability.R",
    "benchmark_example1_practical_fixed_nu.R",
    "benchmark_example1_scaling.R",
    "study_surrogate_refinement.R",
    "summarize_offline_online_costs.R",
    "diagnose_example1_curvature.R",
    "prepare_example2_box_and_refinement.R",
    "check_example2_local_accuracy.R",
    "benchmark_example2_practical.R",
    "benchmark_example2_scaling.R"
  ),
  all = c(
    "build_exact_local_grid.R",
    "fit_spline_surrogate.R",
    "fit_chebyshev_surrogate.R",
    "check_local_accuracy.R",
    "benchmark_example1_reduced.R",
    "benchmark_example1_microbenchmark.R",
    "audit_example1_truncation_stability.R",
    "benchmark_example1_practical_fixed_nu.R",
    "benchmark_example1_scaling.R",
    "study_surrogate_refinement.R",
    "summarize_offline_online_costs.R",
    "diagnose_example1_curvature.R",
    "prepare_example2_box_and_refinement.R",
    "check_example2_local_accuracy.R",
    "benchmark_example2_practical.R",
    "benchmark_example2_scaling.R",
    "assess_spline_penalty_sensitivity.R",
    "illustrate_example2_reuse_weights.R",
    "diagnose_example1_active_weighted_error.R",
    "stress_example1_free_nu.R",
    "diagnose_example1_uncertainty_proxy.R",
    "audit_example1_reference_depth.R",
    "diagnose_clipping.R",
    "diagnose_example1_criterion_derivatives.R"
  )
)

if (!scope %in% names(workflow_groups)) {
  stop("Unknown scope. Use `core` or `all`.", call. = FALSE)
}

for (workflow in workflow_groups[[scope]]) {
  workflow_path <- file.path(workflow_dir, workflow)
  message(sprintf("Running %s", workflow))
  status <- system2(rscript_bin, workflow_path)
  if (!identical(status, 0L)) {
    stop(sprintf("Workflow failed: %s", workflow), call. = FALSE)
  }
}

message(sprintf("Completed `%s` workflow run.", scope))
