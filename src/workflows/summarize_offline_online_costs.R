script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

reference_bundle <- readRDS(result_path("01_exact_local_grid.rds"))
practical_path <- result_path("practical_beta_free_summary.csv")

if (!file.exists(practical_path)) {
  source(file.path(script_dir, "09_practical_beta_free_benchmark.R"), local = FALSE)
}

practical_summary <- read.csv(practical_path, stringsAsFactors = FALSE)

reference_grid_seconds <- system.time({
  grid_rebuild <- build_clayton_local_table(
    box = reference_bundle$primary$box,
    grid = reference_bundle$primary$grid,
    nu = reference_bundle$primary$nu,
    backend = "series_reference",
    tol = 1e-12,
    max_m = 700L,
    max_n = 700L
  )
})[["elapsed"]]

spline_fit_seconds <- system.time({
  fit_spline_surrogate(
    local_table = grid_rebuild,
    box = reference_bundle$primary$box,
    knots = list(x = 11L, q = 9L, t = 9L, degree = 3L),
    penalty = list(x = 1e-6, q = 1e-6, t = 1e-6),
    nu = reference_bundle$primary$nu
  )
})[["elapsed"]]

cheb_table_seconds <- system.time({
  cheb_table <- expand.grid(
    x = chebyshev_nodes(8L, reference_bundle$primary$box$x),
    q = chebyshev_nodes(8L, reference_bundle$primary$box$q),
    t = chebyshev_nodes(8L, reference_bundle$primary$box$t),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  cheb_table$value <- reference_transformed_clayton(
    cheb_table$x,
    cheb_table$q,
    cheb_table$t,
    nu = reference_bundle$primary$nu,
    log = TRUE,
    max_m = 700L,
    max_n = 700L
  )
  cheb_table$nu <- reference_bundle$primary$nu
})[["elapsed"]]

cheb_fit_seconds <- system.time({
  fit_chebyshev_surrogate(
    local_table = cheb_table,
    box = reference_bundle$primary$box,
    degree = list(x = 8L, q = 8L, t = 8L),
    nu = reference_bundle$primary$nu
  )
})[["elapsed"]]

offline_online <- data.frame(
  method = c("reference", "spline", "chebyshev"),
  table_build_seconds = c(0, reference_grid_seconds, cheb_table_seconds),
  fit_seconds = c(0, spline_fit_seconds, cheb_fit_seconds),
  total_offline_seconds = c(0, reference_grid_seconds + spline_fit_seconds, cheb_table_seconds + cheb_fit_seconds),
  practical_online_seconds = c(
    practical_summary$optimization_time[match("reference", practical_summary$method)],
    practical_summary$optimization_time[match("spline", practical_summary$method)],
    practical_summary$optimization_time[match("chebyshev", practical_summary$method)]
  ),
  stringsAsFactors = FALSE
)

compute_break_even <- function(context_name, online_column) {
  ref_time <- offline_online[[online_column]][offline_online$method == "reference"]
  do.call(
    rbind,
    lapply(setdiff(offline_online$method, "reference"), function(method_name) {
      method_time <- offline_online[[online_column]][offline_online$method == method_name]
      offline_time <- offline_online$total_offline_seconds[offline_online$method == method_name]
      savings <- ref_time - method_time
      data.frame(
        context = context_name,
        method = method_name,
        offline_seconds = offline_time,
        online_seconds = method_time,
        per_run_saving_seconds = savings,
        break_even_runs = if (savings > 0) ceiling(offline_time / savings) else Inf,
        stringsAsFactors = FALSE
      )
    })
  )
}

break_even <- compute_break_even("practical_simulation_n120", "practical_online_seconds")

write.csv(offline_online, result_path("offline_online_cost_summary.csv"), row.names = FALSE)
write.csv(break_even, result_path("break_even_summary.csv"), row.names = FALSE)

