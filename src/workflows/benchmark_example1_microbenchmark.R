script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

spline_fit <- readRDS(result_path("02_spline_fit.rds"))
cheb_bundle <- readRDS(result_path("03_chebyshev_fit.rds"))
comparison <- readRDS(result_path("05_reduced_wpl_comparison.rds"))

sim <- comparison$simulation
data <- prepare_reduced_clayton_data(sim$y, sim$coords, sim$X, m = 2L)

benchmark <- runtime_benchmark_summary(
  data = data,
  truth = sim$true,
  spline_fit = spline_fit,
  chebyshev_fit = cheb_bundle$fit,
  optimization_times = 5L
)

benchmark_summary <- rbind(
  cbind(component = "local", method = rownames(benchmark$local), as.data.frame(benchmark$local), row.names = NULL),
  cbind(component = "criterion", method = rownames(benchmark$criterion), as.data.frame(benchmark$criterion), row.names = NULL),
  cbind(component = "optimization", method = rownames(benchmark$optimization), as.data.frame(benchmark$optimization), row.names = NULL)
)
write.csv(benchmark_summary, result_path("example1_microbenchmark_summary.csv"), row.names = FALSE)

