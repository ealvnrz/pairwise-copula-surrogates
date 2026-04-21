script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

nu <- 6
box <- list(
  x = c(-4, 4),
  q = c(0, 25),
  t = c(0, 2.5)
)
grid <- list(x = 17, q = 13, t = 13)

local_table <- build_clayton_local_table(
  box = box,
  grid = grid,
  nu = nu,
  backend = "series_reference",
  tol = 1e-12,
  max_m = 700L,
  max_n = 700L
)

saveRDS(
  list(
    primary = list(
      box = box,
      grid = grid,
      nu = nu,
      local_table = local_table
    ),
    box = box,
    grid = grid,
    nu = nu,
    local_table = local_table
  ),
  result_path("01_exact_local_grid.rds")
)
write.csv(local_table, result_path("01_exact_local_grid.csv"), row.names = FALSE)

