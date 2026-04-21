script_dir <- dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "..", "setup", "project_setup.R"))

reference_bundle <- readRDS(result_path("01_exact_local_grid.rds"))
degree <- list(x = 8L, q = 8L, t = 8L)

cheb_table <- expand.grid(
  x = chebyshev_nodes(degree$x, reference_bundle$primary$box$x),
  q = chebyshev_nodes(degree$q, reference_bundle$primary$box$q),
  t = chebyshev_nodes(degree$t, reference_bundle$primary$box$t),
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

cheb_fit <- fit_chebyshev_surrogate(
  local_table = cheb_table,
  box = reference_bundle$primary$box,
  degree = degree,
  nu = reference_bundle$primary$nu
)

saveRDS(list(table = cheb_table, fit = cheb_fit), result_path("03_chebyshev_fit.rds"))

