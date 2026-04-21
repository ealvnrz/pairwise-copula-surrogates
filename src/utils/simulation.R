coerce_beta_design <- function(beta, n) {
  beta <- as.numeric(beta)
  if (length(beta) < 1L) {
    stop("`beta` must contain at least one coefficient.", call. = FALSE)
  }

  if (length(beta) == 1L) {
    X <- matrix(1, nrow = n, ncol = 1L)
    colnames(X) <- "intercept"
    return(X)
  }

  X <- cbind(1, replicate(length(beta) - 1L, runif(n)))
  colnames(X) <- c("intercept", paste0("x", seq_len(length(beta) - 1L)))
  X
}

beta_parameter_names <- function(p) {
  if (p < 1L) {
    stop("`p` must be positive.", call. = FALSE)
  }
  c("mean", if (p > 1L) paste0("mean", seq_len(p - 1L)))
}

clip_unit_response <- function(y, eps = 1e-6) {
  pmin(pmax(as.numeric(y), eps), 1 - eps)
}

simulate_clayton_beta_example <- function(
    n = 150L,
    beta = 0,
    delta = 6,
    b = 0.35,
    nu = 6,
    seed = 20260409) {
  if (!requireNamespace("GeoModels", quietly = TRUE)) {
    stop("Package `GeoModels` is required for simulation.", call. = FALSE)
  }

  set.seed(seed)
  coords <- cbind(runif(n), runif(n))
  X <- coerce_beta_design(beta, n)
  beta <- as.numeric(beta)

  geo_param <- c(
    list(
      smooth = 0,
      power2 = 4,
      min = 0,
      max = 1,
      scale = b,
      nugget = 0,
      shape = delta,
      nu = nu
    ),
    stats::setNames(as.list(beta), beta_parameter_names(length(beta)))
  )

  sim <- GeoModels::GeoSimCopula(
    coordx = coords,
    corrmodel = "GenWend",
    model = "Beta2",
    param = geo_param,
    copula = "Clayton",
    sparse = TRUE,
    X = X
  )

  list(
    y = clip_unit_response(sim$data),
    coords = coords,
    X = X,
    true = list(beta = beta, delta = delta, b = b, nu = nu)
  )
}

simulate_spatial_factorcopula_example <- function(
    n = 100L,
    R = 60L,
    alpha = 1.2,
    b = 0.30,
    factor_spec = NULL,
    seed = 20260501L) {
  if (!is.finite(n) || n < 2L) {
    stop("`n` must be at least 2.", call. = FALSE)
  }
  if (!is.finite(R) || R < 5L) {
    stop("`R` must be at least 5.", call. = FALSE)
  }
  factor_spec <- coerce_factorcopula_spec(factor_spec)

  set.seed(seed)
  coords <- cbind(runif(n), runif(n))
  dmat <- as.matrix(stats::dist(coords))
  sigma <- spatial_factor_correlation(dmat, alpha = alpha, b = b)
  diag(sigma) <- 1
  chol_factor <- chol(sigma)
  latent <- matrix(stats::rnorm(R * n), nrow = R) %*% chol_factor
  common_factor <- factorcopula_rpareto(R, factor_spec)
  Y <- latent + common_factor
  u_exact <- matrix(
    factorcopula_fw(as.vector(Y), factor_spec),
    nrow = R,
    ncol = n
  )

  list(
    Y = Y,
    coords = coords,
    u_exact = u_exact,
    common_factor = common_factor,
    true = list(alpha = alpha, b = b, factor_spec = factor_spec, R = R)
  )
}
