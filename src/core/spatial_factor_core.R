default_factorcopula_spec <- function(
    pareto_shape = 2.5,
    pareto_scale = 1.0,
    quadrature_n = 40L,
    lookup_grid_size = 4000L,
    lower_tail_target = 1e-6,
    upper_tail_target = 1e-6) {
  list(
    pareto_shape = as.numeric(pareto_shape),
    pareto_scale = as.numeric(pareto_scale),
    quadrature_n = as.integer(quadrature_n),
    lookup_grid_size = as.integer(lookup_grid_size),
    lower_tail_target = as.numeric(lower_tail_target),
    upper_tail_target = as.numeric(upper_tail_target)
  )
}

coerce_factorcopula_spec <- function(fixed_factor = NULL) {
  defaults <- default_factorcopula_spec()
  if (is.null(fixed_factor)) {
    fixed_factor <- defaults
  } else if (is.list(fixed_factor)) {
    fixed_factor <- modifyList(defaults, fixed_factor)
  } else {
    stop("`fixed_factor` must be a list or NULL.", call. = FALSE)
  }
  validate_factorcopula_spec(fixed_factor)
  fixed_factor
}

validate_factorcopula_spec <- function(fixed_factor) {
  required <- c(
    "pareto_shape",
    "pareto_scale",
    "quadrature_n",
    "lookup_grid_size",
    "lower_tail_target",
    "upper_tail_target"
  )
  if (!all(required %in% names(fixed_factor))) {
    stop(sprintf("`fixed_factor` must contain: %s.", paste(required, collapse = ", ")), call. = FALSE)
  }
  if (!is.finite(fixed_factor$pareto_shape) || fixed_factor$pareto_shape <= 1) {
    stop("`pareto_shape` must be greater than 1.", call. = FALSE)
  }
  if (!is.finite(fixed_factor$pareto_scale) || fixed_factor$pareto_scale <= 0) {
    stop("`pareto_scale` must be strictly positive.", call. = FALSE)
  }
  if (!is.finite(fixed_factor$quadrature_n) || fixed_factor$quadrature_n < 8L) {
    stop("`quadrature_n` must be at least 8.", call. = FALSE)
  }
  if (!is.finite(fixed_factor$lookup_grid_size) || fixed_factor$lookup_grid_size < 200L) {
    stop("`lookup_grid_size` must be at least 200.", call. = FALSE)
  }
  if (!is.finite(fixed_factor$lower_tail_target) || fixed_factor$lower_tail_target <= 0 || fixed_factor$lower_tail_target >= 0.1) {
    stop("`lower_tail_target` must lie in (0, 0.1).", call. = FALSE)
  }
  if (!is.finite(fixed_factor$upper_tail_target) || fixed_factor$upper_tail_target <= 0 || fixed_factor$upper_tail_target >= 0.1) {
    stop("`upper_tail_target` must lie in (0, 0.1).", call. = FALSE)
  }
  invisible(TRUE)
}

spatial_factor_correlation <- function(distance, alpha, b, eps = 1e-10) {
  if (any(!is.finite(distance)) || any(distance < 0)) {
    stop("`distance` must be finite and nonnegative.", call. = FALSE)
  }
  if (!is.finite(alpha) || length(alpha) != 1L || alpha <= 0 || alpha > 2) {
    stop("`alpha` must lie in (0, 2].", call. = FALSE)
  }
  if (!is.finite(b) || length(b) != 1L || b <= 0) {
    stop("`b` must be strictly positive.", call. = FALSE)
  }
  pmin(exp(- (distance / b)^alpha), 1 - eps)
}

factorcopula_rpareto <- function(n, fixed_factor) {
  fixed_factor <- coerce_factorcopula_spec(fixed_factor)
  u <- stats::runif(n)
  fixed_factor$pareto_scale / (1 - u)^(1 / fixed_factor$pareto_shape)
}

build_legendre_unit_rule <- function(n) {
  gl <- pracma::gaussLegendre(n, -1, 1)
  list(
    nodes = 0.5 * (gl$x + 1),
    weights = 0.5 * gl$w
  )
}

factorcopula_quadrature <- function(fixed_factor) {
  fixed_factor <- coerce_factorcopula_spec(fixed_factor)
  base_rule <- build_legendre_unit_rule(fixed_factor$quadrature_n)
  y <- 1 - base_rule$nodes
  list(
    nodes = base_rule$nodes,
    v = fixed_factor$pareto_scale / y,
    weights = base_rule$weights * fixed_factor$pareto_shape * y^(fixed_factor$pareto_shape - 1)
  )
}

bivariate_standard_normal_logpdf_factor <- function(z1, z2, rho) {
  validate_correlation(rho, name = "rho")
  denom <- 1 - rho^2
  -log(2 * pi) - 0.5 * log(denom) - (z1^2 - 2 * rho * z1 * z2 + z2^2) / (2 * denom)
}

factorcopula_row_logsumexp <- function(log_mat) {
  row_max <- apply(log_mat, 1, max)
  row_max + log(rowSums(exp(log_mat - row_max)))
}

factorcopula_marginal_cdf_exact <- function(w, fixed_factor) {
  fixed_factor <- coerce_factorcopula_spec(fixed_factor)
  w <- as.numeric(w)
  if (!length(w)) {
    return(numeric(0))
  }
  quad <- factorcopula_quadrature(fixed_factor)
  shifts <- outer(w, quad$v, FUN = "-")
  as.vector(stats::pnorm(shifts) %*% quad$weights)
}

factorcopula_marginal_pdf_exact <- function(w, fixed_factor) {
  fixed_factor <- coerce_factorcopula_spec(fixed_factor)
  w <- as.numeric(w)
  if (!length(w)) {
    return(numeric(0))
  }
  quad <- factorcopula_quadrature(fixed_factor)
  shifts <- outer(w, quad$v, FUN = "-")
  as.vector(stats::dnorm(shifts) %*% quad$weights)
}

factorcopula_find_w_limits <- function(fixed_factor) {
  fixed_factor <- coerce_factorcopula_spec(fixed_factor)
  lower <- -8
  upper <- 12
  lower_target <- fixed_factor$lower_tail_target
  upper_target <- fixed_factor$upper_tail_target

  repeat {
    lower_prob <- factorcopula_marginal_cdf_exact(lower, fixed_factor)
    if (lower_prob <= lower_target || lower <= -30) {
      break
    }
    lower <- lower - 2
  }

  repeat {
    upper_tail <- 1 - factorcopula_marginal_cdf_exact(upper, fixed_factor)
    if (upper_tail <= upper_target || upper >= 80) {
      break
    }
    upper <- upper + 4
  }

  c(lower, upper)
}

.factorcopula_cache <- new.env(parent = emptyenv())

factorcopula_spec_key <- function(fixed_factor) {
  fixed_factor <- coerce_factorcopula_spec(fixed_factor)
  paste(
    fixed_factor$pareto_shape,
    fixed_factor$pareto_scale,
    fixed_factor$quadrature_n,
    fixed_factor$lookup_grid_size,
    fixed_factor$lower_tail_target,
    fixed_factor$upper_tail_target,
    sep = "|"
  )
}

build_factorcopula_lookup <- function(fixed_factor) {
  fixed_factor <- coerce_factorcopula_spec(fixed_factor)
  limits <- factorcopula_find_w_limits(fixed_factor)
  w_grid <- seq(limits[1], limits[2], length.out = fixed_factor$lookup_grid_size)
  cdf_grid <- factorcopula_marginal_cdf_exact(w_grid, fixed_factor)
  pdf_grid <- factorcopula_marginal_pdf_exact(w_grid, fixed_factor)
  keep <- c(TRUE, diff(cdf_grid) > 1e-12)
  cdf_unique <- cdf_grid[keep]
  w_unique <- w_grid[keep]
  list(
    fixed_factor = fixed_factor,
    limits = limits,
    w_grid = w_grid,
    cdf_grid = cdf_grid,
    pdf_grid = pmax(pdf_grid, 1e-300),
    cdf_fun = stats::approxfun(w_grid, cdf_grid, rule = 2),
    pdf_fun = stats::approxfun(w_grid, pmax(pdf_grid, 1e-300), rule = 2),
    q_fun = stats::approxfun(cdf_unique, w_unique, rule = 2, ties = "ordered"),
    p_bounds = range(cdf_unique)
  )
}

get_factorcopula_lookup <- function(fixed_factor, rebuild = FALSE) {
  key <- factorcopula_spec_key(fixed_factor)
  if (rebuild || !exists(key, envir = .factorcopula_cache, inherits = FALSE)) {
    assign(key, build_factorcopula_lookup(fixed_factor), envir = .factorcopula_cache)
  }
  get(key, envir = .factorcopula_cache, inherits = FALSE)
}

factorcopula_fw <- function(w, fixed_factor) {
  lookup <- get_factorcopula_lookup(fixed_factor)
  clip_unit_response(lookup$cdf_fun(as.numeric(w)), eps = min(1e-8, 0.5 * min(diff(lookup$p_bounds))))
}

factorcopula_dw <- function(w, fixed_factor) {
  lookup <- get_factorcopula_lookup(fixed_factor)
  pmax(as.numeric(lookup$pdf_fun(as.numeric(w))), 1e-300)
}

factorcopula_qw <- function(u, fixed_factor) {
  validate_probabilities(u, "u")
  lookup <- get_factorcopula_lookup(fixed_factor)
  clipped_u <- pmin(pmax(as.numeric(u), lookup$p_bounds[1]), lookup$p_bounds[2])
  as.numeric(lookup$q_fun(clipped_u))
}

reference_quantile_factorcopula <- function(w1, w2, rho, fixed_factor, log = TRUE) {
  n <- max(length(w1), length(w2), length(rho))
  w1 <- recycle_to_length(w1, n, "w1")
  w2 <- recycle_to_length(w2, n, "w2")
  rho <- recycle_to_length(rho, n, "rho")
  validate_correlation(rho)
  fixed_factor <- coerce_factorcopula_spec(fixed_factor)
  quad <- factorcopula_quadrature(fixed_factor)

  log_terms <- matrix(0, nrow = n, ncol = length(quad$v))
  for (k in seq_along(quad$v)) {
    log_terms[, k] <- log(quad$weights[k]) +
      bivariate_standard_normal_logpdf_factor(
        w1 - quad$v[k],
        w2 - quad$v[k],
        rho = rho
      )
  }

  joint_log <- factorcopula_row_logsumexp(log_terms)
  out_log <- joint_log - log(factorcopula_dw(w1, fixed_factor)) - log(factorcopula_dw(w2, fixed_factor))
  if (log) out_log else exp(out_log)
}

reference_local_factorcopula <- function(u, v, rho, fixed_factor, log = TRUE) {
  n <- max(length(u), length(v), length(rho))
  u <- recycle_to_length(u, n, "u")
  v <- recycle_to_length(v, n, "v")
  rho <- recycle_to_length(rho, n, "rho")
  validate_probabilities(u, "u")
  validate_probabilities(v, "v")
  validate_correlation(rho)

  w1 <- factorcopula_qw(u, fixed_factor)
  w2 <- factorcopula_qw(v, fixed_factor)
  reference_quantile_factorcopula(w1, w2, rho = rho, fixed_factor = fixed_factor, log = log)
}

transform_factorcopula_inputs <- function(u, v, rho, fixed_factor) {
  n <- max(length(u), length(v), length(rho))
  u <- recycle_to_length(u, n, "u")
  v <- recycle_to_length(v, n, "v")
  rho <- recycle_to_length(rho, n, "rho")
  validate_probabilities(u, "u")
  validate_probabilities(v, "v")
  validate_correlation(rho)

  w1 <- factorcopula_qw(u, fixed_factor)
  w2 <- factorcopula_qw(v, fixed_factor)
  d <- w1 - w2
  data.frame(
    x = 0.5 * (w1 + w2),
    q = d^2,
    t = safe_atanh(rho),
    stringsAsFactors = FALSE
  )
}

inverse_transform_factorcopula_inputs <- function(x, q, t, fixed_factor, branch = c("positive", "negative")) {
  branch <- match.arg(branch)
  n <- max(length(x), length(q), length(t))
  x <- recycle_to_length(x, n, "x")
  q <- recycle_to_length(q, n, "q")
  t <- recycle_to_length(t, n, "t")
  if (any(!is.finite(q)) || any(q < 0)) {
    stop("`q` must be finite and nonnegative.", call. = FALSE)
  }

  d <- sqrt(q)
  if (branch == "negative") {
    d <- -d
  }
  w1 <- x + 0.5 * d
  w2 <- x - 0.5 * d
  data.frame(
    u = factorcopula_fw(w1, fixed_factor),
    v = factorcopula_fw(w2, fixed_factor),
    rho = tanh(t),
    stringsAsFactors = FALSE
  )
}

reference_transformed_factorcopula <- function(x, q, t, fixed_factor, log = TRUE) {
  n <- max(length(x), length(q), length(t))
  x <- recycle_to_length(x, n, "x")
  q <- recycle_to_length(q, n, "q")
  t <- recycle_to_length(t, n, "t")
  if (any(!is.finite(q)) || any(q < 0)) {
    stop("`q` must be finite and nonnegative.", call. = FALSE)
  }
  d <- sqrt(q)
  w1 <- x + 0.5 * d
  w2 <- x - 0.5 * d
  reference_quantile_factorcopula(w1, w2, rho = tanh(t), fixed_factor = fixed_factor, log = log)
}

build_factorcopula_local_table <- function(box, grid, fixed_factor) {
  stopifnot(all(c("x", "q", "t") %in% names(box)))
  stopifnot(all(c("x", "q", "t") %in% names(grid)))
  fixed_factor <- coerce_factorcopula_spec(fixed_factor)

  grid_df <- build_transformed_evaluation_grid(box, grid)
  grid_df$value <- reference_transformed_factorcopula(
    grid_df$x,
    grid_df$q,
    grid_df$t,
    fixed_factor = fixed_factor,
    log = TRUE
  )
  grid_df$pareto_shape <- fixed_factor$pareto_shape
  grid_df$pareto_scale <- fixed_factor$pareto_scale
  grid_df
}
