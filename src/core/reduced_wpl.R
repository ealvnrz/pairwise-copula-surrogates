compute_beta_mean <- function(X, beta) {
  X <- as.matrix(X)
  beta <- as.numeric(beta)
  drop(plogis(X %*% beta))
}

build_nn_pair_design <- function(coords, m = 2L) {
  coords <- as.matrix(coords)
  if (ncol(coords) < 2L) {
    stop("`coords` must contain at least two spatial columns.", call. = FALSE)
  }
  n <- nrow(coords)
  if (n <= 1L) {
    stop("At least two spatial locations are required.", call. = FALSE)
  }
  dmat <- as.matrix(dist(coords[, 1:2, drop = FALSE]))
  diag(dmat) <- Inf
  pair_i <- integer(0)
  pair_j <- integer(0)
  pair_d <- numeric(0)
  for (i in seq_len(n)) {
    ord <- order(dmat[i, ], decreasing = FALSE)
    keep <- ord[seq_len(min(m, n - 1L))]
    pair_i <- c(pair_i, rep.int(i, length(keep)))
    pair_j <- c(pair_j, keep)
    pair_d <- c(pair_d, dmat[i, keep])
  }
  list(
    pair_index = cbind(i = pair_i, j = pair_j),
    distances = pair_d,
    weights = rep(1, length(pair_d)),
    m = m,
    n = n
  )
}

build_unique_nn_pair_design <- function(coords, m = 2L) {
  directed <- build_nn_pair_design(coords, m = m)
  pair_i <- pmin(directed$pair_index[, "i"], directed$pair_index[, "j"])
  pair_j <- pmax(directed$pair_index[, "i"], directed$pair_index[, "j"])
  keys <- paste(pair_i, pair_j, sep = "_")
  keep <- !duplicated(keys)

  list(
    pair_index = cbind(i = pair_i[keep], j = pair_j[keep]),
    distances = directed$distances[keep],
    weights = rep(1, sum(keep)),
    m = m,
    n = directed$n
  )
}

coerce_pair_weights <- function(weights, design) {
  if (is.null(weights)) {
    return(design$weights)
  }
  if (is.matrix(weights)) {
    idx <- design$pair_index
    return(weights[cbind(idx[, "i"], idx[, "j"])])
  }
  recycle_to_length(weights, nrow(design$pair_index), "weights")
}

prepare_pairwise_beta_data <- function(y, coords, X, m = 2L) {
  list(
    y = clip_unit_response(y),
    coords = as.matrix(coords),
    X = as.matrix(X),
    design = build_nn_pair_design(coords, m = m)
  )
}

prepare_pairwise_factor_data <- function(Y, coords, m = 6L, fixed_factor = NULL) {
  Y <- as.matrix(Y)
  if (!is.numeric(Y) || !all(is.finite(Y))) {
    stop("`Y` must be a finite numeric matrix of replicated spatial observations.", call. = FALSE)
  }
  if (nrow(Y) < 5L || ncol(Y) < 2L) {
    stop("`Y` must have at least 5 replicates and 2 spatial sites.", call. = FALSE)
  }
  fixed_factor <- coerce_factorcopula_spec(fixed_factor)
  design <- build_unique_nn_pair_design(coords, m = m)
  idx <- design$pair_index

  u_rank <- apply(Y, 2, function(column) rank(column, ties.method = "average") / (length(column) + 1))
  u_rank <- as.matrix(u_rank)
  w_rank <- matrix(
    factorcopula_qw(as.vector(u_rank), fixed_factor = fixed_factor),
    nrow = nrow(Y),
    ncol = ncol(Y)
  )
  left <- w_rank[, idx[, "i"], drop = FALSE]
  right <- w_rank[, idx[, "j"], drop = FALSE]

  list(
    Y = Y,
    coords = as.matrix(coords),
    fixed_factor = fixed_factor,
    u = u_rank,
    w = w_rank,
    w_left_matrix = left,
    w_right_matrix = right,
    x_pair_matrix = 0.5 * (left + right),
    q_pair_matrix = (left - right)^2,
    design = design,
    n_rep = nrow(Y),
    n_sites = ncol(Y),
    n_pairs = nrow(idx)
  )
}

prepare_reduced_clayton_data <- function(y, coords, X, m = 2L) {
  prepare_pairwise_beta_data(y = y, coords = coords, X = X, m = m)
}

practical_theta_parts <- function(theta, data) {
  theta <- as.numeric(theta)
  p <- ncol(as.matrix(data$X))
  if (length(theta) != p + 2L) {
    stop(sprintf("`theta` must have length %d = p + 2.", p + 2L), call. = FALSE)
  }
  list(
    beta = theta[seq_len(p)],
    delta = theta[p + 1L],
    b = theta[p + 2L]
  )
}

practical_parameter_names <- function(data) {
  p <- ncol(as.matrix(data$X))
  beta_names <- if (!is.null(colnames(data$X))) {
    paste0("beta_", colnames(data$X))
  } else {
    paste0("beta", seq_len(p))
  }
  c(beta_names, "delta", "b")
}

build_clayton_components_common <- function(delta, b, nu, data, beta, weights = NULL, design = NULL) {
  if (!is.finite(delta) || delta <= 0 || !is.finite(b) || b <= 0 || !is.finite(nu) || nu <= 0) {
    return(NULL)
  }

  design <- if (is.null(design)) data$design else design
  weights <- coerce_pair_weights(weights, design)
  idx <- design$pair_index
  mu <- compute_beta_mean(data$X, beta)
  alpha <- mu * delta
  gamma <- (1 - mu) * delta
  u <- clip_unit_response(pbeta(data$y, shape1 = alpha, shape2 = gamma))
  r <- wendland_r4(design$distances, b)

  list(
    theta = c(delta = delta, b = b, nu = nu),
    beta = as.numeric(beta),
    weights = weights,
    design = design,
    index = idx,
    mu = mu,
    alpha = alpha,
    gamma = gamma,
    u = u,
    r = r,
    nu = nu
  )
}

build_reduced_clayton_components <- function(theta, data, fixed_beta, fixed_nu, weights = NULL, design = NULL) {
  theta <- as.numeric(theta)
  if (length(theta) != 2L) {
    stop("`theta` must be `(delta, b)`.", call. = FALSE)
  }
  build_clayton_components_common(
    delta = theta[1],
    b = theta[2],
    nu = fixed_nu,
    data = data,
    beta = fixed_beta,
    weights = weights,
    design = design
  )
}

build_practical_clayton_components <- function(theta, data, fixed_nu, weights = NULL, design = NULL) {
  parts <- practical_theta_parts(theta, data)
  build_clayton_components_common(
    delta = parts$delta,
    b = parts$b,
    nu = fixed_nu,
    data = data,
    beta = parts$beta,
    weights = weights,
    design = design
  )
}

build_factorcopula_components <- function(theta, data, fixed_factor, weights = NULL, design = NULL) {
  theta <- as.numeric(theta)
  if (length(theta) != 2L) {
    stop("`theta` must be `(alpha, b)`.", call. = FALSE)
  }
  alpha <- theta[1]
  b <- theta[2]
  if (!is.finite(alpha) || alpha <= 0 || alpha > 2 || !is.finite(b) || b <= 0) {
    return(NULL)
  }

  design <- if (is.null(design)) data$design else design
  weights <- coerce_pair_weights(weights, design)
  rho_pair <- spatial_factor_correlation(design$distances, alpha = alpha, b = b)

  list(
    theta = c(alpha = alpha, b = b),
    weights = weights,
    weights_vec = rep(weights, each = data$n_rep),
    design = design,
    index = design$pair_index,
    rho_pair = rho_pair,
    t_pair = safe_atanh(rho_pair),
    t_vec = rep(safe_atanh(rho_pair), each = data$n_rep),
    w1_vec = as.vector(data$w_left_matrix),
    w2_vec = as.vector(data$w_right_matrix),
    x_vec = as.vector(data$x_pair_matrix),
    q_vec = as.vector(data$q_pair_matrix),
    fixed_factor = coerce_factorcopula_spec(fixed_factor)
  )
}

build_stress_clayton_components <- function(theta, data, fixed_beta, weights = NULL, design = NULL) {
  theta <- as.numeric(theta)
  if (length(theta) != 3L) {
    stop("`theta` must be `(delta, b, nu)`.", call. = FALSE)
  }
  build_clayton_components_common(
    delta = theta[1],
    b = theta[2],
    nu = theta[3],
    data = data,
    beta = fixed_beta,
    weights = weights,
    design = design
  )
}

pairwise_margin_terms <- function(data, comp) {
  idx <- comp$index
  dbeta(data$y[idx[, "i"]], shape1 = comp$alpha[idx[, "i"]], shape2 = comp$gamma[idx[, "i"]], log = TRUE) +
    dbeta(data$y[idx[, "j"]], shape1 = comp$alpha[idx[, "j"]], shape2 = comp$gamma[idx[, "j"]], log = TRUE)
}

wpl_reference_reduced <- function(theta, data, fixed_beta, fixed_nu, weights = NULL, design = NULL, tol = 1e-12, max_m = 400L, max_n = 400L) {
  comp <- build_reduced_clayton_components(theta, data, fixed_beta, fixed_nu, weights = weights, design = design)
  if (is.null(comp)) {
    return(-Inf)
  }

  idx <- comp$index
  local <- reference_local_clayton(
    comp$u[idx[, "i"]],
    comp$u[idx[, "j"]],
    comp$r,
    nu = comp$nu,
    log = TRUE,
    tol = tol,
    max_m = max_m,
    max_n = max_n
  )
  margins <- pairwise_margin_terms(data, comp)
  sum(comp$weights * (local + margins))
}

wpl_surrogate_reduced <- function(theta, data, surrogate_fit, fixed_beta, fixed_nu, weights = NULL, design = NULL) {
  comp <- build_reduced_clayton_components(theta, data, fixed_beta, fixed_nu, weights = weights, design = design)
  if (is.null(comp)) {
    return(-Inf)
  }

  idx <- comp$index
  local <- predict_local_surrogate(surrogate_fit, comp$u[idx[, "i"]], comp$u[idx[, "j"]], comp$r, nu = comp$nu)
  margins <- pairwise_margin_terms(data, comp)
  sum(comp$weights * (local + margins))
}

wpl_reference_practical <- function(theta, data, fixed_nu, weights = NULL, design = NULL, tol = 1e-12, max_m = 400L, max_n = 400L) {
  comp <- build_practical_clayton_components(theta, data, fixed_nu, weights = weights, design = design)
  if (is.null(comp)) {
    return(-Inf)
  }

  idx <- comp$index
  local <- reference_local_clayton(
    comp$u[idx[, "i"]],
    comp$u[idx[, "j"]],
    comp$r,
    nu = comp$nu,
    log = TRUE,
    tol = tol,
    max_m = max_m,
    max_n = max_n
  )
  margins <- pairwise_margin_terms(data, comp)
  sum(comp$weights * (local + margins))
}

wpl_surrogate_practical <- function(theta, data, surrogate_fit, fixed_nu, weights = NULL, design = NULL) {
  comp <- build_practical_clayton_components(theta, data, fixed_nu, weights = weights, design = design)
  if (is.null(comp)) {
    return(-Inf)
  }

  idx <- comp$index
  local <- predict_local_surrogate(surrogate_fit, comp$u[idx[, "i"]], comp$u[idx[, "j"]], comp$r, nu = comp$nu)
  margins <- pairwise_margin_terms(data, comp)
  sum(comp$weights * (local + margins))
}

wpl_reference_factorcopula <- function(theta, data, fixed_factor, weights = NULL, design = NULL) {
  comp <- build_factorcopula_components(theta, data, fixed_factor = fixed_factor, weights = weights, design = design)
  if (is.null(comp)) {
    return(-1e12)
  }

  local <- reference_quantile_factorcopula(
    comp$w1_vec,
    comp$w2_vec,
    rho = rep(comp$rho_pair, each = data$n_rep),
    fixed_factor = comp$fixed_factor,
    log = TRUE
  )
  if (any(!is.finite(local))) {
    return(-1e12)
  }
  out <- sum(comp$weights_vec * local)
  if (!is.finite(out)) -1e12 else out
}

prepare_factorcopula_surrogate_cache <- function(data, surrogate_fit) {
  x_values <- clip_to_limits(as.vector(data$x_pair_matrix), surrogate_fit$box$x)
  q_values <- clip_to_limits(as.vector(data$q_pair_matrix), surrogate_fit$box$q)
  if (inherits(surrogate_fit, "clayton_spline_surrogate")) {
    basis_x <- make_bs_basis(
      x_values,
      df = surrogate_fit$dims$x,
      degree = surrogate_fit$basis$x$degree,
      boundary_knots = surrogate_fit$basis$x$boundary,
      knots = surrogate_fit$basis$x$knots
    )
    basis_q <- make_bs_basis(
      q_values,
      df = surrogate_fit$dims$q,
      degree = surrogate_fit$basis$q$degree,
      boundary_knots = surrogate_fit$basis$q$boundary,
      knots = surrogate_fit$basis$q$knots
    )
    coef_array <- array(
      surrogate_fit$coefficients,
      dim = c(surrogate_fit$dims$x, surrogate_fit$dims$q, surrogate_fit$dims$t)
    )
    static_array <- array(0, dim = c(data$n_rep, data$n_pairs, surrogate_fit$dims$t))
    for (k in seq_len(surrogate_fit$dims$t)) {
      slice <- coef_array[, , k, drop = TRUE]
      static_array[, , k] <- matrix(
        rowSums((basis_x %*% slice) * basis_q),
        nrow = data$n_rep,
        ncol = data$n_pairs
      )
    }
    return(structure(
      list(
        fit = surrogate_fit,
        static_array = static_array,
        n_rep = data$n_rep,
        n_pairs = data$n_pairs
      ),
      class = "factorcopula_surrogate_cache"
    ))
  }

  if (inherits(surrogate_fit, "clayton_chebyshev_surrogate")) {
    basis_x <- chebyshev_basis(
      affine_to_chebyshev(x_values, surrogate_fit$box$x),
      surrogate_fit$degree$x
    )
    basis_q <- chebyshev_basis(
      affine_to_chebyshev(q_values, surrogate_fit$box$q),
      surrogate_fit$degree$q
    )
    coef_array <- array(
      surrogate_fit$coefficients,
      dim = c(
        surrogate_fit$degree$x + 1L,
        surrogate_fit$degree$q + 1L,
        surrogate_fit$degree$t + 1L
      )
    )
    static_array <- array(0, dim = c(data$n_rep, data$n_pairs, surrogate_fit$degree$t + 1L))
    for (k in seq_len(surrogate_fit$degree$t + 1L)) {
      slice <- coef_array[, , k, drop = TRUE]
      static_array[, , k] <- matrix(
        rowSums((basis_x %*% slice) * basis_q),
        nrow = data$n_rep,
        ncol = data$n_pairs
      )
    }
    return(structure(
      list(
        fit = surrogate_fit,
        static_array = static_array,
        n_rep = data$n_rep,
        n_pairs = data$n_pairs
      ),
      class = "factorcopula_surrogate_cache"
    ))
  }

  stop("Unsupported surrogate class for factor-copula caching.", call. = FALSE)
}

evaluate_factorcopula_surrogate_cached <- function(cache, t_pair) {
  if (!inherits(cache, "factorcopula_surrogate_cache")) {
    stop("`cache` must be created by `prepare_factorcopula_surrogate_cache()`.", call. = FALSE)
  }
  t_pair <- recycle_to_length(t_pair, cache$n_pairs, "t_pair")
  t_pair <- clip_to_limits(t_pair, cache$fit$box$t)

  if (inherits(cache$fit, "clayton_spline_surrogate")) {
    basis_t <- make_bs_basis(
      t_pair,
      df = cache$fit$dims$t,
      degree = cache$fit$basis$t$degree,
      boundary_knots = cache$fit$basis$t$boundary,
      knots = cache$fit$basis$t$knots
    )
  } else if (inherits(cache$fit, "clayton_chebyshev_surrogate")) {
    basis_t <- chebyshev_basis(
      affine_to_chebyshev(t_pair, cache$fit$box$t),
      cache$fit$degree$t
    )
  } else {
    stop("Unsupported surrogate class in factor-copula cache.", call. = FALSE)
  }

  pred_mat <- matrix(0, nrow = cache$n_rep, ncol = cache$n_pairs)
  for (k in seq_len(ncol(basis_t))) {
    pred_mat <- pred_mat + cache$static_array[, , k] * matrix(
      basis_t[, k],
      nrow = cache$n_rep,
      ncol = cache$n_pairs,
      byrow = TRUE
    )
  }
  as.vector(pred_mat)
}

wpl_surrogate_factorcopula <- function(theta, data, surrogate_fit, fixed_factor, weights = NULL, design = NULL) {
  comp <- build_factorcopula_components(theta, data, fixed_factor = fixed_factor, weights = weights, design = design)
  if (is.null(comp)) {
    return(-1e12)
  }

  surrogate_cache <- if (inherits(surrogate_fit, "factorcopula_surrogate_cache")) surrogate_fit else prepare_factorcopula_surrogate_cache(data, surrogate_fit)
  local <- evaluate_factorcopula_surrogate_cached(surrogate_cache, t_pair = comp$t_pair)
  if (any(!is.finite(local))) {
    return(-1e12)
  }
  out <- sum(comp$weights_vec * local)
  if (!is.finite(out)) -1e12 else out
}

wpl_reference_stress <- function(theta, data, fixed_beta, weights = NULL, design = NULL, tol = 1e-12, max_m = 400L, max_n = 400L) {
  comp <- build_stress_clayton_components(theta, data, fixed_beta, weights = weights, design = design)
  if (is.null(comp)) {
    return(-Inf)
  }

  idx <- comp$index
  local <- reference_local_clayton(
    comp$u[idx[, "i"]],
    comp$u[idx[, "j"]],
    comp$r,
    nu = comp$nu,
    log = TRUE,
    tol = tol,
    max_m = max_m,
    max_n = max_n
  )
  margins <- pairwise_margin_terms(data, comp)
  sum(comp$weights * (local + margins))
}

wpl_surrogate_stress <- function(theta, data, surrogate_fit, fixed_beta, weights = NULL, design = NULL) {
  comp <- build_stress_clayton_components(theta, data, fixed_beta, weights = weights, design = design)
  if (is.null(comp)) {
    return(-Inf)
  }

  idx <- comp$index
  local <- predict_local_surrogate(surrogate_fit, comp$u[idx[, "i"]], comp$u[idx[, "j"]], comp$r, nu = comp$nu)
  margins <- pairwise_margin_terms(data, comp)
  sum(comp$weights * (local + margins))
}

optimize_reduced_wpl <- function(objective, start, lower, upper, ...) {
  optim(
    par = start,
    fn = function(par, ...) -objective(par, ...),
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    ...
  )
}
