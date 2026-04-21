logit <- function(p) {
  qlogis(p)
}

inv_logit <- function(x) {
  plogis(x)
}

recycle_to_length <- function(x, n, name) {
  if (length(x) == 1L) {
    return(rep_len(x, n))
  }
  if (length(x) != n) {
    stop(sprintf("`%s` must have length 1 or %d.", name, n), call. = FALSE)
  }
  x
}

validate_probabilities <- function(x, name) {
  if (any(!is.finite(x)) || any(x <= 0 | x >= 1)) {
    stop(sprintf("`%s` must lie strictly inside (0, 1).", name), call. = FALSE)
  }
}

validate_correlation <- function(r, name = "r") {
  if (any(!is.finite(r)) || any(abs(r) >= 1)) {
    stop(sprintf("`%s` must lie strictly inside (-1, 1).", name), call. = FALSE)
  }
}

safe_atanh <- function(r, eps = 1e-10) {
  atanh(pmax(pmin(r, 1 - eps), -1 + eps))
}

wendland_r4 <- function(distance, b) {
  pmax(1 - distance / b, 0)^4
}

affine_to_unit_interval <- function(x, limits) {
  (x - limits[1]) / (limits[2] - limits[1])
}

affine_to_chebyshev <- function(x, limits) {
  2 * affine_to_unit_interval(x, limits) - 1
}

affine_from_chebyshev <- function(z, limits) {
  limits[1] + 0.5 * (z + 1) * (limits[2] - limits[1])
}

clip_to_limits <- function(x, limits) {
  pmin(pmax(x, limits[1]), limits[2])
}

tensor_row_design <- function(...) {
  bases <- list(...)
  n <- nrow(bases[[1]])
  out <- bases[[1]]
  if (length(bases) == 1L) {
    return(out)
  }
  for (k in 2:length(bases)) {
    current <- bases[[k]]
    next_out <- matrix(0, nrow = n, ncol = ncol(out) * ncol(current))
    for (i in seq_len(n)) {
      next_out[i, ] <- as.vector(kronecker(current[i, ], out[i, ]))
    }
    out <- next_out
  }
  out
}

build_second_difference_penalty <- function(dims, lambda) {
  lambda <- modifyList(list(x = 0, q = 0, t = 0), lambda)
  Ix <- diag(dims$x)
  Iq <- diag(dims$q)
  It <- diag(dims$t)
  Dx <- diff(Ix, differences = 2)
  Dq <- diff(Iq, differences = 2)
  Dt <- diff(It, differences = 2)
  Px <- if (nrow(Dx) > 0) crossprod(Dx) else matrix(0, dims$x, dims$x)
  Pq <- if (nrow(Dq) > 0) crossprod(Dq) else matrix(0, dims$q, dims$q)
  Pt <- if (nrow(Dt) > 0) crossprod(Dt) else matrix(0, dims$t, dims$t)
  lambda$x * kronecker(It, kronecker(Iq, Px)) +
    lambda$q * kronecker(It, kronecker(Pq, Ix)) +
    lambda$t * kronecker(Pt, kronecker(Iq, Ix))
}

finite_difference_gradient <- function(fun, x, step = 1e-5, ...) {
  x <- as.numeric(x)
  grad <- numeric(length(x))
  for (k in seq_along(x)) {
    shift <- rep(0, length(x))
    shift[k] <- step
    grad[k] <- (fun(x + shift, ...) - fun(x - shift, ...)) / (2 * step)
  }
  grad
}

finite_difference_hessian <- function(fun, x, step = 1e-4, ...) {
  x <- as.numeric(x)
  p <- length(x)
  hess <- matrix(0, nrow = p, ncol = p)
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      ei <- rep(0, p)
      ej <- rep(0, p)
      ei[i] <- step
      ej[j] <- step
      hess[i, j] <- (
        fun(x + ei + ej, ...) -
          fun(x + ei - ej, ...) -
          fun(x - ei + ej, ...) +
          fun(x - ei - ej, ...)
      ) / (4 * step^2)
    }
  }
  hess
}

transform_clayton_inputs <- function(u, v, r) {
  n <- max(length(u), length(v), length(r))
  u <- recycle_to_length(u, n, "u")
  v <- recycle_to_length(v, n, "v")
  r <- recycle_to_length(r, n, "r")
  validate_probabilities(u, "u")
  validate_probabilities(v, "v")
  validate_correlation(r)

  d <- logit(u) - logit(v)
  data.frame(
    x = 0.5 * (logit(u) + logit(v)),
    q = d^2,
    t = safe_atanh(r),
    stringsAsFactors = FALSE
  )
}

inverse_transform_clayton_inputs <- function(x, q, t, branch = c("positive", "negative")) {
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
  data.frame(
    u = inv_logit(x + 0.5 * d),
    v = inv_logit(x - 0.5 * d),
    r = tanh(t),
    stringsAsFactors = FALSE
  )
}

appell_f4_series_scalar <- function(a, b, c, d, x, y, tol = 1e-12, max_m = 400L, max_n = 400L) {
  if (!all(is.finite(c(a, b, c, d, x, y)))) {
    return(NaN)
  }
  if (abs(x) >= 1 || abs(y) >= 1) {
    return(NaN)
  }
  if (sqrt(abs(x)) + sqrt(abs(y)) >= 1) {
    return(NaN)
  }

  total <- 0
  row_start <- 1

  for (m in 0:max_m) {
    term <- row_start
    row_sum <- term

    if (m < max_m) {
      for (n in 0:(max_n - 1L)) {
        term <- term * ((a + m + n) * (b + m + n)) / ((d + n) * (n + 1)) * y
        row_sum <- row_sum + term
        if (abs(term) <= tol * max(1, abs(row_sum))) {
          break
        }
      }
    }

    total <- total + row_sum

    if (m == max_m) {
      break
    }

    row_start <- row_start * ((a + m) * (b + m)) / ((c + m) * (m + 1)) * x
    if (abs(row_start) <= tol * max(1, abs(total))) {
      break
    }
  }

  total
}

appell_f4_series <- function(a, b, c, d, x, y, tol = 1e-12, max_m = 400L, max_n = 400L) {
  n <- max(length(a), length(b), length(c), length(d), length(x), length(y))
  a <- recycle_to_length(a, n, "a")
  b <- recycle_to_length(b, n, "b")
  c <- recycle_to_length(c, n, "c")
  d <- recycle_to_length(d, n, "d")
  x <- recycle_to_length(x, n, "x")
  y <- recycle_to_length(y, n, "y")

  out <- numeric(n)
  for (i in seq_len(n)) {
    out[i] <- appell_f4_series_scalar(
      a[i], b[i], c[i], d[i], x[i], y[i],
      tol = tol, max_m = max_m, max_n = max_n
    )
  }
  out
}

reference_local_clayton <- function(u, v, r, nu, log = TRUE, tol = 1e-12, max_m = 400L, max_n = 400L) {
  n <- max(length(u), length(v), length(r), length(nu))
  u <- recycle_to_length(u, n, "u")
  v <- recycle_to_length(v, n, "v")
  r <- recycle_to_length(r, n, "r")
  nu <- recycle_to_length(nu, n, "nu")
  validate_probabilities(u, "u")
  validate_probabilities(v, "v")
  validate_correlation(r)
  if (any(!is.finite(nu)) || any(nu <= 0)) {
    stop("`nu` must be strictly positive.", call. = FALSE)
  }

  nu2 <- nu / 2
  a <- nu2 + 1
  rho2 <- r^2
  u_pow <- u^(1 / nu2)
  v_pow <- v^(1 / nu2)
  x <- rho2 * u_pow * v_pow
  y <- rho2 * (1 - u_pow) * (1 - v_pow)
  f4 <- appell_f4_series(a, a, nu2, rep(1, n), x, y, tol = tol, max_m = max_m, max_n = max_n)
  out_log <- a * log1p(-rho2) + log(f4)

  if (log) out_log else exp(out_log)
}

reference_transformed_clayton <- function(x, q, t, nu, log = TRUE, tol = 1e-12, max_m = 400L, max_n = 400L) {
  raw_inputs <- inverse_transform_clayton_inputs(x, q, t, branch = "positive")
  reference_local_clayton(
    raw_inputs$u,
    raw_inputs$v,
    raw_inputs$r,
    nu = nu,
    log = log,
    tol = tol,
    max_m = max_m,
    max_n = max_n
  )
}

geomodels_external_reference_clayton <- function(u, v, r, nu, log = TRUE) {
  if (!requireNamespace("GeoModels", quietly = TRUE)) {
    stop("Package `GeoModels` is required for the external reference evaluator.", call. = FALSE)
  }
  n <- max(length(u), length(v), length(r), length(nu))
  u <- recycle_to_length(u, n, "u")
  v <- recycle_to_length(v, n, "v")
  r <- recycle_to_length(r, n, "r")
  nu <- recycle_to_length(nu, n, "nu")

  out <- numeric(n)
  for (i in seq_len(n)) {
    out[i] <- .C(
      "biv_unif_CopulaClayton_call",
      as.double(u[i]),
      as.double(v[i]),
      as.double(r[i]),
      as.double(nu[i]),
      result = double(1),
      PACKAGE = "GeoModels"
    )$result
  }
  if (log) out else exp(out)
}

build_transformed_evaluation_grid <- function(box, grid) {
  dims <- names(box)
  if (length(dims) == 0L) {
    stop("`box` must contain at least one named coordinate range.", call. = FALSE)
  }
  if (!all(dims %in% names(grid))) {
    stop("`grid` must provide lengths for every coordinate in `box`.", call. = FALSE)
  }

  grid_axes <- lapply(dims, function(dim_name) {
    seq(box[[dim_name]][1], box[[dim_name]][2], length.out = grid[[dim_name]])
  })
  names(grid_axes) <- dims

  expand.grid(
    grid_axes,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}

build_clayton_local_table <- function(box, grid, nu, backend = c("series_reference", "series", "geomodels"), tol = 1e-12, max_m = 400L, max_n = 400L) {
  backend <- match.arg(backend)
  stopifnot(all(c("x", "q", "t") %in% names(box)))
  stopifnot(all(c("x", "q", "t") %in% names(grid)))

  grid_df <- build_transformed_evaluation_grid(box, grid)
  nu_eval <- if ("nu" %in% names(grid_df)) grid_df$nu else recycle_to_length(nu, nrow(grid_df), "nu")
  evaluator <- switch(
    backend,
    series_reference = function(x, q, t, nu_arg) {
      reference_transformed_clayton(x, q, t, nu = nu_arg, log = TRUE, tol = tol, max_m = max_m, max_n = max_n)
    },
    series = function(x, q, t, nu_arg) {
      reference_transformed_clayton(x, q, t, nu = nu_arg, log = TRUE, tol = tol, max_m = max_m, max_n = max_n)
    },
    geomodels = function(x, q, t, nu_arg) {
      uvrt <- inverse_transform_clayton_inputs(x, q, t, branch = "positive")
      geomodels_external_reference_clayton(uvrt$u, uvrt$v, uvrt$r, nu = nu_arg, log = TRUE)
    }
  )

  grid_df$value <- evaluator(grid_df$x, grid_df$q, grid_df$t, nu_eval)
  if (!"nu" %in% names(grid_df)) {
    grid_df$nu <- nu_eval
  }
  grid_df$backend <- backend
  grid_df
}
