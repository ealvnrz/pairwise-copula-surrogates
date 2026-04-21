make_bs_basis <- function(x, df, degree, boundary_knots, knots = NULL) {
  splines::bs(
    x,
    df = if (is.null(knots)) df else NULL,
    degree = degree,
    knots = knots,
    Boundary.knots = boundary_knots,
    intercept = TRUE
  )
}

build_coordinate_penalty <- function(dims, lambda) {
  coord_names <- names(dims)
  lambda_defaults <- setNames(as.list(rep(0, length(coord_names))), coord_names)
  lambda <- modifyList(lambda_defaults, lambda)

  penalties <- lapply(coord_names, function(coord_name) {
    n_basis <- dims[[coord_name]]
    identity_mat <- diag(n_basis)
    diff_mat <- diff(identity_mat, differences = 2)
    penalty_core <- if (nrow(diff_mat) > 0L) crossprod(diff_mat) else matrix(0, n_basis, n_basis)
    factors <- lapply(coord_names, function(name) {
      if (name == coord_name) {
        penalty_core
      } else {
        diag(dims[[name]])
      }
    })
    kron <- factors[[1]]
    if (length(factors) > 1L) {
      for (k in 2:length(factors)) {
        kron <- kronecker(factors[[k]], kron)
      }
    }
    lambda[[coord_name]] * kron
  })

  Reduce(`+`, penalties)
}

chebyshev_basis <- function(z, degree) {
  z <- as.numeric(z)
  out <- matrix(0, nrow = length(z), ncol = degree + 1L)
  out[, 1] <- 1
  if (degree >= 1L) {
    out[, 2] <- z
  }
  if (degree >= 2L) {
    for (k in 2:degree) {
      out[, k + 1L] <- 2 * z * out[, k] - out[, k - 1L]
    }
  }
  out
}

chebyshev_nodes <- function(degree, limits) {
  k <- 0:degree
  z <- cos((2 * k + 1) * pi / (2 * (degree + 1)))
  affine_from_chebyshev(z, limits)
}

collect_prediction_coordinates <- function(fit, x, q = NULL, t = NULL, nu = NULL) {
  coord_values <- list(x = x, q = q, t = t, nu = nu)
  missing_coords <- fit$coord_names[vapply(fit$coord_names, function(name) is.null(coord_values[[name]]), logical(1))]
  if (length(missing_coords) > 0L) {
    stop(sprintf("Missing coordinates for prediction: %s.", paste(missing_coords, collapse = ", ")), call. = FALSE)
  }
  target_length <- max(vapply(fit$coord_names, function(name) length(coord_values[[name]]), integer(1)))
  out <- lapply(fit$coord_names, function(name) {
    clip_to_limits(
      recycle_to_length(coord_values[[name]], target_length, name),
      fit$box[[name]]
    )
  })
  names(out) <- fit$coord_names
  out
}

fit_spline_surrogate <- function(local_table, box, knots, penalty, nu) {
  coord_names <- names(box)
  stopifnot(all(c(coord_names, "value") %in% names(local_table)))
  if (is.null(knots)) {
    knots <- list()
  }
  degree <- if (!is.null(knots$degree)) knots$degree else 3L
  knot_defaults <- setNames(as.list(rep(10L, length(coord_names))), coord_names)
  knot_spec <- modifyList(knot_defaults, knots[coord_names])
  lambda <- if (length(penalty) == 1L) {
    setNames(as.list(rep(penalty, length(coord_names))), coord_names)
  } else {
    modifyList(setNames(as.list(rep(0, length(coord_names))), coord_names), penalty)
  }

  basis_mats <- lapply(coord_names, function(coord_name) {
    make_bs_basis(
      local_table[[coord_name]],
      df = knot_spec[[coord_name]],
      degree = degree,
      boundary_knots = box[[coord_name]]
    )
  })
  names(basis_mats) <- coord_names
  X <- do.call(tensor_row_design, basis_mats)
  dims <- lapply(basis_mats, ncol)
  penalty_matrix <- build_coordinate_penalty(dims, lambda = lambda)
  coef <- solve(crossprod(X) + penalty_matrix, crossprod(X, local_table$value))

  basis_info <- lapply(coord_names, function(coord_name) {
    B <- basis_mats[[coord_name]]
    list(
      knots = attr(B, "knots"),
      boundary = attr(B, "Boundary.knots"),
      degree = degree
    )
  })
  names(basis_info) <- coord_names

  structure(
    list(
      coefficients = coef,
      nu = if ("nu" %in% coord_names) NULL else nu,
      box = box,
      coord_names = coord_names,
      dims = dims,
      penalty = lambda,
      basis = basis_info
    ),
    class = "clayton_spline_surrogate"
  )
}

predict_spline_transformed <- function(fit, x, q, t, nu = NULL) {
  stopifnot(inherits(fit, "clayton_spline_surrogate"))
  coords <- collect_prediction_coordinates(fit, x = x, q = q, t = t, nu = nu)
  basis_mats <- lapply(fit$coord_names, function(coord_name) {
    basis_info <- fit$basis[[coord_name]]
    make_bs_basis(
      coords[[coord_name]],
      df = fit$dims[[coord_name]],
      degree = basis_info$degree,
      boundary_knots = basis_info$boundary,
      knots = basis_info$knots
    )
  })
  X <- do.call(tensor_row_design, basis_mats)
  drop(X %*% fit$coefficients)
}

predict_spline_surrogate <- function(fit, u, v, r, nu = NULL) {
  transformed <- transform_clayton_inputs(u, v, r)
  if ("nu" %in% fit$coord_names) {
    predict_spline_transformed(fit, transformed$x, transformed$q, transformed$t, nu = nu)
  } else {
    predict_spline_transformed(fit, transformed$x, transformed$q, transformed$t)
  }
}

fit_chebyshev_surrogate <- function(local_table, box, degree, nu) {
  coord_names <- names(box)
  stopifnot(all(c(coord_names, "value") %in% names(local_table)))
  degree_defaults <- setNames(as.list(rep(10L, length(coord_names))), coord_names)
  degree <- modifyList(degree_defaults, degree)

  basis_mats <- lapply(coord_names, function(coord_name) {
    chebyshev_basis(
      affine_to_chebyshev(local_table[[coord_name]], box[[coord_name]]),
      degree[[coord_name]]
    )
  })
  names(basis_mats) <- coord_names
  X <- do.call(tensor_row_design, basis_mats)
  coef <- qr.solve(X, local_table$value)

  structure(
    list(
      coefficients = coef,
      nu = if ("nu" %in% coord_names) NULL else nu,
      box = box,
      coord_names = coord_names,
      degree = degree
    ),
    class = "clayton_chebyshev_surrogate"
  )
}

predict_chebyshev_transformed <- function(fit, x, q, t, nu = NULL) {
  stopifnot(inherits(fit, "clayton_chebyshev_surrogate"))
  coords <- collect_prediction_coordinates(fit, x = x, q = q, t = t, nu = nu)
  basis_mats <- lapply(fit$coord_names, function(coord_name) {
    chebyshev_basis(
      affine_to_chebyshev(coords[[coord_name]], fit$box[[coord_name]]),
      fit$degree[[coord_name]]
    )
  })
  X <- do.call(tensor_row_design, basis_mats)
  drop(X %*% fit$coefficients)
}

predict_chebyshev_surrogate <- function(fit, u, v, r, nu = NULL) {
  transformed <- transform_clayton_inputs(u, v, r)
  if ("nu" %in% fit$coord_names) {
    predict_chebyshev_transformed(fit, transformed$x, transformed$q, transformed$t, nu = nu)
  } else {
    predict_chebyshev_transformed(fit, transformed$x, transformed$q, transformed$t)
  }
}

predict_local_surrogate <- function(surrogate_fit, u, v, r, nu = NULL) {
  if (inherits(surrogate_fit, "clayton_spline_surrogate")) {
    return(predict_spline_surrogate(surrogate_fit, u = u, v = v, r = r, nu = nu))
  }
  if (inherits(surrogate_fit, "clayton_chebyshev_surrogate")) {
    return(predict_chebyshev_surrogate(surrogate_fit, u = u, v = v, r = r, nu = nu))
  }
  stop("Unsupported surrogate class.", call. = FALSE)
}

evaluate_transformed_surrogate <- function(fit, x, q, t, nu = NULL) {
  if (inherits(fit, "clayton_spline_surrogate")) {
    return(predict_spline_transformed(fit, x, q, t, nu = nu))
  }
  if (inherits(fit, "clayton_chebyshev_surrogate")) {
    return(predict_chebyshev_transformed(fit, x, q, t, nu = nu))
  }
  stop("Unsupported surrogate class.", call. = FALSE)
}
