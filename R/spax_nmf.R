#' Regularized NMF for IFU Data
#'
#' Fits a non-negative matrix factorization model with an optional smoothness
#' penalty on the recovered spectra.
#'
#' The factorization follows the linear mixing model
#' \deqn{X \approx A S}
#' where rows of \eqn{X} correspond to spectra from individual spatial elements,
#' \eqn{A} contains non-negative component weights, and \eqn{S} contains
#' non-negative spectral components.
#'
#' @param x Numeric matrix-like object with spectra in rows and wavelength
#'   channels in columns.
#' @param k Number of components to estimate.
#' @param lambda_smooth Non-negative weight for the first-difference smoothness
#'   penalty applied to the recovered spectra.
#' @param lambda_spatial Non-negative weight for the grid-based smoothness
#'   penalty applied to the recovered component-weight maps.
#' @param lambda_sparse Non-negative weight for an L1 sparsity penalty applied
#'   to selected component-weight columns.
#' @param spatial_sigma Positive kernel width in spaxels for the Gaussian
#'   spatial weighting. Nearby spaxels receive stronger coupling than distant
#'   ones.
#' @param smooth_components Optional component indices that should receive the
#'   spatial smoothness prior when \code{lambda_spatial} is scalar.
#' @param sparse_components Optional component indices that should receive the
#'   sparsity prior when \code{lambda_sparse} is scalar.
#' @param lr Learning rate passed to \code{torch::optim_adam()}.
#' @param niter Number of optimization iterations.
#' @param nx Optional number of x-axis pixels. Required for spatial smoothing
#'   when it cannot be inferred from \code{x}.
#' @param ny Optional number of y-axis pixels. Required for spatial smoothing
#'   when it cannot be inferred from \code{x}.
#' @param center Logical; if \code{TRUE}, center wavelength channels before
#'   fitting.
#' @param scale Logical; if \code{TRUE}, scale wavelength channels before
#'   fitting.
#' @param cuda Logical; if \code{TRUE}, request CUDA execution.
#' @param verbose Logical; if \code{TRUE}, emit progress every 100 iterations.
#' @param metadata Optional metadata list to carry through the fitted object,
#'   reconstructed cubes, and prediction outputs.
#'
#' @return A list with class \code{"spax_nmf"} containing:
#' \itemize{
#'   \item \code{spatial}: component-weight matrix with one column per component.
#'   \item \code{weights}: alias of \code{spatial}.
#'   \item \code{spectra}: component spectra with one row per component.
#'   \item \code{basis}: component spectra arranged as wavelength-by-component.
#'   \item \code{coef}: component weights arranged as component-by-spaxel.
#'   \item \code{reconstruction}: reconstructed matrix \eqn{AS}.
#'   \item \code{fitted}: alias of \code{reconstruction}.
#'   \item \code{loss}: objective history over the optimization.
#'   \item \code{metadata}: optional metadata carried from the input matrix or
#'     provided explicitly.
#'   \item \code{center}: centering values used during preprocessing, or
#'     \code{FALSE}.
#'   \item \code{scale}: scaling values used during preprocessing, or
#'     \code{FALSE}.
#'   \item \code{call}: matched function call.
#' }
#'
#' @details
#' This implementation is aimed at integral-field spectroscopy data after
#' reshaping the spatial dimensions into a two-dimensional matrix of spaxels by
#' wavelength. The objective function is
#' \deqn{\|X - AS\|^2 + \lambda_S \|D S\|^2 + \sum_c \lambda_{A,c} \sum_{i,j} w_{ij}(A_{ic} - A_{jc})^2 + \sum_c \lambda_{C,c}\|A_{:,c}\|_1}
#' where the second term penalizes rapid channel-to-channel variations in each
#' component spectrum, the third penalizes sharp changes between nearby spatial
#' pixels using Gaussian kernel weights \eqn{w_{ij}} when \code{nx} and
#' \code{ny} are available, and the fourth promotes compact sparse spatial
#' components when desired.
#'
#' For the current package workflow, a typical sequence is:
#' \enumerate{
#'   \item reshape a cube with [cube_to_matrix()];
#'   \item fit the model with \code{spax_nmf()};
#'   \item recover component maps with [component_map()].
#' }
#'
#' @references
#' Lee, D. D., and Seung, H. S. (1999). Learning the parts of objects by
#' non-negative matrix factorization. \emph{Nature}, 401, 788-791.
#'
#' @examples
#' cube <- simulate_ifu_cube(nx = 12, ny = 12, n_wave = 80)
#' fit <- spax_nmf(cube$matrix, k = 3, niter = 150, lr = 0.03)
#' print(fit)
#'
#' map1 <- component_map(fit, nx = cube$nx, ny = cube$ny, component = 1)
#' dim(map1)
#'
#' @export
spax_nmf <- function(x,
                           k = 3,
                           lambda_smooth = 0.01,
                           lr = 0.02,
                           niter = 2000,
                           center = FALSE,
                           scale = FALSE,
                           cuda = FALSE,
                           verbose = FALSE,
                           metadata = NULL,
                           lambda_spatial = 0,
                           lambda_sparse = 0,
                           spatial_sigma = 1,
                           smooth_components = NULL,
                           sparse_components = NULL,
                           nx = NULL,
                           ny = NULL) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("The 'torch' package must be installed to use spax_nmf().", call. = FALSE)
  }

  input_metadata <- merge_metadata(extract_matrix_metadata(x), metadata)
  x <- validate_input_matrix(x)

  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1L) {
    stop("'k' must be a positive integer.", call. = FALSE)
  }
  k <- as.integer(k)

  if (k > min(dim(x))) {
    stop("'k' must not exceed the smaller input dimension.", call. = FALSE)
  }
  if (!is.numeric(lambda_smooth) || length(lambda_smooth) != 1L || is.na(lambda_smooth) ||
      lambda_smooth < 0) {
    stop("'lambda_smooth' must be a non-negative number.", call. = FALSE)
  }
  if (!is.numeric(spatial_sigma) || length(spatial_sigma) != 1L || is.na(spatial_sigma) ||
      spatial_sigma <= 0) {
    stop("'spatial_sigma' must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(lr) || length(lr) != 1L || is.na(lr) || lr <= 0) {
    stop("'lr' must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(niter) || length(niter) != 1L || is.na(niter) || niter < 1) {
    stop("'niter' must be a positive integer.", call. = FALSE)
  }
  niter <- as.integer(niter)

  lambda_spatial_weights <- resolve_component_penalty(
    lambda_spatial,
    k = k,
    active_components = smooth_components,
    name = "lambda_spatial"
  )
  lambda_sparse_weights <- resolve_component_penalty(
    lambda_sparse,
    k = k,
    active_components = sparse_components,
    name = "lambda_sparse"
  )

  spatial_dims <- resolve_spatial_dims(x, nx = nx, ny = ny, lambda_spatial = lambda_spatial_weights)
  x <- base::scale(x, center = center, scale = scale)
  cen <- attr(x, "scaled:center")
  sc <- attr(x, "scaled:scale")

  if (any(!is.finite(x))) {
    stop("Preprocessing produced non-finite values.", call. = FALSE)
  }
  if (any(x < 0)) {
    warning(
      "Input contains negative values after preprocessing. Interpretability of ",
      "non-negative components may be reduced.",
      call. = FALSE
    )
  }

  device <- select_torch_device(cuda)
  X <- torch::torch_tensor(as.matrix(x), device = device)

  n_spectra <- nrow(x)
  n_wave <- ncol(x)
  spatial_graph <- build_spatial_graph(
    spatial_dims,
    n_expected = n_spectra,
    spatial_sigma = spatial_sigma
  )

  A <- torch::torch_rand(n_spectra, k, device = device, requires_grad = TRUE)
  S <- torch::torch_rand(k, n_wave, device = device, requires_grad = TRUE)
  opt <- torch::optim_adam(list(A, S), lr = lr)

  loss_history <- numeric(niter)

  for (i in seq_len(niter)) {
    opt$zero_grad()

    Xhat <- A$matmul(S)
    recon <- torch::torch_mean((X - Xhat)^2)

    if (n_wave > 1L) {
      diff <- S[, 2:n_wave] - S[, 1:(n_wave - 1L)]
      smooth <- torch::torch_mean(diff^2)
    } else {
      smooth <- torch::torch_tensor(0, device = device)
    }

    spatial_smooth <- spatial_penalty(
      A,
      spatial_graph,
      device = device,
      component_weights = lambda_spatial_weights
    )
    sparse_abundance <- sparse_penalty(
      A,
      device = device,
      component_weights = lambda_sparse_weights
    )
    loss <- recon + lambda_smooth * smooth + spatial_smooth + sparse_abundance
    loss$backward()
    opt$step()

    A$data()$clamp_(min = 0)
    S$data()$clamp_(min = 0)

    loss_history[i] <- as.numeric(loss$item())
    if (verbose && (i %% 100L == 0L || i == 1L || i == niter)) {
      message(sprintf("Iteration %d/%d, loss = %.6f", i, niter, loss_history[i]))
    }
  }

  spatial <- as.matrix(A$detach()$cpu())
  spectra <- as.matrix(S$detach()$cpu())
  reconstruction <- spatial %*% spectra

  fit <- list(
    spatial = spatial,
    weights = spatial,
    abundance = spatial,
    spectra = spectra,
    basis = t(spectra),
    coef = t(spatial),
    reconstruction = reconstruction,
    fitted = reconstruction,
    loss = loss_history,
    metadata = input_metadata,
    spatial_dims = spatial_dims,
    lambda_smooth = lambda_smooth,
    lambda_spatial = lambda_spatial,
    lambda_spatial_weights = lambda_spatial_weights,
    lambda_sparse = lambda_sparse,
    lambda_sparse_weights = lambda_sparse_weights,
    spatial_sigma = spatial_sigma,
    smooth_components = which(lambda_spatial_weights > 0),
    sparse_components = which(lambda_sparse_weights > 0),
    center = if (is.null(cen)) FALSE else cen,
    scale = if (is.null(sc)) FALSE else sc,
    call = match.call()
  )

  class(fit) <- "spax_nmf"
  fit
}

#' Reshape a Spectral Cube Into a Matrix
#'
#' Converts a three-dimensional array with dimensions \code{nx x ny x n_wave}
#' into a matrix suitable for [spax_nmf()], with one row per spaxel.
#'
#' @param cube Numeric 3D array, or a list-like object containing a cube in
#'   \code{$imDat} or \code{$cube}.
#'
#' @return A numeric matrix with \code{nx * ny} rows and \code{n_wave} columns.
#' @examples
#' cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
#' x <- cube_to_matrix(cube$cube)
#' dim(x)
#' @export
cube_to_matrix <- function(cube) {
  metadata <- extract_cube_metadata(cube)
  cube <- resolve_cube(cube)

  if (!is.array(cube) || length(dim(cube)) != 3L) {
    stop("'cube' must be a numeric 3D array.", call. = FALSE)
  }
  if (!is.numeric(cube) || any(!is.finite(cube))) {
    stop("'cube' must contain finite numeric values.", call. = FALSE)
  }

  dims <- dim(cube)
  x <- matrix(cube, nrow = dims[1L] * dims[2L], ncol = dims[3L], byrow = FALSE)
  attr(x, "spax_nmf_metadata") <- metadata
  attr(x, "spectral_unmix_dim") <- dims[1:2]
  x
}

#' Reshape a Matrix Back Into a Cube
#'
#' Reconstructs a \code{nx x ny x n_wave} cube from a matrix with one row per
#' spaxel.
#'
#' @param x Numeric matrix with one row per spatial pixel.
#' @param nx Number of x-axis pixels.
#' @param ny Number of y-axis pixels.
#' @param metadata Optional metadata list to attach to the returned cube.
#'
#' @return A numeric 3D array.
#' @examples
#' cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
#' rebuilt <- matrix_to_cube(cube$matrix, cube$nx, cube$ny)
#' dim(rebuilt)
#' @export
matrix_to_cube <- function(x, nx, ny, metadata = NULL) {
  x <- validate_input_matrix(x)

  if (!is.numeric(nx) || length(nx) != 1L || nx < 1L ||
      !is.numeric(ny) || length(ny) != 1L || ny < 1L) {
    stop("'nx' and 'ny' must be positive integers.", call. = FALSE)
  }

  nx <- as.integer(nx)
  ny <- as.integer(ny)

  if (nrow(x) != nx * ny) {
    stop("nrow(x) must equal nx * ny.", call. = FALSE)
  }

  cube <- array(x, dim = c(nx, ny, ncol(x)))
  cube_metadata <- merge_metadata(extract_matrix_metadata(x), metadata)
  if (length(cube_metadata) > 0L) {
    attr(cube, "spax_nmf_metadata") <- cube_metadata
  }
  cube
}

#' Extract a Component Map
#'
#' Reshapes one component-weight column from a fitted model into an image plane.
#'
#' @param fit Result from [spax_nmf()].
#' @param nx Number of x-axis pixels.
#' @param ny Number of y-axis pixels.
#' @param component Component index to extract.
#'
#' @return A numeric matrix of size \code{nx x ny}.
#' @examples
#' cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
#' fit <- spax_nmf(cube$matrix, k = 3, niter = 50)
#' img <- component_map(fit, cube$nx, cube$ny, 2)
#' dim(img)
#' @export
component_map <- function(fit, nx, ny, component = 1) {
  if (!inherits(fit, "spax_nmf")) {
    stop("'fit' must be a result from spax_nmf().", call. = FALSE)
  }
  if (!is.numeric(nx) || length(nx) != 1L || nx < 1L ||
      !is.numeric(ny) || length(ny) != 1L || ny < 1L) {
    stop("'nx' and 'ny' must be positive integers.", call. = FALSE)
  }
  if (!is.numeric(component) || length(component) != 1L || component < 1L ||
      component > ncol(fit$spatial)) {
    stop("'component' is out of bounds.", call. = FALSE)
  }

  component <- as.integer(component)
  matrix(fit$spatial[, component], nrow = as.integer(nx), ncol = as.integer(ny))
}

#' Extract a Component Spectrum
#'
#' Returns one recovered component spectrum from a fitted model.
#'
#' @param fit Result from [spax_nmf()].
#' @param component Component index to extract.
#'
#' @return A numeric vector.
#' @examples
#' \dontrun{
#' cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
#' fit <- spax_nmf(cube$matrix, k = 3, niter = 50)
#' spec1 <- component_spectrum(fit, 1)
#' }
#' @export
component_spectrum <- function(fit, component = 1) {
  if (!inherits(fit, "spax_nmf")) {
    stop("'fit' must be a result from spax_nmf().", call. = FALSE)
  }
  if (!is.numeric(component) || length(component) != 1L || component < 1L ||
      component > nrow(fit$spectra)) {
    stop("'component' is out of bounds.", call. = FALSE)
  }

  fit$spectra[as.integer(component), ]
}

#' Reconstruct a Single Component
#'
#' Builds the contribution from one component, either as a spaxel-by-wavelength
#' matrix or, when dimensions are supplied, as a spectral cube.
#'
#' @param fit Result from [spax_nmf()].
#' @param component Component index to reconstruct.
#' @param nx Optional x-axis size for cube output.
#' @param ny Optional y-axis size for cube output.
#'
#' @return A numeric matrix, or a numeric 3D array when \code{nx} and
#'   \code{ny} are supplied.
#' @examples
#' \dontrun{
#' cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
#' fit <- spax_nmf(cube$matrix, k = 3, niter = 50)
#' comp1 <- component_reconstruction(fit, 1)
#' comp1_cube <- component_reconstruction(fit, 1, nx = cube$nx, ny = cube$ny)
#' }
#' @export
component_reconstruction <- function(fit, component = 1, nx = NULL, ny = NULL) {
  if (!inherits(fit, "spax_nmf")) {
    stop("'fit' must be a result from spax_nmf().", call. = FALSE)
  }
  if (!is.numeric(component) || length(component) != 1L || component < 1L ||
      component > nrow(fit$spectra)) {
    stop("'component' is out of bounds.", call. = FALSE)
  }

  component <- as.integer(component)
  reconstruction <- fit$spatial[, component, drop = FALSE] %*%
    fit$spectra[component, , drop = FALSE]

  if (is.null(nx) && is.null(ny)) {
    return(reconstruction)
  }
  if (is.null(nx) || is.null(ny)) {
    stop("Supply both 'nx' and 'ny' to return a cube.", call. = FALSE)
  }

  matrix_to_cube(reconstruction, nx = nx, ny = ny, metadata = fit$metadata)
}

#' Simulate a Small IFU-like Cube
#'
#' Generates a synthetic cube with smooth continuum structure and simple
#' emission features. This is intended for examples, demos, and documentation.
#'
#' @param nx Number of x-axis pixels.
#' @param ny Number of y-axis pixels.
#' @param n_wave Number of wavelength channels.
#' @param noise Standard deviation of additive Gaussian noise.
#'
#' @return A list with the simulated \code{cube}, its matrix form, the true
#' component \code{spectra}, component-weight maps \code{abundances}, wavelength grid,
#' and spatial dimensions.
#' @examples
#' demo <- simulate_ifu_cube()
#' str(demo, max.level = 1)
#' @export
simulate_ifu_cube <- function(nx = 20, ny = 20, n_wave = 120, noise = 0.01) {
  if (!is.numeric(nx) || !is.numeric(ny) || !is.numeric(n_wave) || !is.numeric(noise)) {
    stop("All arguments must be numeric.", call. = FALSE)
  }

  nx <- as.integer(nx)
  ny <- as.integer(ny)
  n_wave <- as.integer(n_wave)

  wavelength <- seq(3600, 7200, length.out = n_wave)
  xgrid <- seq(-1, 1, length.out = nx)
  ygrid <- seq(-1, 1, length.out = ny)
  xy <- expand.grid(x = xgrid, y = ygrid)

  g1 <- exp(-((xy$x + 0.35)^2 + (xy$y + 0.1)^2) / 0.10)
  g2 <- exp(-((xy$x - 0.25)^2 + (xy$y - 0.3)^2) / 0.18)
  disk <- exp(-sqrt(xy$x^2 + xy$y^2) / 0.9)

  spectra <- rbind(
    0.8 + 0.6 * (wavelength / max(wavelength))^-1.2,
    0.15 + 0.25 * exp(-0.5 * ((wavelength - 5007) / 40)^2) +
      0.18 * exp(-0.5 * ((wavelength - 6563) / 55)^2),
    0.25 + 0.35 * exp(-0.5 * ((wavelength - 4300) / 180)^2)
  )

  abundances <- cbind(
    0.6 * disk + 0.25 * g1,
    0.8 * g1,
    0.7 * g2
  )

  matrix_data <- abundances %*% spectra
  matrix_data <- matrix_data + stats::rnorm(length(matrix_data), sd = noise)
  matrix_data[matrix_data < 0] <- 0

  list(
    cube = matrix_to_cube(matrix_data, nx = nx, ny = ny),
    matrix = matrix_data,
    spectra = spectra,
    abundances = lapply(seq_len(ncol(abundances)), function(i) {
      matrix(abundances[, i], nrow = nx, ncol = ny)
    }),
    wavelength = wavelength,
    nx = nx,
    ny = ny
  )
}

#' Run the Built-in IFU Demonstration
#'
#' Generates a synthetic cube, fits the model, and optionally plots the result.
#'
#' @param k Number of components.
#' @param nx Number of x-axis pixels in the synthetic cube.
#' @param ny Number of y-axis pixels in the synthetic cube.
#' @param n_wave Number of wavelength channels.
#' @param niter Number of optimization iterations.
#' @param plot Logical; if \code{TRUE}, display a compact base-R summary plot.
#'
#' @return A list with the simulated data and fitted model.
#' @examples
#' \dontrun{
#' demo <- demo_ifu_unmix(plot = TRUE, niter = 150)
#' }
#' @export
demo_ifu_unmix <- function(k = 3, nx = 18, ny = 18, n_wave = 120, niter = 300, plot = interactive()) {
  demo <- simulate_ifu_cube(nx = nx, ny = ny, n_wave = n_wave)
  fit <- spax_nmf(demo$matrix, k = k, niter = niter)

  if (isTRUE(plot)) {
    plot(fit, type = "spectra", wavelength = demo$wavelength)
    plot(fit, type = "maps", nx = demo$nx, ny = demo$ny)
  }

  list(data = demo, fit = fit)
}

#' @export
print.spax_nmf <- function(x, ...) {
  cat("<spax_nmf>\n", sep = "")
  cat(sprintf("  spaxels: %d\n", nrow(x$spatial)))
  cat(sprintf("  wavelength channels: %d\n", ncol(x$spectra)))
  cat(sprintf("  components: %d\n", nrow(x$spectra)))
  cat(sprintf("  final loss: %.6f\n", utils::tail(x$loss, 1L)))
  invisible(x)
}

#' Return Stored Metadata
#'
#' Retrieves metadata carried by matrices, cubes, or fitted
#' \code{"spax_nmf"} objects.
#'
#' @param object Matrix, cube, or fitted object.
#'
#' @return A metadata list, or \code{NULL} when no metadata is present.
#' @examples
#' demo <- simulate_ifu_cube()
#' cube_metadata(demo$cube)
#' @export
cube_metadata <- function(object) {
  if (inherits(object, "spax_nmf")) {
    return(object$metadata)
  }

  metadata <- attr(object, "spax_nmf_metadata", exact = TRUE)
  if (is.null(metadata)) {
    metadata <- attr(object, "spectral_unmix_metadata", exact = TRUE)
  }
  metadata
}

#' Summarize a SpaxNMF Fit
#'
#' @param object Result from [spax_nmf()].
#' @param ... Unused.
#'
#' @return A list with class \code{"summary.spax_nmf"}.
#' @examples
#' fake <- list(
#'   spatial = matrix(runif(12), 4, 3),
#'   spectra = matrix(runif(15), 3, 5),
#'   loss = c(3, 2, 1)
#' )
#' class(fake) <- "spax_nmf"
#' summary(fake)
#' @export
summary.spax_nmf <- function(object, ...) {
  if (!inherits(object, "spax_nmf")) {
    stop("'object' must be a result from spax_nmf().", call. = FALSE)
  }

  summary_object <- list(
    n_spaxels = nrow(object$spatial),
    n_wave = ncol(object$spectra),
    rank = nrow(object$spectra),
    final_loss = utils::tail(object$loss, 1L),
    n_iter = length(object$loss)
  )

  class(summary_object) <- "summary.spax_nmf"
  summary_object
}

#' @export
print.summary.spax_nmf <- function(x, ...) {
  cat("SpaxNMF summary\n")
  cat(sprintf("  spaxels: %d\n", x$n_spaxels))
  cat(sprintf("  wavelength channels: %d\n", x$n_wave))
  cat(sprintf("  components: %d\n", x$rank))
  cat(sprintf("  iterations: %d\n", x$n_iter))
  cat(sprintf("  final loss: %.6f\n", x$final_loss))
  invisible(x)
}

#' NMF-style Basis Matrix
#'
#' Returns the component spectra arranged as wavelength-by-component, similar to
#' a basis matrix in NMF software interfaces.
#'
#' @param object Result from [spax_nmf()].
#' @param ... Unused.
#'
#' @return A numeric matrix with wavelengths in rows and components in columns.
#' @examples
#' fake <- list(spectra = matrix(runif(15), 3, 5))
#' class(fake) <- "spax_nmf"
#' basis(fake)
#' @export
basis <- function(object, ...) {
  UseMethod("basis")
}

#' @export
basis.spax_nmf <- function(object, ...) {
  if (!inherits(object, "spax_nmf")) {
    stop("'object' must be a result from spax_nmf().", call. = FALSE)
  }

  t(object$spectra)
}

#' Extract Component-Weight Coefficients
#'
#' Returns the component weights arranged as component-by-spaxel, similar to the
#' coefficient matrix used by NMF interfaces.
#'
#' @param object Result from [spax_nmf()].
#' @param ... Unused.
#'
#' @return A numeric matrix with components in rows and spaxels in columns.
#' @examples
#' fake <- list(spatial = matrix(runif(12), 4, 3))
#' class(fake) <- "spax_nmf"
#' coef(fake)
#' @export
coef.spax_nmf <- function(object, ...) {
  if (!inherits(object, "spax_nmf")) {
    stop("'object' must be a result from spax_nmf().", call. = FALSE)
  }

  t(object$spatial)
}

#' Return the Fitted Reconstruction
#'
#' @param object Result from [spax_nmf()].
#' @param nx Optional x-axis size for cube output.
#' @param ny Optional y-axis size for cube output.
#' @param ... Unused.
#'
#' @return A numeric matrix, or a numeric cube when \code{nx} and \code{ny} are
#' supplied.
#' @examples
#' fake <- list(reconstruction = matrix(runif(20), 4, 5))
#' class(fake) <- "spax_nmf"
#' fitted(fake)
#' @export
fitted.spax_nmf <- function(object, nx = NULL, ny = NULL, ...) {
  if (!inherits(object, "spax_nmf")) {
    stop("'object' must be a result from spax_nmf().", call. = FALSE)
  }

  if (is.null(nx) && is.null(ny)) {
    return(object$reconstruction)
  }
  if (is.null(nx) || is.null(ny)) {
    stop("Supply both 'nx' and 'ny' to return a cube.", call. = FALSE)
  }

  matrix_to_cube(object$reconstruction, nx = nx, ny = ny, metadata = object$metadata)
}

#' Compute Residuals
#'
#' @param object Result from [spax_nmf()].
#' @param x Numeric matrix originally used for the fit.
#' @param nx Optional x-axis size for cube output.
#' @param ny Optional y-axis size for cube output.
#' @param ... Unused.
#'
#' @return A numeric matrix, or a numeric cube when \code{nx} and \code{ny} are
#' supplied.
#' @examples
#' fake <- list(reconstruction = matrix(runif(20), 4, 5))
#' class(fake) <- "spax_nmf"
#' x <- matrix(runif(20), 4, 5)
#' residuals(fake, x = x)
#' @export
residuals.spax_nmf <- function(object, x, nx = NULL, ny = NULL, ...) {
  if (!inherits(object, "spax_nmf")) {
    stop("'object' must be a result from spax_nmf().", call. = FALSE)
  }

  x <- validate_input_matrix(x)
  if (nrow(x) != nrow(object$reconstruction) || ncol(x) != ncol(object$reconstruction)) {
    stop("'x' must have the same dimensions as object$reconstruction.", call. = FALSE)
  }
  residual <- x - object$reconstruction

  if (is.null(nx) && is.null(ny)) {
    return(residual)
  }
  if (is.null(nx) || is.null(ny)) {
    stop("Supply both 'nx' and 'ny' to return a cube.", call. = FALSE)
  }

  matrix_to_cube(residual, nx = nx, ny = ny, metadata = object$metadata)
}

#' Predict From a Fitted SpaxNMF Model
#'
#' Uses the fitted component spectra to estimate component weights for new
#' data, or returns the in-sample reconstruction when \code{newdata} is
#' omitted.
#'
#' @param object Result from [spax_nmf()].
#' @param newdata Optional new matrix or cube.
#' @param type Prediction target: \code{"reconstruction"}, \code{"spatial"},
#'   or \code{"cube"}.
#' @param nx Optional x-axis size for cube output.
#' @param ny Optional y-axis size for cube output.
#' @param lambda_spatial Non-negative weight for the spatial smoothness penalty
#'   when estimating component weights for new data. If \code{NULL}, the fitted
#'   object's spatial penalties are reused.
#' @param lambda_sparse Non-negative weight for the sparsity penalty when
#'   estimating component weights for new data. If \code{NULL}, the fitted
#'   object's sparsity penalties are reused.
#' @param spatial_sigma Positive kernel width in spaxels for the Gaussian
#'   spatial weighting used during component-weight estimation.
#' @param smooth_components Optional component indices that should receive the
#'   spatial smoothness prior when \code{lambda_spatial} is scalar.
#' @param sparse_components Optional component indices that should receive the
#'   sparsity prior when \code{lambda_sparse} is scalar.
#' @param lr Learning rate used when estimating component weights for new data.
#' @param niter Number of optimization iterations for new data.
#' @param cuda Logical; if \code{TRUE}, request CUDA execution.
#' @param verbose Logical; if \code{TRUE}, emit progress every 100 iterations.
#' @param ... Unused.
#'
#' @return A numeric matrix or cube, depending on \code{type}.
#' @examples
#' \dontrun{
#' demo <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
#' fit <- spax_nmf(demo$matrix, k = 3, niter = 50)
#' fitted_x <- predict(fit)
#' new_weights <- predict(fit, newdata = demo$cube, type = "spatial")
#' }
#' @export
predict.spax_nmf <- function(object,
                                   newdata = NULL,
                                   type = c("reconstruction", "spatial", "cube"),
                                   nx = NULL,
                                   ny = NULL,
                                   lr = 0.05,
                                   niter = 500,
                                   cuda = FALSE,
                                   verbose = FALSE,
                                   lambda_spatial = NULL,
                                   lambda_sparse = NULL,
                                   spatial_sigma = NULL,
                                   smooth_components = NULL,
                                   sparse_components = NULL,
                                   ...) {
  if (!inherits(object, "spax_nmf")) {
    stop("'object' must be a result from spax_nmf().", call. = FALSE)
  }

  type <- match.arg(type)

  if (is.null(newdata)) {
    if (type == "spatial") {
      return(object$spatial)
    }
    if (type == "reconstruction") {
      return(object$reconstruction)
    }
    if (is.null(nx) || is.null(ny)) {
      stop("'nx' and 'ny' are required when type = 'cube'.", call. = FALSE)
    }
    return(matrix_to_cube(object$reconstruction, nx = nx, ny = ny, metadata = object$metadata))
  }

  prediction_metadata <- extract_matrix_metadata(newdata)
  cube_dims <- infer_cube_dims(newdata)
  if (!is.null(cube_dims)) {
    nx <- cube_dims[1L]
    ny <- cube_dims[2L]
    newdata <- cube_to_matrix(newdata)
    prediction_metadata <- merge_metadata(prediction_metadata, extract_matrix_metadata(newdata))
  }

  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("The 'torch' package must be installed to use predict().", call. = FALSE)
  }

  newdata <- validate_input_matrix(newdata)
  if (ncol(newdata) != ncol(object$spectra)) {
    stop("newdata must have the same number of columns as the fitted spectra.", call. = FALSE)
  }
  if (!is.numeric(lr) || length(lr) != 1L || is.na(lr) || lr <= 0) {
    stop("'lr' must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(niter) || length(niter) != 1L || is.na(niter) || niter < 1L) {
    stop("'niter' must be a positive integer.", call. = FALSE)
  }
  if (is.null(spatial_sigma)) {
    spatial_sigma <- object$spatial_sigma
    if (is.null(spatial_sigma)) {
      spatial_sigma <- 1
    }
  }
  if (!is.numeric(spatial_sigma) || length(spatial_sigma) != 1L || is.na(spatial_sigma) ||
      spatial_sigma <= 0) {
    stop("'spatial_sigma' must be a positive number.", call. = FALSE)
  }

  if (is.null(lambda_spatial)) {
    if (!is.null(object$lambda_spatial_weights)) {
      lambda_spatial_weights <- object$lambda_spatial_weights
    } else {
      lambda_spatial_weights <- resolve_component_penalty(
        object$lambda_spatial,
        k = nrow(object$spectra),
        name = "lambda_spatial"
      )
    }
  } else {
    lambda_spatial_weights <- resolve_component_penalty(
      lambda_spatial,
      k = nrow(object$spectra),
      active_components = smooth_components,
      name = "lambda_spatial"
    )
  }

  if (is.null(lambda_sparse)) {
    if (!is.null(object$lambda_sparse_weights)) {
      lambda_sparse_weights <- object$lambda_sparse_weights
    } else {
      lambda_sparse_weights <- rep(0, nrow(object$spectra))
    }
  } else {
    lambda_sparse_weights <- resolve_component_penalty(
      lambda_sparse,
      k = nrow(object$spectra),
      active_components = sparse_components,
      name = "lambda_sparse"
    )
  }

  spatial_dims <- resolve_spatial_dims(newdata, nx = nx, ny = ny, lambda_spatial = lambda_spatial_weights)
  spatial_graph <- build_spatial_graph(
    spatial_dims,
    n_expected = nrow(newdata),
    spatial_sigma = spatial_sigma
  )

  device <- select_torch_device(cuda)
  X <- torch::torch_tensor(as.matrix(newdata), device = device)
  S_fixed <- torch::torch_tensor(object$spectra, device = device)
  A <- torch::torch_rand(nrow(newdata), nrow(object$spectra), device = device, requires_grad = TRUE)
  opt <- torch::optim_adam(list(A), lr = lr)

  for (i in seq_len(as.integer(niter))) {
    opt$zero_grad()
    Xhat <- A$matmul(S_fixed)
    recon <- torch::torch_mean((X - Xhat)^2)
    spatial_smooth <- spatial_penalty(
      A,
      spatial_graph,
      device = device,
      component_weights = lambda_spatial_weights
    )
    sparse_abundance <- sparse_penalty(
      A,
      device = device,
      component_weights = lambda_sparse_weights
    )
    loss <- recon + spatial_smooth + sparse_abundance
    loss$backward()
    opt$step()
    A$data()$clamp_(min = 0)

    if (verbose && (i %% 100L == 0L || i == 1L || i == niter)) {
      message(sprintf("Prediction iteration %d/%d, loss = %.6f", i, niter, as.numeric(loss$item())))
    }
  }

  spatial <- as.matrix(A$detach()$cpu())
  reconstruction <- spatial %*% object$spectra

  if (type == "spatial") {
    return(spatial)
  }
  if (type == "reconstruction") {
    return(reconstruction)
  }

  if (is.null(nx) || is.null(ny)) {
    stop("'nx' and 'ny' are required when type = 'cube'.", call. = FALSE)
  }

  matrix_to_cube(
    reconstruction,
    nx = nx,
    ny = ny,
    metadata = merge_metadata(prediction_metadata, object$metadata)
  )
}

#' Plot a Spectral Unmixing Result
#'
#' Convenience plots for component spectra, spatial maps, and optimization loss.
#'
#' @param x Result from [spax_nmf()].
#' @param type Plot type: \code{"spectra"}, \code{"maps"}, or \code{"loss"}.
#' @param nx Number of x-axis pixels when plotting maps.
#' @param ny Number of y-axis pixels when plotting maps.
#' @param components Optional vector of component indices to display.
#' @param wavelength Optional wavelength vector for spectra.
#' @param ... Additional graphical arguments passed to base plotting functions.
#'
#' @return Invisibly returns \code{x}.
#' @examples
#' \dontrun{
#' cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
#' fit <- spax_nmf(cube$matrix, k = 3, niter = 50)
#' plot(fit, type = "spectra", wavelength = cube$wavelength)
#' plot(fit, type = "maps", nx = cube$nx, ny = cube$ny)
#' plot(fit, type = "loss")
#' }
#' @export
plot.spax_nmf <- function(x,
                                type = c("spectra", "maps", "loss"),
                                nx = NULL,
                                ny = NULL,
                                components = NULL,
                                wavelength = NULL,
                                ...) {
  if (!inherits(x, "spax_nmf")) {
    stop("'x' must be a result from spax_nmf().", call. = FALSE)
  }

  type <- match.arg(type)

  if (is.null(components)) {
    components <- seq_len(nrow(x$spectra))
  }
  components <- as.integer(components)
  if (any(is.na(components)) || any(components < 1L) || any(components > nrow(x$spectra))) {
    stop("'components' contains invalid indices.", call. = FALSE)
  }

  if (type == "loss") {
    graphics::plot(
      seq_along(x$loss),
      x$loss,
      type = "l",
      xlab = "Iteration",
      ylab = "Loss",
      main = "Optimization history",
      ...
    )
    return(invisible(x))
  }

  if (type == "spectra") {
    if (is.null(wavelength)) {
      wavelength <- seq_len(ncol(x$spectra))
      xlab <- "Wavelength channel"
    } else {
      xlab <- "Wavelength"
    }

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(length(components), 1L), mar = c(3, 4, 2, 1))

    for (component in components) {
      graphics::plot(
        wavelength,
        component_spectrum(x, component),
        type = "l",
        lwd = 2,
        xlab = xlab,
        ylab = "Flux",
        main = sprintf("Component spectrum %d", component),
        ...
      )
    }

    return(invisible(x))
  }

  if (is.null(nx) || is.null(ny)) {
    stop("'nx' and 'ny' are required when type = 'maps'.", call. = FALSE)
  }

  n_panels <- length(components)
  n_col <- min(3L, n_panels)
  n_row <- ceiling(n_panels / n_col)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(n_row, n_col), mar = c(3, 3, 2, 1))

  for (component in components) {
    graphics::image(
      component_map(x, nx = nx, ny = ny, component = component),
      xlab = "x",
      ylab = "y",
      main = sprintf("Component map %d", component),
      col = grDevices::hcl.colors(32, "YlOrRd", rev = TRUE),
      ...
    )
  }

  invisible(x)
}

#' Plot Data Versus Reconstruction
#'
#' Plots observed and reconstructed spectra for selected spaxels.
#'
#' @param fit Result from [spax_nmf()].
#' @param x Numeric matrix used as input to the fit.
#' @param pixels Optional vector of spaxel indices to plot.
#' @param n Number of random spaxels to plot when \code{pixels} is omitted.
#' @param wavelength Optional wavelength vector.
#' @param ... Additional graphical arguments passed to \code{plot()}.
#'
#' @return Invisibly returns \code{fit}.
#' @examples
#' \dontrun{
#' cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
#' fit <- spax_nmf(cube$matrix, k = 3, niter = 50)
#' plot_reconstruction(fit, cube$matrix, n = 4, wavelength = cube$wavelength)
#' }
#' @export
plot_reconstruction <- function(fit, x, pixels = NULL, n = 6, wavelength = NULL, ...) {
  if (!inherits(fit, "spax_nmf")) {
    stop("'fit' must be a result from spax_nmf().", call. = FALSE)
  }

  x <- validate_input_matrix(x)

  if (nrow(x) != nrow(fit$reconstruction) || ncol(x) != ncol(fit$reconstruction)) {
    stop("'x' must have the same dimensions as fit$reconstruction.", call. = FALSE)
  }

  if (is.null(pixels)) {
    n <- min(as.integer(n), nrow(x))
    pixels <- sort(sample.int(nrow(x), n))
  } else {
    pixels <- as.integer(pixels)
  }
  if (any(is.na(pixels)) || any(pixels < 1L) || any(pixels > nrow(x))) {
    stop("'pixels' contains invalid indices.", call. = FALSE)
  }

  if (is.null(wavelength)) {
    wavelength <- seq_len(ncol(x))
    xlab <- "Wavelength channel"
  } else {
    xlab <- "Wavelength"
  }

  n_panels <- length(pixels)
  n_col <- min(2L, n_panels)
  n_row <- ceiling(n_panels / n_col)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(n_row, n_col), mar = c(3, 4, 2, 1))

  for (pixel in pixels) {
    graphics::plot(
      wavelength,
      x[pixel, ],
      type = "l",
      lwd = 2,
      col = "black",
      lty = 1,
      xlab = xlab,
      ylab = "Flux",
      main = sprintf("Spaxel %d", pixel),
      ...
    )
    graphics::lines(
      wavelength,
      fit$reconstruction[pixel, ],
      col = "red",
      lwd = 2,
      lty = 2
    )
    graphics::legend(
      "topright",
      legend = c("data", "reconstruction"),
      col = c("black", "red"),
      lty = c(1, 2),
      bty = "n"
    )
  }

  invisible(fit)
}

select_torch_device <- function(cuda) {
  if (!isTRUE(cuda)) {
    return(torch::torch_device("cpu"))
  }

  if (!torch::cuda_is_available()) {
    stop("CUDA was requested but torch reports that CUDA is not available.", call. = FALSE)
  }

  torch::torch_device("cuda:0")
}

validate_input_matrix <- function(x) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("'x' must be a numeric matrix or data frame.", call. = FALSE)
  }
  if (nrow(x) < 2L || ncol(x) < 2L) {
    stop("'x' must have at least two rows and two columns.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("'x' must contain only finite numeric values.", call. = FALSE)
  }

  x
}

resolve_cube <- function(cube) {
  if (is.list(cube) && !is.null(cube$imDat)) {
    return(cube$imDat)
  }
  if (is.list(cube) && !is.null(cube$cube)) {
    return(cube$cube)
  }

  cube
}

infer_cube_dims <- function(x) {
  cube <- resolve_cube(x)
  if (!is.array(cube) || length(dim(cube)) != 3L) {
    return(NULL)
  }

  dim(cube)[1:2]
}

extract_cube_metadata <- function(x) {
  if (inherits(x, "spax_nmf")) {
    return(x$metadata)
  }

  metadata <- attr(x, "spax_nmf_metadata", exact = TRUE)
  if (is.null(metadata)) {
    metadata <- attr(x, "spectral_unmix_metadata", exact = TRUE)
  }
  if (!is.null(metadata)) {
    return(metadata)
  }

  if (is.list(x)) {
    keep <- setdiff(names(x), c("imDat", "cube"))
    if (length(keep) > 0L) {
      return(x[keep])
    }
  }

  NULL
}

extract_matrix_metadata <- function(x) {
  metadata <- attr(x, "spax_nmf_metadata", exact = TRUE)
  if (is.null(metadata)) {
    metadata <- attr(x, "spectral_unmix_metadata", exact = TRUE)
  }
  dims <- attr(x, "spectral_unmix_dim", exact = TRUE)

  if (!is.null(dims)) {
    metadata <- merge_metadata(metadata, list(dim_xy = dims))
  }

  metadata
}

merge_metadata <- function(x, y) {
  if (is.null(x) || length(x) == 0L) {
    return(y)
  }
  if (is.null(y) || length(y) == 0L) {
    return(x)
  }

  utils::modifyList(x, y)
}

resolve_spatial_dims <- function(x, nx = NULL, ny = NULL, lambda_spatial = 0) {
  if (!is.null(nx) || !is.null(ny)) {
    if (is.null(nx) || is.null(ny)) {
      stop("Supply both 'nx' and 'ny' when providing spatial dimensions.", call. = FALSE)
    }
    if (!is.numeric(nx) || length(nx) != 1L || nx < 1L ||
        !is.numeric(ny) || length(ny) != 1L || ny < 1L) {
      stop("'nx' and 'ny' must be positive integers.", call. = FALSE)
    }
    return(c(as.integer(nx), as.integer(ny)))
  }

  dims <- attr(x, "spectral_unmix_dim", exact = TRUE)
  if (!is.null(dims)) {
    return(as.integer(dims))
  }

  metadata <- extract_matrix_metadata(x)
  if (!is.null(metadata$dim_xy)) {
    return(as.integer(metadata$dim_xy))
  }

  if (has_positive_penalty(lambda_spatial)) {
    stop(
      "Spatial smoothing requires grid dimensions. Supply 'nx' and 'ny', ",
      "or fit from an object produced by cube_to_matrix().",
      call. = FALSE
    )
  }

  NULL
}

build_spatial_graph <- function(spatial_dims, n_expected = NULL, spatial_sigma = 1) {
  if (is.null(spatial_dims)) {
    return(NULL)
  }

  nx <- as.integer(spatial_dims[1L])
  ny <- as.integer(spatial_dims[2L])
  if (!is.null(n_expected) && nx * ny != n_expected) {
    stop("Spatial dimensions do not match the number of spectra.", call. = FALSE)
  }
  if (!is.numeric(spatial_sigma) || length(spatial_sigma) != 1L || is.na(spatial_sigma) ||
      spatial_sigma <= 0) {
    stop("'spatial_sigma' must be a positive number.", call. = FALSE)
  }

  ids <- matrix(seq_len(nx * ny), nrow = nx, ncol = ny)
  radius <- max(1L, as.integer(ceiling(3 * spatial_sigma)))
  edges_i <- integer()
  edges_j <- integer()
  weights <- numeric()

  for (dx in 0:radius) {
    for (dy in (-radius):radius) {
      if (dx == 0L && dy <= 0L) {
        next
      }
      if (dx >= nx || abs(dy) >= ny) {
        next
      }

      dist <- sqrt(dx^2 + dy^2)
      if (dist > radius) {
        next
      }

      x_from <- seq_len(nx - dx)
      x_to <- x_from + dx

      if (dy >= 0L) {
        y_from <- seq_len(ny - dy)
      } else {
        y_from <- seq.int(1L - dy, ny)
      }
      y_to <- y_from + dy

      edge_i <- as.vector(ids[x_from, y_from, drop = FALSE])
      edge_j <- as.vector(ids[x_to, y_to, drop = FALSE])

      edges_i <- c(edges_i, edge_i)
      edges_j <- c(edges_j, edge_j)
      weights <- c(weights, rep(exp(-0.5 * (dist / spatial_sigma)^2), length(edge_i)))
    }
  }

  if (!length(weights)) {
    return(NULL)
  }

  weights <- weights / mean(weights)

  list(
    edges = cbind(as.integer(edges_i), as.integer(edges_j)),
    weights = as.numeric(weights),
    sigma = spatial_sigma
  )
}

spatial_penalty <- function(A, spatial_graph, device, component_weights = NULL) {
  if (is.null(spatial_graph)) {
    return(torch::torch_tensor(0, device = device))
  }

  edges <- spatial_graph$edges
  if (is.null(edges) || nrow(edges) == 0L) {
    return(torch::torch_tensor(0, device = device))
  }

  diff <- A[edges[, 1L], , drop = FALSE] - A[edges[, 2L], , drop = FALSE]
  edge_weights <- torch::torch_tensor(
    matrix(spatial_graph$weights, ncol = 1L),
    device = device
  )
  weighted_diff <- edge_weights * (diff^2)

  if (is.null(component_weights)) {
    return(torch::torch_mean(weighted_diff))
  }

  component_weights <- as.numeric(component_weights)
  if (!length(component_weights) || !any(component_weights > 0)) {
    return(torch::torch_tensor(0, device = device))
  }

  comp_weights <- torch::torch_tensor(matrix(component_weights, ncol = 1L), device = device)
  torch::torch_mean(weighted_diff$matmul(comp_weights))
}

sparse_penalty <- function(A, device, component_weights = NULL) {
  if (is.null(component_weights)) {
    return(torch::torch_tensor(0, device = device))
  }

  component_weights <- as.numeric(component_weights)
  if (!length(component_weights) || !any(component_weights > 0)) {
    return(torch::torch_tensor(0, device = device))
  }

  comp_weights <- torch::torch_tensor(matrix(component_weights, ncol = 1L), device = device)
  torch::torch_mean(A$matmul(comp_weights))
}

has_positive_penalty <- function(x) {
  is.numeric(x) && length(x) > 0L && any(is.finite(x) & x > 0)
}

normalize_component_indices <- function(components, k, name) {
  if (is.null(components)) {
    return(NULL)
  }

  components <- as.integer(components)
  components <- unique(components[!is.na(components)])
  if (!length(components)) {
    return(integer())
  }
  if (any(components < 1L | components > k)) {
    stop(sprintf("'%s' contains invalid component indices.", name), call. = FALSE)
  }

  sort(components)
}

resolve_component_penalty <- function(value, k, active_components = NULL, name = "penalty") {
  if (!is.numeric(value) || !length(value) || any(is.na(value)) || any(value < 0)) {
    stop(sprintf("'%s' must be a non-negative numeric scalar or length-k vector.", name), call. = FALSE)
  }

  active_components <- normalize_component_indices(
    active_components,
    k = k,
    name = paste0(name, "_components")
  )

  if (length(value) == 1L) {
    weights <- rep(0, k)
    if (!value) {
      return(weights)
    }

    if (is.null(active_components)) {
      active_components <- seq_len(k)
    }
    if (!length(active_components)) {
      return(weights)
    }

    weights[active_components] <- as.numeric(value) / length(active_components)
    return(weights)
  }

  if (length(value) != k) {
    stop(sprintf("'%s' must have length 1 or k = %d.", name, k), call. = FALSE)
  }
  if (!is.null(active_components)) {
    stop(
      sprintf(
        "Supply either a scalar '%s' with selected components or a length-k vector, not both.",
        name
      ),
      call. = FALSE
    )
  }

  as.numeric(value)
}
