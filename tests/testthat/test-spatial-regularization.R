if (!exists("spectral_unmix", mode = "function")) {
  candidate_paths <- c(
    file.path(getwd(), "R", "spectral_unmix.R"),
    file.path(getwd(), "..", "..", "R", "spectral_unmix.R")
  )
  source_path <- candidate_paths[file.exists(candidate_paths)][1]
  if (is.na(source_path)) {
    stop("Could not locate R/spectral_unmix.R for source-based tests.", call. = FALSE)
  }
  sys.source(source_path, envir = environment())
}

set_all_seeds <- function(seed) {
  set.seed(seed)
  torch::torch_manual_seed(seed)
}

normalize_map <- function(x) {
  scale <- max(x, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) {
    return(matrix(0, nrow = nrow(x), ncol = ncol(x)))
  }

  x / scale
}

mean_map_roughness <- function(maps) {
  mean(vapply(maps, function(map) {
    map <- normalize_map(map)
    dx <- if (nrow(map) > 1L) diff(map, differences = 1L, lag = 1L) else 0
    dy <- if (ncol(map) > 1L) t(diff(t(map), differences = 1L, lag = 1L)) else 0
    mean(dx^2) + mean(dy^2)
  }, numeric(1)))
}

fit_maps <- function(fit, nx, ny) {
  lapply(seq_len(ncol(fit$spatial)), function(j) {
    matrix(fit$spatial[, j], nrow = nx, ncol = ny)
  })
}

normalize_spectrum <- function(x) {
  scale <- max(x, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) {
    return(rep(0, length(x)))
  }

  x / scale
}

best_spectral_permutation <- function(fit, truth) {
  k <- nrow(truth)
  perms <- expand.grid(rep(list(seq_len(k)), k))
  perms <- as.matrix(perms)
  perms <- perms[apply(perms, 1, function(z) length(unique(z)) == k), , drop = FALSE]

  fit_norm <- t(apply(fit$spectra, 1, normalize_spectrum))
  truth_norm <- t(apply(truth, 1, normalize_spectrum))

  scores <- apply(perms, 1, function(perm) {
    sum(vapply(seq_len(k), function(i) {
      mean((fit_norm[perm[i], ] - truth_norm[i, ])^2)
    }, numeric(1)))
  })

  perms[which.min(scores), ]
}

reorder_fit <- function(fit, perm) {
  fit$spatial <- fit$spatial[, perm, drop = FALSE]
  fit$abundance <- fit$spatial
  fit$spectra <- fit$spectra[perm, , drop = FALSE]
  fit$basis <- t(fit$spectra)
  fit$coef <- t(fit$spatial)
  fit$reconstruction <- fit$spatial %*% fit$spectra
  fit$fitted <- fit$reconstruction
  fit
}

true_maps <- function(demo) {
  lapply(demo$abundances, as.matrix)
}

simulate_realistic_ifu_cube <- function(nx = 10, ny = 10, n_wave = 60, noise = 0.03) {
  wave <- seq(3600, 9200, length.out = n_wave)

  absorption_line <- function(center, width, depth) {
    depth * exp(-0.5 * ((wave - center) / width)^2)
  }

  emission_line <- function(center, width, amp) {
    amp * exp(-0.5 * ((wave - center) / width)^2)
  }

  hot_cont <- 1.2 * (wave / 4000)^-1.35
  hot <- hot_cont -
    absorption_line(3934, 25, 0.18) -
    absorption_line(4102, 35, 0.14) -
    absorption_line(4341, 40, 0.16) -
    absorption_line(4861, 45, 0.12)

  old_cont <- 0.85 * (wave / 5500)^-0.25 + 0.06 * exp(-0.5 * ((wave - 6000) / 850)^2)
  old <- old_cont -
    absorption_line(3934, 28, 0.08) -
    absorption_line(3968, 30, 0.07) -
    absorption_line(4305, 55, 0.12) -
    absorption_line(5175, 65, 0.10) -
    absorption_line(5892, 45, 0.08)

  nebular <- 0.03 +
    emission_line(3727, 18, 0.35) +
    emission_line(4861, 18, 0.28) +
    emission_line(4959, 15, 0.24) +
    emission_line(5007, 15, 0.42) +
    emission_line(6563, 18, 0.34) +
    emission_line(6583, 16, 0.12) +
    emission_line(6716, 15, 0.09) +
    emission_line(6731, 15, 0.09)

  spectra <- rbind(hot, old, nebular)
  spectra <- pmax(spectra, 1e-4)
  spectra <- spectra / apply(spectra, 1, max)

  xgrid <- seq(-1, 1, length.out = nx)
  ygrid <- seq(-1, 1, length.out = ny)
  xy <- expand.grid(x = xgrid, y = ygrid)

  bulge <- exp(-((xy$x + 0.05)^2 + (xy$y + 0.05)^2) / 0.22)
  arm_a <- exp(-((xy$x + 0.42)^2 + (xy$y - 0.25)^2) / 0.08)
  arm_b <- exp(-((xy$x - 0.35)^2 + (xy$y + 0.30)^2) / 0.10)
  ring <- exp(-((sqrt(xy$x^2 + xy$y^2) - 0.45)^2) / 0.030)

  abundances <- cbind(
    0.70 * bulge + 0.20 * ring,
    0.45 * bulge + 0.30 * arm_b,
    0.80 * arm_a + 0.55 * ring
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
    wavelength = wave,
    nx = nx,
    ny = ny
  )
}

simulate_compact_ifu_cube <- function(nx = 12, ny = 12, n_wave = 70, noise = 0.015) {
  wave <- seq(4200, 8800, length.out = n_wave)

  gaussian_line <- function(center, width, amp) {
    amp * exp(-0.5 * ((wave - center) / width)^2)
  }

  old_pop <- 0.95 * (wave / 5500)^-0.20 -
    gaussian_line(4305, 60, 0.08) -
    gaussian_line(5175, 70, 0.09) -
    gaussian_line(5892, 45, 0.06)

  disk_pop <- 0.75 * (wave / 5000)^-0.55 +
    gaussian_line(6563, 22, 0.10) +
    gaussian_line(5007, 18, 0.08) +
    gaussian_line(3727, 16, 0.06)

  gc_pop <- 1.20 * (wave / 4300)^-1.00 -
    gaussian_line(4341, 35, 0.14) -
    gaussian_line(4861, 38, 0.10) +
    gaussian_line(5007, 16, 0.04)

  spectra <- rbind(old_pop, disk_pop, pmax(gc_pop, 1e-4))
  spectra <- pmax(spectra, 1e-4)
  spectra <- spectra / apply(spectra, 1, max)

  xy <- expand.grid(
    x = seq(-1, 1, length.out = nx),
    y = seq(-1, 1, length.out = ny)
  )

  bulge <- exp(-((xy$x)^2 + (xy$y)^2) / 0.18)
  disk <- exp(-((sqrt((xy$x + 0.15)^2 + (xy$y - 0.05)^2) - 0.45)^2) / 0.03)

  gc_map <- matrix(0, nrow = nx, ncol = ny)
  gc_centers <- rbind(c(3, 9), c(8, 4), c(10, 10))
  for (ii in seq_len(nrow(gc_centers))) {
    gc_map[gc_centers[ii, 1], gc_centers[ii, 2]] <- gc_map[gc_centers[ii, 1], gc_centers[ii, 2]] + 1
  }

  abundances <- cbind(
    0.90 * bulge,
    0.55 * disk + 0.20 * bulge,
    as.vector(gc_map)
  )

  matrix_data <- abundances %*% spectra
  matrix_data <- matrix_data + stats::rnorm(length(matrix_data), sd = noise)
  matrix_data[matrix_data < 0] <- 0

  list(
    matrix = matrix_data,
    spectra = spectra,
    compact_map = gc_map / sum(gc_map),
    nx = nx,
    ny = ny
  )
}

top_flux_fraction <- function(map, n = 4L) {
  x <- sort(as.numeric(map), decreasing = TRUE)
  denom <- sum(x)
  if (!is.finite(denom) || denom <= 0) {
    return(0)
  }

  sum(x[seq_len(min(n, length(x)))]) / denom
}

test_that("lambda_spatial = 0 reduces to the original smooth torch NMF objective", {
  skip_if_not_installed("torch")

  demo <- simulate_realistic_ifu_cube(nx = 7, ny = 6, n_wave = 40, noise = 0.01)

  set_all_seeds(11)
  fit_plain <- spectral_unmix(
    demo$matrix,
    k = 3,
    lambda_smooth = 0.02,
    niter = 60,
    lr = 0.03
  )

  set_all_seeds(11)
  fit_spatial_zero <- spectral_unmix(
    demo$matrix,
    k = 3,
    lambda_smooth = 0.02,
    lambda_spatial = 0,
    spatial_sigma = 1,
    nx = demo$nx,
    ny = demo$ny,
    niter = 60,
    lr = 0.03
  )

  expect_equal(fit_plain$spatial, fit_spatial_zero$spatial, tolerance = 1e-6)
  expect_equal(fit_plain$spectra, fit_spatial_zero$spectra, tolerance = 1e-6)
  expect_equal(fit_plain$reconstruction, fit_spatial_zero$reconstruction, tolerance = 1e-6)
  expect_equal(fit_plain$loss, fit_spatial_zero$loss, tolerance = 1e-6)
})

test_that("spatial regularization improves smooth-map recovery on smooth synthetic IFU data", {
  skip_if_not_installed("torch")

  demo <- simulate_realistic_ifu_cube(nx = 10, ny = 10, n_wave = 60, noise = 0.03)

  set_all_seeds(23)
  fit_vanilla <- spectral_unmix(
    demo$matrix,
    k = 3,
    lambda_smooth = 0.01,
    lambda_spatial = 0,
    niter = 100,
    lr = 0.03
  )

  set_all_seeds(23)
  fit_spatial <- spectral_unmix(
    demo$matrix,
    k = 3,
    lambda_smooth = 0.01,
    lambda_spatial = 0.01,
    spatial_sigma = 1,
    nx = demo$nx,
    ny = demo$ny,
    niter = 120,
    lr = 0.03
  )

  rough_true <- mean_map_roughness(true_maps(demo))
  rough_vanilla <- mean_map_roughness(fit_maps(fit_vanilla, demo$nx, demo$ny))
  rough_spatial <- mean_map_roughness(fit_maps(fit_spatial, demo$nx, demo$ny))

  mse_vanilla <- mean((demo$matrix - fit_vanilla$reconstruction)^2)
  mse_spatial <- mean((demo$matrix - fit_spatial$reconstruction)^2)

  expect_lt(abs(rough_spatial - rough_true), abs(rough_vanilla - rough_true))
  expect_lte(mse_spatial, mse_vanilla * 1.25)
})

test_that("sparse component priors improve compact-source recovery", {
  skip_if_not_installed("torch")

  demo <- simulate_compact_ifu_cube()

  set_all_seeds(101)
  fit_all_smooth <- spectral_unmix(
    demo$matrix,
    k = 3,
    lambda_smooth = 0.01,
    lambda_spatial = 0.04,
    nx = demo$nx,
    ny = demo$ny,
    niter = 180,
    lr = 0.03
  )

  set_all_seeds(101)
  fit_mixed <- spectral_unmix(
    demo$matrix,
    k = 3,
    lambda_smooth = 0.01,
    lambda_spatial = 0.04,
    smooth_components = 1:2,
    lambda_sparse = 0.02,
    sparse_components = 3,
    nx = demo$nx,
    ny = demo$ny,
    niter = 180,
    lr = 0.03
  )

  fit_all_smooth <- reorder_fit(fit_all_smooth, best_spectral_permutation(fit_all_smooth, demo$spectra))
  fit_mixed <- reorder_fit(fit_mixed, best_spectral_permutation(fit_mixed, demo$spectra))

  compact_smooth <- matrix(fit_all_smooth$spatial[, 3], nrow = demo$nx, ncol = demo$ny)
  compact_mixed <- matrix(fit_mixed$spatial[, 3], nrow = demo$nx, ncol = demo$ny)

  truth_frac <- top_flux_fraction(demo$compact_map)
  smooth_frac <- top_flux_fraction(compact_smooth)
  mixed_frac <- top_flux_fraction(compact_mixed)

  mse_smooth <- mean((demo$matrix - fit_all_smooth$reconstruction)^2)
  mse_mixed <- mean((demo$matrix - fit_mixed$reconstruction)^2)

  expect_lt(abs(mixed_frac - truth_frac), abs(smooth_frac - truth_frac))
  expect_gt(mixed_frac, smooth_frac)
  expect_lte(mse_mixed, mse_smooth * 1.20)
})
