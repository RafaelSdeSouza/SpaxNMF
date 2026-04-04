if (!exists("spax_nmf", mode = "function")) {
  candidate_paths <- c(
    file.path(getwd(), "R", "spax_nmf.R"),
    file.path(getwd(), "..", "..", "R", "spax_nmf.R")
  )
  source_path <- candidate_paths[file.exists(candidate_paths)][1]
  if (is.na(source_path)) {
    stop("Could not locate R/spax_nmf.R for source-based tests.", call. = FALSE)
  }
  sys.source(source_path, envir = environment())
}

set_all_seeds <- function(seed) {
  set.seed(seed)
  if (requireNamespace("torch", quietly = TRUE)) {
    torch::torch_manual_seed(seed)
  }
}

display_positive <- function(x) {
  s <- max(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) {
    return(rep(0, length(x)))
  }

  x / s
}

best_assignment <- function(dist_mat) {
  k <- nrow(dist_mat)
  perms <- expand.grid(rep(list(seq_len(k)), k))
  perms <- as.matrix(perms)
  perms <- perms[apply(perms, 1, function(z) length(unique(z)) == k), , drop = FALSE]

  scores <- apply(perms, 1, function(p) {
    sum(vapply(seq_len(k), function(i) dist_mat[i, p[i]], numeric(1)))
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

spectral_distance <- function(truth, fit) {
  outer(
    seq_len(nrow(truth)),
    seq_len(nrow(fit)),
    Vectorize(function(i, j) {
      mean((display_positive(truth[i, ]) - display_positive(fit[j, ]))^2)
    })
  )
}

fit_best_spax <- function(x, seeds, ...) {
  best_fit <- NULL
  best_mse <- Inf

  for (seed in seeds) {
    set_all_seeds(seed)
    fit <- spax_nmf(x, ...)
    mse <- mean((x - fit$reconstruction)^2)
    if (mse < best_mse) {
      best_fit <- fit
      best_mse <- mse
    }
  }

  best_fit
}

extract_reference_fit <- function(fit) {
  spatial <- as.matrix(NMF::basis(fit))
  spectra <- as.matrix(NMF::coef(fit))

  list(
    spatial = spatial,
    abundance = spatial,
    spectra = spectra,
    basis = t(spectra),
    coef = t(spatial),
    reconstruction = spatial %*% spectra,
    fitted = spatial %*% spectra
  )
}

fit_best_reference <- function(x, rank, seeds) {
  best_fit <- NULL
  best_mse <- Inf

  for (seed in seeds) {
    set.seed(seed)
    fit_raw <- suppressMessages(suppressWarnings(
      NMF::nmf(x, rank = rank, method = "lee", nrun = 1, seed = seed)
    ))
    fit <- extract_reference_fit(fit_raw)
    mse <- mean((x - fit$reconstruction)^2)
    if (mse < best_mse) {
      best_fit <- fit
      best_mse <- mse
    }
  }

  best_fit
}

spectra_mse <- function(truth, fit) {
  mean((t(apply(truth, 1, display_positive)) - t(apply(fit, 1, display_positive)))^2)
}

test_that("standalone spax_nmf recovers an exact low-rank toy factorization", {
  toy <- simulate_ifu_cube(nx = 8, ny = 8, n_wave = 60, noise = 0)
  fit_seeds <- c(11, 23, 37, 41, 53, 67)

  fit <- fit_best_spax(
    toy$matrix,
    seeds = fit_seeds,
    k = 3,
    lambda_smooth = 0,
    lambda_spatial = 0,
    lambda_sparse = 0,
    solver = "adam",
    niter = 2000,
    lr = 0.03,
    tol = 0
  )

  fit <- reorder_fit(fit, best_assignment(spectral_distance(toy$spectra, fit$spectra)))

  recon_mse <- mean((toy$matrix - fit$reconstruction)^2)
  spec_mse <- spectra_mse(toy$spectra, fit$spectra)

  expect_identical(fit$solver, "adam")
  expect_lt(recon_mse, 1e-5)
  expect_lt(spec_mse, 5e-2)
})

test_that("standalone spax_nmf is competitive with NMF package on the same toy problem", {
  skip_if_not_installed("NMF")

  toy <- simulate_ifu_cube(nx = 8, ny = 8, n_wave = 60, noise = 0)
  fit_seeds <- c(11, 23, 37, 41, 53, 67)

  fit_spax <- fit_best_spax(
    toy$matrix,
    seeds = fit_seeds,
    k = 3,
    lambda_smooth = 0,
    lambda_spatial = 0,
    lambda_sparse = 0,
    solver = "adam",
    niter = 2000,
    lr = 0.03,
    tol = 0
  )
  fit_ref <- fit_best_reference(toy$matrix, rank = 3, seeds = fit_seeds)

  fit_spax <- reorder_fit(fit_spax, best_assignment(spectral_distance(toy$spectra, fit_spax$spectra)))
  fit_ref <- reorder_fit(fit_ref, best_assignment(spectral_distance(toy$spectra, fit_ref$spectra)))

  recon_spax <- mean((toy$matrix - fit_spax$reconstruction)^2)
  recon_ref <- mean((toy$matrix - fit_ref$reconstruction)^2)
  spec_spax <- spectra_mse(toy$spectra, fit_spax$spectra)
  spec_ref <- spectra_mse(toy$spectra, fit_ref$spectra)

  expect_lte(recon_spax, recon_ref * 1.25 + 1e-6)
  expect_lte(spec_spax, spec_ref * 2 + 1e-6)
})
