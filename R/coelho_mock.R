coelho_stellar_subset_internal <- function() {
  subset_path <- system.file("extdata", "coelho-stellar-subset.rds", package = "SpaxNMF")
  if (!nzchar(subset_path)) {
    subset_path <- file.path("inst", "extdata", "coelho-stellar-subset.rds")
  }

  if (!file.exists(subset_path)) {
    stop(
      "Bundled Coelho subset was not found in inst/extdata. ",
      "Expected coelho-stellar-subset.rds.",
      call. = FALSE
    )
  }

  readRDS(subset_path)
}

coelho_class_prototypes_internal <- function(classes = c("hot_star", "a_f_star", "solar_like", "cool_dwarf"),
                                             n_members_per_class = 10,
                                             wavelength_min = 1500,
                                             wavelength_max = 9500,
                                             wave_step = 4,
                                             summary_method = c("median", "mean")) {
  lib <- coelho_stellar_subset_internal()
  summary_method <- match.arg(summary_method)

  keep_wave <- which(lib$wavelength >= wavelength_min & lib$wavelength <= wavelength_max)
  if (!length(keep_wave)) {
    stop("No wavelength channels satisfy the requested range.", call. = FALSE)
  }
  if (wave_step > 1L) {
    keep_wave <- keep_wave[seq(1L, length(keep_wave), by = as.integer(wave_step))]
  }

  wave <- lib$wavelength[keep_wave]
  spectra_lib <- lib$matrix[, keep_wave, drop = FALSE]
  metadata <- lib$metadata

  classes <- unique(as.character(classes))
  if (!length(classes)) {
    stop("'classes' must contain at least one class label.", call. = FALSE)
  }

  missing_classes <- setdiff(classes, unique(metadata$type))
  if (length(missing_classes)) {
    stop(
      "Unknown Coelho classes: ",
      paste(missing_classes, collapse = ", "),
      call. = FALSE
    )
  }

  prototypes <- matrix(0, nrow = length(classes), ncol = length(wave))
  rownames(prototypes) <- classes
  selected_members <- vector("list", length(classes))

  for (ii in seq_along(classes)) {
    class_name <- classes[ii]
    pool <- which(metadata$type == class_name)
    n_draw <- min(as.integer(n_members_per_class), length(pool))
    chosen <- sample(pool, size = n_draw, replace = FALSE)
    chosen_spectra <- spectra_lib[chosen, , drop = FALSE]

    prototype <- if (summary_method == "median") {
      apply(chosen_spectra, 2, stats::median)
    } else {
      colMeans(chosen_spectra)
    }

    scale <- max(prototype)
    if (!is.finite(scale) || scale <= 0) {
      stop("Coelho class prototype has non-positive scale.", call. = FALSE)
    }

    prototypes[ii, ] <- prototype / scale
    selected_members[[ii]] <- transform(
      metadata[chosen, , drop = FALSE],
      prototype_class = class_name,
      prototype_summary = summary_method
    )
  }

  list(
    wavelength = wave,
    spectra = prototypes,
    metadata = do.call(rbind, selected_members),
    classes = classes
  )
}

gaussian_blur_map_internal <- function(map, sigma = 0) {
  if (!is.numeric(sigma) || length(sigma) != 1L || is.na(sigma) || sigma < 0) {
    stop("'sigma' must be a non-negative number.", call. = FALSE)
  }
  if (sigma == 0) {
    return(map)
  }

  radius <- max(1L, as.integer(ceiling(3 * sigma)))
  offsets <- seq.int(-radius, radius)
  kernel <- stats::dnorm(offsets, mean = 0, sd = sigma)
  kernel <- kernel / sum(kernel)

  blur_along_axis <- function(mat, by_row = TRUE) {
    out <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
    if (by_row) {
      for (ii in seq_len(nrow(mat))) {
        for (jj in seq_len(ncol(mat))) {
          idx <- pmin(pmax(jj + offsets, 1L), ncol(mat))
          out[ii, jj] <- sum(mat[ii, idx] * kernel)
        }
      }
    } else {
      for (ii in seq_len(nrow(mat))) {
        for (jj in seq_len(ncol(mat))) {
          idx <- pmin(pmax(ii + offsets, 1L), 1L * nrow(mat))
          out[ii, jj] <- sum(mat[idx, jj] * kernel)
        }
      }
    }
    out
  }

  blur_along_axis(blur_along_axis(map, by_row = TRUE), by_row = FALSE)
}

simulate_coelho_recovery_cube <- function(nx = 18,
                                          ny = 18,
                                          classes = c("hot_star", "a_f_star", "solar_like", "cool_dwarf"),
                                          n_members_per_class = 10,
                                          wavelength_min = 1500,
                                          wavelength_max = 9500,
                                          wave_step = 4,
                                          summary_method = c("median", "mean"),
                                          psf_sigma = 1.1,
                                          noise = 0.003) {
  proto <- coelho_class_prototypes_internal(
    classes = classes,
    n_members_per_class = n_members_per_class,
    wavelength_min = wavelength_min,
    wavelength_max = wavelength_max,
    wave_step = wave_step,
    summary_method = summary_method
  )

  nx <- as.integer(nx)
  ny <- as.integer(ny)
  xgrid <- seq(-1, 1, length.out = nx)
  ygrid <- seq(-1, 1, length.out = ny)
  xy <- expand.grid(x = xgrid, y = ygrid)

  r_main <- sqrt((xy$x + 0.03)^2 + (xy$y + 0.02)^2)
  bulge <- exp(-(r_main^2) / 0.10)
  disk <- exp(-((xy$x / 1.05)^2 + (xy$y / 0.72)^2) / 0.70)
  ring <- exp(-((r_main - 0.56)^2) / 0.018)
  arm_a <- exp(-((xy$x + 0.42)^2 + (xy$y - 0.18)^2) / 0.055)
  arm_b <- exp(-((xy$x - 0.30)^2 + (xy$y + 0.26)^2) / 0.070)
  clump_a <- exp(-((xy$x + 0.12)^2 + (xy$y + 0.52)^2) / 0.025)
  clump_b <- exp(-((xy$x - 0.55)^2 + (xy$y - 0.10)^2) / 0.030)

  envelopes <- list(
    hot_star = 0.80 * ring + 0.70 * arm_a + 0.15 * clump_a,
    a_f_star = 0.50 * ring + 0.70 * arm_b + 0.10 * clump_b,
    solar_like = 0.95 * disk + 0.22 * bulge,
    cool_dwarf = 0.90 * bulge + 0.10 * disk,
    cool_giant = 0.55 * bulge + 0.30 * ring + 0.20 * clump_b
  )

  n_spaxels <- nx * ny
  abundance_matrix <- matrix(0, nrow = n_spaxels, ncol = length(proto$classes))
  colnames(abundance_matrix) <- proto$classes
  abundance_list <- vector("list", length(proto$classes))
  names(abundance_list) <- proto$classes

  for (ii in seq_along(proto$classes)) {
    class_name <- proto$classes[ii]
    map <- envelopes[[class_name]]
    if (is.null(map)) {
      stop("No spatial envelope is defined for class '", class_name, "'.", call. = FALSE)
    }

    map <- matrix(map, nrow = nx, ncol = ny)
    map <- gaussian_blur_map_internal(map, sigma = psf_sigma)
    map <- map / max(map)
    abundance_list[[ii]] <- map
    abundance_matrix[, ii] <- as.vector(map)
  }

  matrix_data <- abundance_matrix %*% proto$spectra
  matrix_data <- matrix_data + stats::rnorm(length(matrix_data), sd = noise)
  matrix_data[matrix_data < 0] <- 0

  list(
    matrix = matrix_data,
    wavelength = proto$wavelength,
    spectra = proto$spectra,
    abundances = abundance_list,
    metadata = proto$metadata,
    classes = proto$classes,
    nx = nx,
    ny = ny
  )
}

simulate_coelho_ifu_cube <- function(nx = 18,
                                     ny = 18,
                                     classes = c("hot_star", "a_f_star", "solar_like", "cool_dwarf"),
                                     n_members_per_class = 10,
                                     wavelength_min = 1500,
                                     wavelength_max = 9500,
                                     wave_step = 4,
                                     noise = 0.004) {
  lib <- coelho_stellar_subset_internal()

  keep_wave <- which(lib$wavelength >= wavelength_min & lib$wavelength <= wavelength_max)
  if (!length(keep_wave)) {
    stop("No wavelength channels satisfy the requested range.", call. = FALSE)
  }
  if (wave_step > 1L) {
    keep_wave <- keep_wave[seq(1L, length(keep_wave), by = as.integer(wave_step))]
  }

  wave <- lib$wavelength[keep_wave]
  spectra_lib <- lib$matrix[, keep_wave, drop = FALSE]
  metadata <- lib$metadata

  classes <- unique(as.character(classes))
  if (!length(classes)) {
    stop("'classes' must contain at least one class label.", call. = FALSE)
  }
  missing_classes <- setdiff(classes, unique(metadata$type))
  if (length(missing_classes)) {
    stop(
      "Unknown Coelho classes: ",
      paste(missing_classes, collapse = ", "),
      call. = FALSE
    )
  }

  xgrid <- seq(-1, 1, length.out = as.integer(nx))
  ygrid <- seq(-1, 1, length.out = as.integer(ny))
  xy <- expand.grid(x = xgrid, y = ygrid)

  r_main <- sqrt((xy$x + 0.05)^2 + (xy$y - 0.02)^2)
  bulge <- exp(-(r_main^2) / 0.12)
  disk <- exp(-((xy$x / 1.05)^2 + (xy$y / 0.75)^2) / 0.65)
  ring <- exp(-((r_main - 0.58)^2) / 0.020)
  arm_a <- exp(-((xy$x + 0.42)^2 + (xy$y - 0.18)^2) / 0.060)
  arm_b <- exp(-((xy$x - 0.28)^2 + (xy$y + 0.28)^2) / 0.080)

  envelopes <- list(
    hot_star = 0.55 * ring + 0.65 * arm_a + 0.20 * disk,
    a_f_star = 0.45 * ring + 0.45 * arm_b + 0.20 * disk,
    solar_like = 0.85 * disk + 0.30 * bulge,
    cool_dwarf = 0.85 * bulge + 0.20 * disk,
    cool_giant = 0.55 * bulge + 0.35 * ring + 0.15 * disk
  )

  n_spaxels <- as.integer(nx) * as.integer(ny)
  matrix_data <- matrix(0, nrow = n_spaxels, ncol = length(wave))
  class_spectra <- matrix(0, nrow = length(classes), ncol = length(wave))
  rownames(class_spectra) <- classes
  abundance_list <- vector("list", length(classes))
  names(abundance_list) <- classes
  member_info <- vector("list", length(classes))

  for (ii in seq_along(classes)) {
    class_name <- classes[ii]
    env <- envelopes[[class_name]]
    env <- env / max(env)

    pool <- which(metadata$type == class_name)
    n_draw <- min(as.integer(n_members_per_class), length(pool))
    chosen <- sample(pool, size = n_draw, replace = FALSE)

    class_map <- numeric(n_spaxels)
    class_weighted_spectrum <- numeric(length(wave))
    class_total_weight <- 0
    class_members <- vector("list", length(chosen))

    for (jj in seq_along(chosen)) {
      idx <- chosen[jj]
      center_id <- sample.int(n_spaxels, size = 1L, prob = env + 1e-8)
      cx <- xy$x[center_id] + stats::rnorm(1L, sd = 0.05)
      cy <- xy$y[center_id] + stats::rnorm(1L, sd = 0.05)
      sx <- stats::runif(1L, 0.10, 0.22)
      sy <- stats::runif(1L, 0.10, 0.24)
      local <- env * exp(-0.5 * (((xy$x - cx) / sx)^2 + ((xy$y - cy) / sy)^2))

      local_max <- max(local)
      if (!is.finite(local_max) || local_max <= 0) {
        next
      }
      local <- local / local_max
      amplitude <- stats::runif(1L, 0.35, 1.10)
      local <- amplitude * local

      spectrum <- spectra_lib[idx, ]
      total_weight <- sum(local)

      matrix_data <- matrix_data + tcrossprod(local, spectrum)
      class_map <- class_map + local
      class_weighted_spectrum <- class_weighted_spectrum + total_weight * spectrum
      class_total_weight <- class_total_weight + total_weight

      class_members[[jj]] <- data.frame(
        class = class_name,
        id = metadata$id[idx],
        filename = metadata$filename[idx],
        teff = metadata$teff[idx],
        logg = metadata$logg[idx],
        feh = metadata$feh[idx],
        amplitude = amplitude,
        total_weight = total_weight,
        stringsAsFactors = FALSE
      )
    }

    if (class_total_weight <= 0) {
      stop("Coelho semi-synthetic generator produced an empty class map.", call. = FALSE)
    }

    class_spectra[ii, ] <- class_weighted_spectrum / class_total_weight
    abundance_list[[ii]] <- matrix(class_map, nrow = as.integer(nx), ncol = as.integer(ny))
    member_info[[ii]] <- do.call(rbind, class_members)
  }

  matrix_data <- matrix_data + stats::rnorm(length(matrix_data), sd = noise)
  matrix_data[matrix_data < 0] <- 0

  list(
    matrix = matrix_data,
    wavelength = wave,
    spectra = class_spectra,
    abundances = abundance_list,
    metadata = do.call(rbind, member_info),
    classes = classes,
    nx = as.integer(nx),
    ny = as.integer(ny)
  )
}
