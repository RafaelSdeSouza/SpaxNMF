source("R/spectral_unmix.R")

dir.create("site/images", recursive = TRUE, showWarnings = FALSE)

library(FITSio)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

set_all_seeds <- function(seed) {
  set.seed(seed)
  torch::torch_manual_seed(seed)
}

read_env_integer <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(as.integer(default))
  }

  parsed <- suppressWarnings(as.integer(value))
  if (!is.finite(parsed) || is.na(parsed) || parsed < 1L) {
    stop(sprintf("Environment variable %s must be a positive integer.", name), call. = FALSE)
  }

  parsed
}

read_env_numeric <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(as.numeric(default))
  }

  parsed <- suppressWarnings(as.numeric(value))
  if (!is.finite(parsed) || is.na(parsed)) {
    stop(sprintf("Environment variable %s must be numeric.", name), call. = FALSE)
  }

  parsed
}

read_env_path <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }

  value
}

read_env_string <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }

  value
}

parse_fits_header <- function(header_lines) {
  kv <- FITSio::parseHdr(header_lines)
  setNames(kv[seq(2, length(kv), by = 2)], kv[seq(1, length(kv), by = 2)])
}

odd_window <- function(n, max_k = 51L) {
  k <- min(as.integer(max_k), as.integer(n))
  if (k %% 2L == 0L) {
    k <- k - 1L
  }
  if (k < 3L) {
    return(3L)
  }

  k
}

read_env_seeds <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(as.integer(default))
  }

  pieces <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  parsed <- suppressWarnings(as.integer(pieces))
  if (!length(parsed) || any(!is.finite(parsed)) || any(is.na(parsed))) {
    stop(sprintf("Environment variable %s must be a comma-separated integer list.", name), call. = FALSE)
  }

  unique(parsed)
}

read_env_indices <- function(name, default = NULL) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }

  pieces <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  parsed <- suppressWarnings(as.integer(pieces))
  if (!length(parsed) || any(!is.finite(parsed)) || any(is.na(parsed))) {
    stop(sprintf("Environment variable %s must be a comma-separated integer list.", name), call. = FALSE)
  }

  unique(parsed)
}

fit_best_model <- function(x, seeds, ...) {
  best_fit <- NULL
  best_mse <- Inf

  for (seed in seeds) {
    set_all_seeds(seed)
    fit <- spectral_unmix(x, ...)
    mse <- mean((x - fit$reconstruction)^2)
    if (mse < best_mse) {
      best_fit <- fit
      best_mse <- mse
    }
  }

  best_fit
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
  fit$spectra <- fit$spectra[perm, , drop = FALSE]
  fit$spatial <- fit$spatial[, perm, drop = FALSE]
  fit$abundance <- fit$spatial
  fit$basis <- t(fit$spectra)
  fit$coef <- t(fit$spatial)
  fit$reconstruction <- fit$spatial %*% fit$spectra
  fit$fitted <- fit$reconstruction
  fit
}

display_positive <- function(x) {
  s <- max(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) {
    return(rep(0, length(x)))
  }

  x / s
}

maps_long <- function(maps, model, component_names) {
  bind_rows(lapply(seq_along(maps), function(i) {
    map <- maps[[i]]
    expand.grid(ix = seq_len(nrow(map)), iy = seq_len(ncol(map))) |>
      mutate(
        component = component_names[i],
        value = as.vector(map)
      )
  })) |>
    group_by(component) |>
    mutate(value = display_positive(value)) |>
    ungroup() |>
    mutate(
      model = model,
      component = factor(component, levels = component_names)
    )
}

residual_long <- function(map, model) {
  expand.grid(ix = seq_len(nrow(map)), iy = seq_len(ncol(map))) |>
    mutate(
      component = factor("Residual", levels = "Residual"),
      model = model,
      value = as.vector(map)
    )
}

residual_map <- function(x, reconstruction, nx, ny) {
  matrix(rowMeans((x - reconstruction)^2), nrow = nx, ncol = ny)
}

pick_default_cube <- function() {
  candidates <- c(
    Sys.getenv("SPAXNMF_REAL_CUBE", unset = ""),
    "/Users/rd23aag/Documents/GitHub/sagui_capivara_MaNGA/manga-7443-12703-LOGCUBE.fits",
    "/Users/rd23aag/Documents/GitHub/path_signatures/manga-7443-12703-LOGCUBE.fits",
    "/Users/rd23aag/Documents/GitHub/Test_S2D/manga-7443-12703-LOGCUBE.fits",
    "/Users/rd23aag/Documents/GitHub/HUB_2026/Niro/NGC1512_cube_mockGCs.fits"
  )
  candidates <- candidates[nzchar(candidates)]
  existing <- candidates[file.exists(candidates)]
  if (!length(existing)) {
    stop(
      "No local MaNGA LOGCUBE file was found. Set SPAXNMF_REAL_CUBE to a local path.",
      call. = FALSE
    )
  }

  existing[1]
}

read_ifu_cube <- function(path,
                          wave_min = 3600,
                          wave_max = 9000,
                          wave_step = 8,
                          min_valid_frac = 0.5) {
  wave_hdu <- try(readFITS(path, hdu = 6), silent = TRUE)

  if (!inherits(wave_hdu, "try-error")) {
    flux <- readFITS(path, hdu = 1)$imDat
    ivar <- readFITS(path, hdu = 2)$imDat
    mask <- readFITS(path, hdu = 3)$imDat
    wave <- as.numeric(wave_hdu$imDat)
    cube_kind <- "MaNGA LOGCUBE"

    keep_wave <- which(wave >= wave_min & wave <= wave_max)
    keep_wave <- keep_wave[seq(1L, length(keep_wave), by = as.integer(wave_step))]

    valid_voxel <- is.finite(flux[, , keep_wave, drop = FALSE]) &
      is.finite(ivar[, , keep_wave, drop = FALSE]) &
      (ivar[, , keep_wave, drop = FALSE] > 0) &
      (mask[, , keep_wave, drop = FALSE] == 0)
  } else {
    data_hdu <- readFITS(path, hdu = 1)
    stat_hdu <- try(readFITS(path, hdu = 2), silent = TRUE)
    flux <- data_hdu$imDat
    stat <- if (inherits(stat_hdu, "try-error")) NULL else stat_hdu$imDat
    hdr <- parse_fits_header(data_hdu$header)

    crval3 <- as.numeric(hdr[["CRVAL3"]])
    crpix3 <- as.numeric(hdr[["CRPIX3"]])
    cd3_3 <- suppressWarnings(as.numeric(hdr[["CD3_3"]]))
    if (!is.finite(cd3_3)) {
      cd3_3 <- suppressWarnings(as.numeric(hdr[["CDELT3"]]))
    }
    if (!is.finite(crval3) || !is.finite(crpix3) || !is.finite(cd3_3)) {
      stop("Could not reconstruct the wavelength axis from the FITS header.", call. = FALSE)
    }

    wave <- crval3 + (seq_len(dim(flux)[3]) - crpix3) * cd3_3
    cube_kind <- "Generic IFU cube"

    keep_wave <- which(wave >= wave_min & wave <= wave_max)
    keep_wave <- keep_wave[seq(1L, length(keep_wave), by = as.integer(wave_step))]

    valid_voxel <- is.finite(flux[, , keep_wave, drop = FALSE])
    if (!is.null(stat)) {
      valid_voxel <- valid_voxel &
        is.finite(stat[, , keep_wave, drop = FALSE]) &
        (stat[, , keep_wave, drop = FALSE] > 0)
    }
  }

  valid2d <- apply(valid_voxel, c(1, 2), mean) > min_valid_frac
  idx <- which(valid2d, arr.ind = TRUE)
  xr <- range(idx[, 1])
  yr <- range(idx[, 2])

  flux_crop <- flux[xr[1]:xr[2], yr[1]:yr[2], keep_wave, drop = FALSE]
  valid_crop <- valid_voxel[xr[1]:xr[2], yr[1]:yr[2], , drop = FALSE]

  flux_crop[!valid_crop] <- 0
  flux_crop[flux_crop < 0] <- 0

  cube_matrix <- cube_to_matrix(flux_crop)
  scale <- stats::quantile(cube_matrix[cube_matrix > 0], probs = 0.995, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) {
    scale <- 1
  }
  cube_matrix <- cube_matrix / scale

  list(
    path = path,
    cube_kind = cube_kind,
    cube = flux_crop / scale,
    matrix = cube_matrix,
    wavelength = wave[keep_wave],
    valid_spaxels = matrix(apply(valid_crop, c(1, 2), mean) > min_valid_frac, nrow = dim(flux_crop)[1]),
    flux_map = apply(flux_crop, c(1, 2), sum),
    nx = dim(flux_crop)[1],
    ny = dim(flux_crop)[2]
  )
}

fit_pca_reconstruction <- function(x, k) {
  pca <- stats::prcomp(x, center = TRUE, scale. = FALSE, rank. = k)
  scores <- pca$x[, seq_len(k), drop = FALSE]
  loadings <- pca$rotation[, seq_len(k), drop = FALSE]
  recon <- scores %*% t(loadings)
  recon <- sweep(recon, 2, pca$center, `+`)
  list(fit = pca, reconstruction = recon)
}

estimate_outer_sky <- function(x,
                               wavelength,
                               flux_map,
                               valid_spaxels,
                               nx,
                               ny,
                               radius_quantile = 0.8,
                               flux_quantile = 0.25,
                               line_quantile = 0.9,
                               min_spaxels = 200L) {
  coords <- expand.grid(ix = seq_len(nx), iy = seq_len(ny))
  valid <- as.vector(valid_spaxels)
  flux <- as.vector(flux_map)
  flux[!is.finite(flux)] <- 0

  weights <- pmax(flux, 0)
  if (!(sum(weights[valid]) > 0)) {
    weights <- rep(1, length(weights))
  }

  cx <- stats::weighted.mean(coords$ix[valid], w = weights[valid])
  cy <- stats::weighted.mean(coords$iy[valid], w = weights[valid])
  radius <- sqrt((coords$ix - cx)^2 + (coords$iy - cy)^2)

  outer_cut <- stats::quantile(radius[valid], probs = radius_quantile, na.rm = TRUE)
  outer <- valid & radius >= outer_cut
  if (!any(outer)) {
    outer <- valid
  }

  flux_cut <- stats::quantile(flux[outer], probs = flux_quantile, na.rm = TRUE)
  selected <- which(outer & flux <= flux_cut)

  if (length(selected) < min_spaxels) {
    fallback <- which(valid)
    fallback <- fallback[order(-radius[fallback], flux[fallback])]
    selected <- unique(c(selected, fallback[seq_len(min(length(fallback), min_spaxels))]))
  }

  template <- apply(x[selected, , drop = FALSE], 2, stats::median, na.rm = TRUE)
  smooth_k <- odd_window(length(template), max_k = 61L)
  continuum <- stats::runmed(template, k = smooth_k, endrule = "median")
  line_template <- pmax(template - continuum, 0)

  line_power <- sum(line_template^2, na.rm = TRUE)
  full_power <- sum(template^2, na.rm = TRUE)

  if (is.finite(line_power) && is.finite(full_power) && line_power > 1e-4 * full_power) {
    cutoff <- stats::quantile(line_template[line_template > 0], probs = line_quantile, na.rm = TRUE)
    line_mask <- line_template >= cutoff & line_template > 0
    fit_template <- line_template
    fit_mode <- "line_guided"
  } else {
    cutoff <- stats::quantile(template, probs = 0.8, na.rm = TRUE)
    line_mask <- template >= cutoff
    fit_template <- template
    fit_mode <- "full_template"
  }

  if (sum(line_mask) < 8L) {
    top_idx <- order(fit_template, decreasing = TRUE)[seq_len(min(12L, length(fit_template)))]
    line_mask <- rep(FALSE, length(fit_template))
    line_mask[top_idx] <- TRUE
  }

  fit_slice <- fit_template[line_mask]
  fit_norm <- sum(fit_slice^2)
  amplitudes <- if (fit_norm > 0) {
    pmax(0, as.vector(x[, line_mask, drop = FALSE] %*% fit_slice) / fit_norm)
  } else {
    rep(0, nrow(x))
  }

  sky_model <- amplitudes %o% template
  clean_matrix <- x - sky_model
  clean_matrix[clean_matrix < 0] <- 0

  selected_mask <- matrix(FALSE, nrow = nx, ncol = ny)
  selected_mask[selected] <- TRUE
  fraction <- rowSums(sky_model) / pmax(rowSums(x), 1e-8)
  fraction[!is.finite(fraction)] <- 0

  list(
    template = template,
    continuum = continuum,
    line_template = line_template,
    line_mask = line_mask,
    fit_mode = fit_mode,
    selected_idx = selected,
    selected_mask = selected_mask,
    amplitudes = amplitudes,
    amplitude_map = matrix(amplitudes, nrow = nx, ncol = ny),
    fraction_map = matrix(fraction, nrow = nx, ncol = ny),
    clean_matrix = clean_matrix,
    sky_model = sky_model,
    centroid = c(cx, cy),
    example_spaxel = if (length(selected)) selected[which.max(amplitudes[selected])] else integer()
  )
}

pick_example_spaxels <- function(x, nx, ny, valid_spaxels, n = 3, preferred = integer()) {
  totals <- rowSums(x)
  coords <- arrayInd(seq_len(length(totals)), .dim = c(nx, ny))
  valid_rows <- which(as.vector(valid_spaxels))
  candidates <- valid_rows[order(totals[valid_rows], decreasing = TRUE)][seq_len(min(50, length(valid_rows)))]

  chosen <- unique(preferred[preferred %in% valid_rows])
  candidates <- setdiff(candidates, chosen)
  target_n <- min(n, length(unique(c(chosen, candidates))))

  if (!length(chosen) && length(candidates)) {
    chosen <- c(chosen, candidates[1])
  }

  while (length(chosen) < target_n) {
    remaining <- setdiff(candidates, chosen)
    if (!length(remaining)) {
      break
    }
    best <- remaining[which.max(vapply(remaining, function(idx) {
      coord <- coords[idx, ]
      min(vapply(chosen, function(sel) {
        sum((coord - coords[sel, ])^2)
      }, numeric(1)))
    }, numeric(1)))]
    chosen <- c(chosen, best)
  }

  chosen[seq_len(min(n, length(chosen)))]
}

wave_step <- read_env_integer("SPAXNMF_REAL_WAVE_STEP", 8L)
wave_min <- read_env_numeric("SPAXNMF_REAL_WAVE_MIN", 3600)
wave_max <- read_env_numeric("SPAXNMF_REAL_WAVE_MAX", 9000)
min_valid_frac <- read_env_numeric("SPAXNMF_REAL_MIN_VALID_FRAC", 0.5)
sky_mode <- tolower(read_env_string("SPAXNMF_REAL_SKY_MODE", "none"))
sky_radius_quantile <- read_env_numeric("SPAXNMF_REAL_SKY_RADIUS_Q", 0.8)
sky_flux_quantile <- read_env_numeric("SPAXNMF_REAL_SKY_FLUX_Q", 0.25)
sky_line_quantile <- read_env_numeric("SPAXNMF_REAL_SKY_LINE_Q", 0.9)
sky_min_spaxels <- read_env_integer("SPAXNMF_REAL_SKY_MIN_SPAXELS", 200L)

real_cube <- read_ifu_cube(
  pick_default_cube(),
  wave_min = wave_min,
  wave_max = wave_max,
  wave_step = wave_step,
  min_valid_frac = min_valid_frac
)

if (sky_mode == "outer") {
  real_cube$matrix_raw <- real_cube$matrix
  real_cube$sky <- estimate_outer_sky(
    real_cube$matrix_raw,
    wavelength = real_cube$wavelength,
    flux_map = real_cube$flux_map,
    valid_spaxels = real_cube$valid_spaxels,
    nx = real_cube$nx,
    ny = real_cube$ny,
    radius_quantile = sky_radius_quantile,
    flux_quantile = sky_flux_quantile,
    line_quantile = sky_line_quantile,
    min_spaxels = sky_min_spaxels
  )
  real_cube$matrix <- real_cube$sky$clean_matrix
  real_cube$cube <- matrix_to_cube(real_cube$matrix, nx = real_cube$nx, ny = real_cube$ny)
  real_cube$flux_map <- matrix(rowSums(real_cube$matrix), nrow = real_cube$nx, ncol = real_cube$ny)
  sky_note <- sprintf(
    "outer-spaxel sky subtraction using %d spaxels (%s)",
    length(real_cube$sky$selected_idx),
    real_cube$sky$fit_mode
  )
} else if (sky_mode %in% c("none", "")) {
  real_cube$sky <- NULL
  sky_note <- "no explicit sky subtraction"
} else {
  stop("Unsupported SPAXNMF_REAL_SKY_MODE. Use 'none' or 'outer'.", call. = FALSE)
}

k <- read_env_integer("SPAXNMF_REAL_K", 6L)
fit_seeds <- read_env_seeds("SPAXNMF_REAL_SEEDS", c(11L, 23L, 37L, 41L, 53L))
fit_niter <- read_env_integer("SPAXNMF_REAL_NITER", 1000L)
fit_lr <- read_env_numeric("SPAXNMF_REAL_LR", 0.02)
lambda_smooth <- read_env_numeric("SPAXNMF_REAL_LAMBDA_SMOOTH", 0.001)
lambda_spatial <- read_env_numeric("SPAXNMF_REAL_LAMBDA_SPATIAL", 0.002)
lambda_sparse <- read_env_numeric("SPAXNMF_REAL_LAMBDA_SPARSE", 0)
spatial_sigma <- read_env_numeric("SPAXNMF_REAL_SPATIAL_SIGMA", 1)
smooth_components <- read_env_indices("SPAXNMF_REAL_SMOOTH_COMPONENTS", NULL)
sparse_components <- read_env_indices("SPAXNMF_REAL_SPARSE_COMPONENTS", NULL)
output_path <- read_env_path("SPAXNMF_REAL_OUTPUT", "site/images/manga-real-diagnostic.png")
component_names <- paste0("Component ", seq_len(k))

prior_note <- if (!is.null(sparse_components) && length(sparse_components)) {
  smooth_label <- if (is.null(smooth_components) || !length(smooth_components)) {
    "all others"
  } else {
    paste(smooth_components, collapse = ",")
  }
  sprintf(
    "mixed priors: smooth {%s}, sparse {%s}, lambda_sparse = %.3f",
    smooth_label,
    paste(sparse_components, collapse = ","),
    lambda_sparse
  )
} else {
  "smooth-only spatial prior"
}

fit_pca <- fit_pca_reconstruction(real_cube$matrix, k = k)
fit_vanilla <- fit_best_model(
  real_cube$matrix,
  seeds = fit_seeds,
  k = k,
  lambda_smooth = lambda_smooth,
  lambda_spatial = 0,
  lambda_sparse = 0,
  niter = fit_niter,
  lr = fit_lr
)
fit_spatial <- fit_best_model(
  real_cube$matrix,
  seeds = fit_seeds,
  k = k,
  lambda_smooth = lambda_smooth,
  lambda_spatial = lambda_spatial,
  lambda_sparse = lambda_sparse,
  spatial_sigma = spatial_sigma,
  smooth_components = smooth_components,
  sparse_components = sparse_components,
  nx = real_cube$nx,
  ny = real_cube$ny,
  niter = fit_niter,
  lr = fit_lr
)

dist_spatial <- outer(
  seq_len(nrow(fit_vanilla$spectra)),
  seq_len(nrow(fit_spatial$spectra)),
  Vectorize(function(i, j) {
    mean((display_positive(fit_vanilla$spectra[i, ]) - display_positive(fit_spatial$spectra[j, ]))^2)
  })
)
fit_spatial <- reorder_fit(fit_spatial, best_assignment(dist_spatial))

spectra_df <- bind_rows(
  as.data.frame(t(fit_vanilla$spectra)) |>
    setNames(component_names) |>
    mutate(wavelength = real_cube$wavelength) |>
    pivot_longer(cols = all_of(component_names), names_to = "component", values_to = "value") |>
    mutate(model = "Vanilla NMF", value = ave(value, component, FUN = display_positive)),
  as.data.frame(t(fit_spatial$spectra)) |>
    setNames(component_names) |>
    mutate(wavelength = real_cube$wavelength) |>
    pivot_longer(cols = all_of(component_names), names_to = "component", values_to = "value") |>
    mutate(model = "Spatial NMF", value = ave(value, component, FUN = display_positive))
) |>
  mutate(
    model = factor(model, levels = c("Vanilla NMF", "Spatial NMF")),
    component = factor(component, levels = component_names)
  )

maps_df <- bind_rows(
  maps_long(lapply(seq_len(ncol(fit_vanilla$spatial)), function(i) {
    matrix(fit_vanilla$spatial[, i], nrow = real_cube$nx, ncol = real_cube$ny)
  }), "Vanilla NMF", component_names),
  maps_long(lapply(seq_len(ncol(fit_spatial$spatial)), function(i) {
    matrix(fit_spatial$spatial[, i], nrow = real_cube$nx, ncol = real_cube$ny)
  }), "Spatial NMF", component_names)
) |>
  mutate(model = factor(model, levels = c("Vanilla NMF", "Spatial NMF")))

residual_df <- bind_rows(
  residual_long(residual_map(real_cube$matrix, fit_pca$reconstruction, real_cube$nx, real_cube$ny), "PCA"),
  residual_long(residual_map(real_cube$matrix, fit_vanilla$reconstruction, real_cube$nx, real_cube$ny), "Vanilla NMF"),
  residual_long(residual_map(real_cube$matrix, fit_spatial$reconstruction, real_cube$nx, real_cube$ny), "Spatial NMF")
) |>
  mutate(
    model = factor(model, levels = c("PCA", "Vanilla NMF", "Spatial NMF"))
  )

example_spaxels <- pick_example_spaxels(
  real_cube$matrix,
  nx = real_cube$nx,
  ny = real_cube$ny,
  valid_spaxels = real_cube$valid_spaxels,
  n = 3,
  preferred = if (is.null(real_cube$sky)) integer() else real_cube$sky$example_spaxel
)

example_df <- bind_rows(lapply(seq_along(example_spaxels), function(ii) {
  idx <- example_spaxels[ii]
  coord <- arrayInd(idx, .dim = c(real_cube$nx, real_cube$ny))
  tibble(
    wavelength = real_cube$wavelength,
    observed = real_cube$matrix[idx, ],
    pca = fit_pca$reconstruction[idx, ],
    vanilla = fit_vanilla$reconstruction[idx, ],
    spatial = fit_spatial$reconstruction[idx, ],
    spaxel = sprintf("Spaxel (%d, %d)", coord[1], coord[2])
  ) |>
    pivot_longer(
      cols = c("observed", "pca", "vanilla", "spatial"),
      names_to = "model",
      values_to = "value"
    )
}))

pca_mse <- mean((real_cube$matrix - fit_pca$reconstruction)^2)
vanilla_mse <- mean((real_cube$matrix - fit_vanilla$reconstruction)^2)
spatial_mse <- mean((real_cube$matrix - fit_spatial$reconstruction)^2)

base_panel_theme <- theme(
  plot.title = element_text(face = "bold", size = 16),
  axis.title = element_text(face = "bold"),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(color = "grey87", linewidth = 0.35),
  panel.border = element_rect(color = "grey84", fill = NA, linewidth = 0.7),
  strip.text = element_text(face = "bold"),
  plot.margin = margin(6, 6, 6, 6)
)

model_pal <- c(
  "observed" = "#111111",
  "pca" = "#2563EB",
  "vanilla" = "#F59E0B",
  "spatial" = "#0F766E",
  "Vanilla NMF" = "#F59E0B",
  "Spatial NMF" = "#0F766E"
)

map_pal <- c("#FFF7BC", "#FEC44F", "#FD8D3C", "#E6550D", "#A63603")

p_examples <- ggplot(example_df, aes(wavelength, value, color = model)) +
  geom_line(linewidth = 0.95) +
  facet_wrap(~ spaxel, nrow = 1) +
  scale_color_manual(
    values = model_pal[c("observed", "pca", "vanilla", "spatial")],
    labels = c("Observed", "PCA", "Vanilla NMF", "Spatial NMF")
  ) +
  scale_x_continuous(
    breaks = c(4000, 5000, 6000, 7000, 8000),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = if (is.null(real_cube$sky)) "Selected real-spaxel reconstructions" else "Selected sky-cleaned spaxel reconstructions",
    subtitle = sprintf(
      "k = %d | PCA MSE = %.6f | Vanilla = %.6f | Spatial = %.6f",
      k, pca_mse, vanilla_mse, spatial_mse
    ),
    x = expression(lambda ~ "[" * A * "]"),
    y = "Scaled flux",
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  base_panel_theme +
  theme(legend.position = "top")

p_sky <- NULL
if (!is.null(real_cube$sky)) {
  sky_df <- tibble(
    wavelength = real_cube$wavelength,
    template = display_positive(real_cube$sky$template),
    continuum = display_positive(real_cube$sky$continuum),
    sky_lines = display_positive(real_cube$sky$line_template)
  ) |>
    pivot_longer(
      cols = c("template", "continuum", "sky_lines"),
      names_to = "series",
      values_to = "value"
    )

  sky_map_df <- expand.grid(ix = seq_len(real_cube$nx), iy = seq_len(real_cube$ny)) |>
    mutate(
      value = as.vector(real_cube$sky$fraction_map),
      selected = as.vector(real_cube$sky$selected_mask)
    )

  p_sky_spec <- ggplot(sky_df, aes(wavelength, value, color = series)) +
    geom_line(linewidth = 0.95) +
    scale_color_manual(
      values = c(
        "template" = "#111111",
        "continuum" = "#F59E0B",
        "sky_lines" = "#0F766E"
      ),
      labels = c("Outer-spaxel median", "Smooth continuum", "Sky-line template")
    ) +
    scale_x_continuous(
      breaks = c(4000, 5000, 6000, 7000, 8000),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = "Estimated sky template from faint outer spaxels",
      subtitle = sprintf(
        "selected spaxels = %d | line channels = %d | mode = %s",
        length(real_cube$sky$selected_idx),
        sum(real_cube$sky$line_mask),
        real_cube$sky$fit_mode
      ),
      x = expression(lambda ~ "[" * A * "]"),
      y = "Norm. flux",
      color = NULL
    ) +
    theme_minimal(base_size = 13) +
    base_panel_theme +
    theme(legend.position = "top")

  p_sky_map <- ggplot(sky_map_df, aes(ix, iy, fill = value)) +
    geom_tile() +
    geom_point(
      data = subset(sky_map_df, selected),
      aes(ix, iy),
      inherit.aes = FALSE,
      color = "#111111",
      alpha = 0.18,
      size = 0.12
    ) +
    scale_fill_gradientn(colors = map_pal, guide = "none") +
    coord_equal() +
    labs(
      title = "Estimated sky weight map",
      subtitle = "Points mark the faint outer spaxels used to build the template",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "grey84", fill = NA, linewidth = 0.7),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", size = 16),
      plot.margin = margin(6, 6, 6, 6)
    )

  p_sky <- p_sky_spec | p_sky_map
}

p_spectra <- ggplot(spectra_df, aes(wavelength, value, color = model)) +
  geom_line(linewidth = 1.0) +
  facet_wrap(~ component, nrow = 1) +
  scale_color_manual(values = model_pal[c("Vanilla NMF", "Spatial NMF")]) +
  scale_x_continuous(
    breaks = c(4000, 5000, 6000, 7000, 8000),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, 1.05),
    labels = number_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Recovered basis spectra on the IFU cube",
    subtitle = sprintf(
      "Spatial NMF uses the weighted spatial prior; %s; component order is matched to vanilla NMF by spectral similarity",
      prior_note
    ),
    x = expression(lambda ~ "[" * A * "]"),
    y = "Norm. flux",
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  base_panel_theme +
  theme(legend.position = "top")

p_maps <- ggplot(maps_df, aes(ix, iy, fill = value)) +
  geom_tile() +
  facet_grid(model ~ component) +
  scale_fill_gradientn(colors = map_pal, limits = c(0, 1), guide = "none") +
  coord_equal() +
  labs(
    title = "Recovered component-weight maps",
    subtitle = sprintf(
      "%s | crop = %dx%d | wave bins = %d | restarts = %d | niter = %d",
      basename(real_cube$path),
      real_cube$nx,
      real_cube$ny,
      length(real_cube$wavelength),
      length(fit_seeds),
      fit_niter
    ),
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey84", fill = NA, linewidth = 0.7),
    strip.text = element_text(face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    plot.margin = margin(6, 6, 6, 6)
  )

p_resid <- ggplot(residual_df, aes(ix, iy, fill = value)) +
  geom_tile() +
  facet_grid(. ~ model) +
  scale_fill_gradientn(
    colors = map_pal,
    limits = c(0, stats::quantile(residual_df$value, probs = 0.99, na.rm = TRUE)),
    oob = scales::squish,
    guide = "none"
  ) +
  coord_equal() +
  labs(
    title = "Per-spaxel residual energy",
    subtitle = sprintf(
      "Secondary diagnostic at fixed k = %d; PCA is expected to minimize residual energy among rank-k linear models",
      k
    ),
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey84", fill = NA, linewidth = 0.7),
    strip.text = element_text(face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    plot.margin = margin(6, 6, 6, 6)
  )

fig <- if (is.null(p_sky)) {
  p_examples / p_spectra / p_maps / p_resid +
    plot_layout(heights = c(0.9, 0.9, 1.4, 0.8))
} else {
  p_examples / p_sky / p_spectra / p_maps / p_resid +
    plot_layout(heights = c(0.9, 0.8, 0.9, 1.4, 0.8))
}

fig <- fig +
  plot_annotation(
    title = "Real IFU cube diagnostic: PCA versus vanilla NMF versus spatial NMF",
    subtitle = sprintf(
      "%s with wavelength thinning, footprint cropping, and tuned NMF defaults (k = %d, %s, %s)",
      real_cube$cube_kind,
      k,
      sky_note,
      prior_note
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 20),
      plot.subtitle = element_text(size = 12, color = "grey25")
    )
  )

ggsave(
  filename = output_path,
  plot = fig,
  width = 13,
  height = 16,
  dpi = 320,
  bg = "white"
)
