source("R/spax_nmf.R")

dir.create("site/images", recursive = TRUE, showWarnings = FALSE)
dir.create("site/data", recursive = TRUE, showWarnings = FALSE)

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

set_all_seeds <- function(seed) {
  set.seed(seed)
  torch::torch_manual_seed(seed)
}

read_env_path <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }

  value
}

display_positive <- function(x) {
  x <- as.numeric(x)
  s <- max(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) {
    return(rep(0, length(x)))
  }

  x / s
}

normalize_map <- function(map) {
  s <- max(map, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) {
    return(matrix(0, nrow = nrow(map), ncol = ncol(map)))
  }

  map / s
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
  fit$reconstruction <- fit$spatial %*% fit$spectra
  fit
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
    spectra = spectra,
    reconstruction = spatial %*% spectra
  )
}

fit_best_reference <- function(x, rank, seeds) {
  if (!requireNamespace("NMF", quietly = TRUE)) {
    stop("The 'NMF' package is required for the toy benchmark.", call. = FALSE)
  }

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

spectral_distance <- function(truth, fit) {
  outer(
    seq_len(nrow(truth)),
    seq_len(nrow(fit)),
    Vectorize(function(i, j) {
      mean((display_positive(truth[i, ]) - display_positive(fit[j, ]))^2)
    })
  )
}

maps_long <- function(maps, method, component_names) {
  bind_rows(lapply(seq_along(maps), function(i) {
    map <- normalize_map(maps[[i]])
    expand.grid(ix = seq_len(nrow(map)), iy = seq_len(ncol(map))) |>
      mutate(
        component = component_names[i],
        value = as.vector(map)
      )
  })) |>
    mutate(
      method = factor(method, levels = c("Truth", "SpaxNMF", "NMF package")),
      component = factor(component, levels = component_names)
    )
}

spectra_long <- function(spectra, wave, method, component_names) {
  as.data.frame(t(spectra)) |>
    setNames(component_names) |>
    mutate(wavelength = wave) |>
    pivot_longer(
      cols = all_of(component_names),
      names_to = "component",
      values_to = "value"
    ) |>
    group_by(component) |>
    mutate(value = display_positive(value)) |>
    ungroup() |>
    mutate(
      method = factor(method, levels = c("Truth", "SpaxNMF", "NMF package")),
      component = factor(component, levels = component_names)
    )
}

component_names <- paste("Component", 1:3)
toy_raw <- simulate_ifu_cube(nx = 12, ny = 12, n_wave = 80, noise = 0)
toy <- list(
  matrix = toy_raw$matrix,
  cube = toy_raw$cube,
  wavelength = toy_raw$wavelength,
  spectra = toy_raw$spectra,
  maps = toy_raw$abundances,
  nx = toy_raw$nx,
  ny = toy_raw$ny
)
k <- length(component_names)
fit_seeds <- c(11, 23, 37, 41, 53, 67)

fit_spax <- fit_best_spax(
  toy$matrix,
  seeds = fit_seeds,
  k = k,
  lambda_smooth = 0,
  lambda_spatial = 0,
  lambda_sparse = 0,
  solver = "adam",
  niter = 2000,
  lr = 0.03,
  tol = 0
)

fit_ref <- fit_best_reference(
  toy$matrix,
  rank = k,
  seeds = fit_seeds
)

fit_spax <- reorder_fit(fit_spax, best_assignment(spectral_distance(toy$spectra, fit_spax$spectra)))
fit_ref <- reorder_fit(fit_ref, best_assignment(spectral_distance(toy$spectra, fit_ref$spectra)))

truth_maps_mat <- do.call(cbind, lapply(toy$maps, as.vector))
spax_maps_mat <- fit_spax$spatial
ref_maps_mat <- fit_ref$spatial

spectra_mse <- function(truth, fit) {
  mean((t(apply(truth, 1, display_positive)) - t(apply(fit, 1, display_positive)))^2)
}

map_mse <- function(truth, fit, nx, ny) {
  truth_maps <- lapply(seq_len(ncol(truth)), function(i) {
    normalize_map(matrix(truth[, i], nrow = nx, ncol = ny))
  })
  fit_maps <- lapply(seq_len(ncol(fit)), function(i) {
    normalize_map(matrix(fit[, i], nrow = nx, ncol = ny))
  })

  mean(vapply(seq_len(length(truth_maps)), function(i) {
    mean((truth_maps[[i]] - fit_maps[[i]])^2)
  }, numeric(1)))
}

summary_df <- tibble(
  method = c("SpaxNMF", "NMF package"),
  solver = c(fit_spax$solver, "lee"),
  reconstruction_mse = c(
    mean((toy$matrix - fit_spax$reconstruction)^2),
    mean((toy$matrix - fit_ref$reconstruction)^2)
  ),
  spectra_mse = c(
    spectra_mse(toy$spectra, fit_spax$spectra),
    spectra_mse(toy$spectra, fit_ref$spectra)
  ),
  map_mse = c(
    map_mse(truth_maps_mat, spax_maps_mat, toy$nx, toy$ny),
    map_mse(truth_maps_mat, ref_maps_mat, toy$nx, toy$ny)
  )
)

spectra_df <- bind_rows(
  spectra_long(toy$spectra, toy$wavelength, "Truth", component_names),
  spectra_long(fit_spax$spectra, toy$wavelength, "SpaxNMF", component_names),
  spectra_long(fit_ref$spectra, toy$wavelength, "NMF package", component_names)
)

maps_df <- bind_rows(
  maps_long(toy$maps, "Truth", component_names),
  maps_long(lapply(seq_len(ncol(fit_spax$spatial)), function(i) {
    matrix(fit_spax$spatial[, i], nrow = toy$nx, ncol = toy$ny)
  }), "SpaxNMF", component_names),
  maps_long(lapply(seq_len(ncol(fit_ref$spatial)), function(i) {
    matrix(fit_ref$spatial[, i], nrow = toy$nx, ncol = toy$ny)
  }), "NMF package", component_names)
)

metrics_df <- summary_df |>
  mutate(
    solver = factor(solver, levels = c("lee", "adam")),
    reconstruction_mse = signif(reconstruction_mse, 4),
    spectra_mse = signif(spectra_mse, 4),
    map_mse = signif(map_mse, 4)
  ) |>
  pivot_longer(
    cols = c("reconstruction_mse", "spectra_mse", "map_mse"),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(
    metric = factor(
      metric,
      levels = c("reconstruction_mse", "spectra_mse", "map_mse"),
      labels = c("Reconstruction MSE", "Spectra MSE", "Map MSE")
    )
  )

component_pal <- c(
  "Component 1" = "#2563EB",
  "Component 2" = "#8B5CF6",
  "Component 3" = "#F59E0B"
)

map_pal <- c("#FFF7BC", "#FEC44F", "#FD8D3C", "#E6550D", "#A63603")

base_panel_theme <- theme(
  plot.title = element_text(face = "bold", size = 16),
  axis.title = element_text(face = "bold"),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(color = "grey87", linewidth = 0.35),
  panel.border = element_rect(color = "grey84", fill = NA, linewidth = 0.7),
  strip.text = element_text(face = "bold"),
  plot.margin = margin(6, 6, 6, 6)
)

p_spectra <- ggplot(spectra_df, aes(wavelength, value, color = component)) +
  geom_line(linewidth = 1.05) +
  facet_grid(method ~ .) +
  scale_color_manual(values = component_pal) +
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
    title = "Controlled low-rank spectra recovery",
    subtitle = "Bundled exact low-rank toy cube with all spatial and sparsity penalties turned off.",
    x = expression(lambda ~ "[" * A * "]"),
    y = "Norm. flux",
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  base_panel_theme +
  theme(legend.position = "top")

p_maps <- ggplot(maps_df, aes(ix, iy, fill = value)) +
  geom_tile() +
  facet_grid(method ~ component) +
  scale_fill_gradientn(colors = map_pal, limits = c(0, 1), guide = "none") +
  coord_equal() +
  labs(
    title = "Recovered component-weight maps",
    subtitle = "Columns are matched to the true components after fitting.",
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

p_metrics <- ggplot(metrics_df, aes(metric, method, fill = value)) +
  geom_tile(color = "white", linewidth = 0.9) +
  geom_text(aes(label = format(value, digits = 3, scientific = TRUE)), fontface = "bold", size = 4) +
  scale_fill_gradient(low = "#F3F4F6", high = "#1D4ED8", guide = "none") +
  labs(
    title = "Toy benchmark summary",
    subtitle = "Lower is better for all three metrics.",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey84", fill = NA, linewidth = 0.7),
    axis.text.x = element_text(angle = 18, hjust = 1),
    plot.title = element_text(face = "bold", size = 16),
    plot.margin = margin(6, 6, 6, 6)
  )

fig <- p_spectra / (p_maps | p_metrics) +
  plot_layout(heights = c(1.1, 1.2), widths = c(1.8, 1.0)) +
  plot_annotation(
    title = "Standalone NMF sanity check: SpaxNMF versus NMF package",
    subtitle = "Built-in exact low-rank toy cube with shared non-negative spectra and spatially varying weights",
    theme = theme(
      plot.title = element_text(face = "bold", size = 20),
      plot.subtitle = element_text(size = 12, color = "grey25")
    )
  )

output_path <- read_env_path("SPAXNMF_TOY_OUTPUT", "site/images/toy-nmf-benchmark.png")
summary_output <- read_env_path("SPAXNMF_TOY_SUMMARY_OUTPUT", "site/data/toy-nmf-benchmark-summary.csv")

ggsave(
  filename = output_path,
  plot = fig,
  width = 14,
  height = 11,
  dpi = 320,
  bg = "white"
)

write.csv(summary_df, summary_output, row.names = FALSE)
