source("R/spectral_unmix.R")
source("R/coelho_mock.R")

dir.create("site/images", recursive = TRUE, showWarnings = FALSE)

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

set_all_seeds <- function(seed) {
  set.seed(seed)
  torch::torch_manual_seed(seed)
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

display_positive <- function(x) {
  s <- max(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) {
    return(rep(0, length(x)))
  }

  x / s
}

mean_map_roughness <- function(maps) {
  mean(vapply(maps, function(map) {
    scale <- max(map, na.rm = TRUE)
    if (is.finite(scale) && scale > 0) {
      map <- map / scale
    }
    dx <- if (nrow(map) > 1L) diff(map) else 0
    dy <- if (ncol(map) > 1L) t(diff(t(map))) else 0
    mean(dx^2) + mean(dy^2)
  }, numeric(1)))
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

pretty_class <- function(x) {
  tools::toTitleCase(gsub("_", " ", x, fixed = TRUE))
}

spectra_long <- function(spectra, wave, model, component_names) {
  as.data.frame(t(spectra)) |>
    setNames(component_names) |>
    mutate(wavelength = wave) |>
    pivot_longer(
      cols = all_of(component_names),
      names_to = "component",
      values_to = "value"
    ) |>
    mutate(
      model = model,
      value = ave(value, component, FUN = display_positive),
      component = factor(component, levels = component_names)
    )
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

residual_map <- function(x, reconstruction, nx, ny) {
  matrix(rowMeans((x - reconstruction)^2), nrow = nx, ncol = ny)
}

demo <- simulate_coelho_recovery_cube(
  nx = 16,
  ny = 16,
  n_members_per_class = 12,
  wave_step = 4,
  summary_method = "median",
  psf_sigma = 1.1,
  noise = 0.003
)

component_names <- pretty_class(demo$classes)
k <- nrow(demo$spectra)

fit_vanilla <- fit_best_model(
  demo$matrix,
  seeds = c(11, 23, 37),
  k = k,
  lambda_smooth = 0.001,
  lambda_spatial = 0,
  niter = 1200,
  lr = 0.03
)

fit_spatial <- fit_best_model(
  demo$matrix,
  seeds = c(11, 23, 37),
  k = k,
  lambda_smooth = 0.001,
  lambda_spatial = 0.002,
  spatial_sigma = 1,
  nx = demo$nx,
  ny = demo$ny,
  niter = 1200,
  lr = 0.03
)

dist_vanilla <- outer(
  seq_len(nrow(demo$spectra)),
  seq_len(nrow(fit_vanilla$spectra)),
  Vectorize(function(i, j) {
    mean((display_positive(demo$spectra[i, ]) - display_positive(fit_vanilla$spectra[j, ]))^2)
  })
)

dist_spatial <- outer(
  seq_len(nrow(demo$spectra)),
  seq_len(nrow(fit_spatial$spectra)),
  Vectorize(function(i, j) {
    mean((display_positive(demo$spectra[i, ]) - display_positive(fit_spatial$spectra[j, ]))^2)
  })
)

fit_vanilla <- reorder_fit(fit_vanilla, best_assignment(dist_vanilla))
fit_spatial <- reorder_fit(fit_spatial, best_assignment(dist_spatial))

spectra_df <- bind_rows(
  spectra_long(demo$spectra, demo$wavelength, "Truth", component_names),
  spectra_long(fit_vanilla$spectra, demo$wavelength, "Vanilla NMF", component_names),
  spectra_long(fit_spatial$spectra, demo$wavelength, "Spatial NMF", component_names)
) |>
  mutate(model = factor(model, levels = c("Truth", "Vanilla NMF", "Spatial NMF")))

maps_df <- bind_rows(
  maps_long(demo$abundances, "Truth", component_names),
  maps_long(lapply(seq_len(ncol(fit_vanilla$spatial)), function(i) {
    matrix(fit_vanilla$spatial[, i], nrow = demo$nx, ncol = demo$ny)
  }), "Vanilla NMF", component_names),
  maps_long(lapply(seq_len(ncol(fit_spatial$spatial)), function(i) {
    matrix(fit_spatial$spatial[, i], nrow = demo$nx, ncol = demo$ny)
  }), "Spatial NMF", component_names)
) |>
  mutate(model = factor(model, levels = c("Truth", "Vanilla NMF", "Spatial NMF")))

residual_df <- bind_rows(
  maps_long(
    list(residual_map(demo$matrix, fit_vanilla$reconstruction, demo$nx, demo$ny)),
    "Vanilla NMF",
    "Residual"
  ),
  maps_long(
    list(residual_map(demo$matrix, fit_spatial$reconstruction, demo$nx, demo$ny)),
    "Spatial NMF",
    "Residual"
  )
) |>
  mutate(
    component = factor("Residual", levels = "Residual"),
    model = factor(model, levels = c("Vanilla NMF", "Spatial NMF"))
  )

mse_vanilla <- mean((demo$matrix - fit_vanilla$reconstruction)^2)
mse_spatial <- mean((demo$matrix - fit_spatial$reconstruction)^2)
rough_true <- mean_map_roughness(demo$abundances)
rough_vanilla <- mean_map_roughness(lapply(seq_len(ncol(fit_vanilla$spatial)), function(i) {
  matrix(fit_vanilla$spatial[, i], nrow = demo$nx, ncol = demo$ny)
}))
rough_spatial <- mean_map_roughness(lapply(seq_len(ncol(fit_spatial$spatial)), function(i) {
  matrix(fit_spatial$spatial[, i], nrow = demo$nx, ncol = demo$ny)
}))

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
  "Truth" = "#111111",
  "Vanilla NMF" = "#F59E0B",
  "Spatial NMF" = "#0F766E"
)

map_pal <- c("#FFF7BC", "#FEC44F", "#FD8D3C", "#E6550D", "#A63603")

p_spectra <- ggplot(spectra_df, aes(wavelength, value, color = model)) +
  geom_line(linewidth = 1.05) +
  facet_wrap(~ component, nrow = 1) +
  scale_color_manual(values = model_pal) +
  scale_x_continuous(
    breaks = c(2000, 4000, 6000, 8000),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0, 1.05),
    labels = number_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Recovered spectra on the low-rank Coelho recovery cube",
    subtitle = "The cube is generated from exact class prototypes estimated from many Coelho spectra, then blurred spatially with a PSF-like kernel",
    x = expression(lambda ~ "[" * A * "]"),
    y = "Norm. flux",
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  base_panel_theme +
  theme(
    legend.position = "top"
  )

p_maps <- ggplot(maps_df, aes(ix, iy, fill = value)) +
  geom_tile() +
  facet_grid(model ~ component) +
  scale_fill_gradientn(colors = map_pal, limits = c(0, 1), guide = "none") +
  coord_equal() +
  labs(
    title = "Abundance maps",
    subtitle = paste0(
      "Truth roughness = ", number(rough_true, accuracy = 0.001),
      " | Vanilla = ", number(rough_vanilla, accuracy = 0.001),
      " | Spatial = ", number(rough_spatial, accuracy = 0.001)
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
  scale_fill_gradientn(colors = map_pal, limits = c(0, 1), guide = "none") +
  coord_equal() +
  labs(
    title = "Per-spaxel residual energy",
    subtitle = paste0(
      "Vanilla MSE = ", number(mse_vanilla, accuracy = 0.0001),
      " | Spatial MSE = ", number(mse_spatial, accuracy = 0.0001)
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

fig <- p_spectra / p_maps / p_resid +
  plot_layout(heights = c(1.0, 1.5, 0.8)) +
  plot_annotation(
    title = "Vanilla versus spatially regularized torch NMF",
    subtitle = "Low-rank Coelho recovery benchmark with exact class prototypes and PSF-smoothed abundance maps",
    theme = theme(
      plot.title = element_text(face = "bold", size = 20),
      plot.subtitle = element_text(size = 12, color = "grey25")
    )
  )

ggsave(
  filename = "site/images/spatial-vs-vanilla-diagnostic.png",
  plot = fig,
  width = 13,
  height = 14,
  dpi = 320,
  bg = "white"
)
