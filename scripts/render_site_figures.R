source("R/spectral_unmix.R")

dir.create("site/images", recursive = TRUE, showWarnings = FALSE)

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggrepel)
library(scales)

nmf_multiplicative <- function(x, k, niter = 400, seed = 42) {
  set.seed(seed)
  n <- nrow(x)
  p <- ncol(x)
  w <- matrix(stats::runif(n * k), nrow = n, ncol = k)
  h <- matrix(stats::runif(k * p), nrow = k, ncol = p)
  eps <- 1e-8

  for (i in seq_len(niter)) {
    h <- h * ((t(w) %*% x) / (t(w) %*% w %*% h + eps))
    w <- w * ((x %*% t(h)) / (w %*% h %*% t(h) + eps))
  }

  list(w = w, h = h, reconstruction = w %*% h)
}

display_signed <- function(x) {
  x <- as.numeric(x)
  s <- max(abs(x), na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  x / s
}

display_positive <- function(x) {
  x <- as.numeric(x)
  s <- max(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  x / s
}

chi2_like <- function(x, y) {
  mean((x - y)^2, na.rm = TRUE)
}

best_assignment <- function(dist_mat) {
  comps <- colnames(dist_mat)
  refs <- rownames(dist_mat)
  n_ref <- nrow(dist_mat)

  perms <- expand.grid(rep(list(seq_len(ncol(dist_mat))), n_ref))
  perms <- as.matrix(perms)
  perms <- perms[apply(perms, 1, function(z) length(unique(z)) == n_ref), , drop = FALSE]

  scores <- apply(perms, 1, function(p) {
    sum(mapply(function(rr, jj) dist_mat[rr, comps[jj]], refs, p))
  })

  best <- perms[which.min(scores), ]

  data.frame(
    reference = refs,
    component = comps[best],
    distance = mapply(function(rr, cc) dist_mat[rr, cc], refs, comps[best]),
    stringsAsFactors = FALSE
  )
}

demo <- simulate_ifu_cube(nx = 16, ny = 16, n_wave = 100, noise = 0.005)

spatial <- do.call(
  cbind,
  lapply(demo$abundances, function(x) as.vector(x))
)

fit <- list(
  spatial = spatial,
  abundance = spatial,
  spectra = demo$spectra,
  basis = t(demo$spectra),
  coef = t(spatial),
  reconstruction = demo$matrix,
  fitted = demo$matrix,
  loss = c(seq(0.2, 0.02, length.out = 20), seq(0.019, 0.01, length.out = 20)),
  metadata = list(source = "simulate_ifu_cube"),
  center = FALSE,
  scale = FALSE,
  call = quote(simulate_ifu_cube())
)
class(fit) <- "spectral_unmix"

grDevices::png(
  filename = "site/images/simulated-spectra.png",
  width = 1400,
  height = 900,
  res = 180
)
plot(fit, type = "spectra", wavelength = demo$wavelength)
grDevices::dev.off()

grDevices::png(
  filename = "site/images/simulated-maps.png",
  width = 1400,
  height = 900,
  res = 180
)
plot(fit, type = "maps", nx = demo$nx, ny = demo$ny)
grDevices::dev.off()

grDevices::png(
  filename = "site/images/simulated-reconstruction.png",
  width = 1400,
  height = 900,
  res = 180
)
plot_reconstruction(fit, demo$matrix, n = 4, wavelength = demo$wavelength)
grDevices::dev.off()

wave <- demo$wavelength
k <- nrow(demo$spectra)

fit_pca <- stats::prcomp(
  demo$matrix,
  center = TRUE,
  scale. = FALSE,
  rank. = k
)

fit_nmf <- nmf_multiplicative(demo$matrix, k = k, niter = 500)

true_disp <- t(apply(demo$spectra, 1, display_positive))
rownames(true_disp) <- paste0("True", seq_len(k))

eig <- fit_pca$rotation[, seq_len(k), drop = FALSE]
flip <- sign(colSums(eig))
flip[flip == 0] <- 1
eig <- sweep(eig, 2, flip, `*`)
eig_disp <- apply(eig, 2, display_signed)
if (is.vector(eig_disp)) {
  eig_disp <- matrix(eig_disp, ncol = 1)
}
colnames(eig_disp) <- paste0("PC", seq_len(k))

nmf_disp <- t(apply(fit_nmf$h, 1, display_positive))
rownames(nmf_disp) <- paste0("NMF", seq_len(k))

dist_pca <- matrix(
  NA_real_,
  nrow = k,
  ncol = k,
  dimnames = list(rownames(true_disp), paste0("PC", seq_len(k)))
)

dist_nmf <- matrix(
  NA_real_,
  nrow = k,
  ncol = k,
  dimnames = list(rownames(true_disp), paste0("NMF", seq_len(k)))
)

for (j in seq_len(k)) {
  vpc <- eig_disp[, j]
  vnmf <- nmf_disp[j, ]

  for (rr in rownames(true_disp)) {
    ref <- true_disp[rr, ]
    dist_pca[rr, j] <- min(chi2_like(vpc, ref), chi2_like(-vpc, ref))
    dist_nmf[rr, j] <- chi2_like(vnmf, ref)
  }
}

assign_pca <- best_assignment(dist_pca)
assign_nmf <- best_assignment(dist_nmf)

ref_pal <- c(
  "True1" = "#2563EB",
  "True2" = "#F59E0B",
  "True3" = "#10B981"
)

pca_pal <- setNames(ref_pal[assign_pca$reference], assign_pca$component)
nmf_pal <- setNames(ref_pal[assign_nmf$reference], assign_nmf$component)

x_rng <- c(min(wave), max(wave))
x_breaks <- pretty(wave, n = 4)

base_panel_theme <- theme(
  plot.title = element_text(face = "bold", size = 17),
  axis.title = element_text(face = "bold"),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(color = "grey87", linewidth = 0.35),
  panel.border = element_rect(color = "grey84", fill = NA, linewidth = 0.7),
  plot.margin = margin(6, 6, 6, 6)
)

common_x_scale <- scale_x_continuous(
  limits = x_rng,
  breaks = x_breaks,
  expand = expansion(mult = c(0, 0))
)

true_df <- as.data.frame(t(true_disp))
colnames(true_df) <- rownames(true_disp)
true_df$wavelength <- wave
true_df <- true_df |>
  pivot_longer(
    cols = starts_with("True"),
    names_to = "reference",
    values_to = "value"
  )

p1 <- ggplot(true_df, aes(wavelength, value, color = reference)) +
  geom_line(linewidth = 1.15) +
  scale_color_manual(
    values = ref_pal,
    labels = c("Component 1", "Component 2", "Component 3")
  ) +
  common_x_scale +
  scale_y_continuous(
    limits = c(0, 1.05),
    labels = number_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "True synthetic spectra",
    x = expression(lambda ~ "[" * A * "]"),
    y = "Norm. flux",
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  base_panel_theme +
  theme(
    legend.position = c(0.60, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = NA, color = NA)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

pca_df <- as.data.frame(eig_disp)
colnames(pca_df) <- paste0("PC", seq_len(k))
pca_df$wavelength <- wave
pca_df <- pca_df |>
  pivot_longer(
    cols = starts_with("PC"),
    names_to = "component",
    values_to = "value"
  ) |>
  left_join(assign_pca, by = "component")

pca_annot <- pca_df |>
  group_by(component) |>
  filter(wavelength == max(wavelength)) |>
  slice(1) |>
  ungroup() |>
  distinct(component, .keep_all = TRUE) |>
  mutate(
    label = component,
    xlab = x_rng[2] - 0.04 * diff(x_rng)
  )

p2 <- ggplot(pca_df, aes(wavelength, value, color = component)) +
  geom_hline(yintercept = 0, color = "grey82", linewidth = 0.45) +
  geom_line(linewidth = 0.95) +
  geom_text_repel(
    data = pca_annot,
    aes(x = xlab, label = label),
    direction = "y",
    hjust = 1,
    nudge_x = 0,
    segment.color = scales::alpha("grey45", 0.5),
    segment.size = 0.35,
    size = 3.05,
    box.padding = 0.22,
    point.padding = 0.1,
    min.segment.length = 0,
    seed = 1,
    show.legend = FALSE
  ) +
  scale_color_manual(values = pca_pal, guide = "none") +
  common_x_scale +
  coord_cartesian(ylim = c(-1.05, 1.05), clip = "on") +
  labs(
    title = "PCA eigenspectra",
    x = expression(lambda ~ "[" * A * "]"),
    y = "Norm. flux"
  ) +
  theme_minimal(base_size = 14) +
  base_panel_theme

nmf_df <- as.data.frame(t(nmf_disp))
colnames(nmf_df) <- rownames(nmf_disp)
nmf_df$wavelength <- wave
nmf_df <- nmf_df |>
  pivot_longer(
    cols = starts_with("NMF"),
    names_to = "component",
    values_to = "value"
  ) |>
  left_join(assign_nmf, by = "component")

nmf_annot <- nmf_df |>
  group_by(component) |>
  filter(wavelength == max(wavelength)) |>
  slice(1) |>
  ungroup() |>
  distinct(component, .keep_all = TRUE) |>
  mutate(
    label = component,
    xlab = x_rng[2] - 0.04 * diff(x_rng)
  )

p3 <- ggplot(nmf_df, aes(wavelength, value, color = component)) +
  geom_line(linewidth = 0.95) +
  geom_text_repel(
    data = nmf_annot,
    aes(x = xlab, label = label),
    direction = "y",
    hjust = 1,
    nudge_x = 0,
    segment.color = scales::alpha("grey45", 0.5),
    segment.size = 0.35,
    size = 3.05,
    box.padding = 0.22,
    point.padding = 0.1,
    min.segment.length = 0,
    seed = 2,
    show.legend = FALSE
  ) +
  scale_color_manual(values = nmf_pal, guide = "none") +
  common_x_scale +
  coord_cartesian(ylim = c(0, 1.05), clip = "on") +
  labs(
    title = "NMF eigenspectra",
    x = expression(lambda ~ "[" * A * "]"),
    y = "Norm. flux"
  ) +
  theme_minimal(base_size = 14) +
  base_panel_theme

dist_df <- rbind(
  data.frame(
    reference = rep(rownames(dist_pca), times = ncol(dist_pca)),
    component = rep(colnames(dist_pca), each = nrow(dist_pca)),
    method = "PCA",
    distance = as.vector(dist_pca),
    stringsAsFactors = FALSE
  ),
  data.frame(
    reference = rep(rownames(dist_nmf), times = ncol(dist_nmf)),
    component = rep(colnames(dist_nmf), each = nrow(dist_nmf)),
    method = "NMF",
    distance = as.vector(dist_nmf),
    stringsAsFactors = FALSE
  )
)

dist_df <- dist_df |>
  mutate(
    method = factor(method, levels = c("PCA", "NMF")),
    reference = factor(reference, levels = rev(rownames(true_disp)))
  ) |>
  group_by(method) |>
  mutate(
    closeness = 1 - distance / max(distance, na.rm = TRUE),
    txt_col = ifelse(closeness > 0.55, "white", "#111111"),
    label = number(distance, accuracy = 0.01)
  ) |>
  ungroup()

reference_ann <- data.frame(
  method = factor("PCA", levels = c("PCA", "NMF")),
  reference = factor(rev(rownames(true_disp)), levels = rev(rownames(true_disp))),
  x = 0.42,
  label = c("Component 3", "Component 2", "Component 1"),
  col = unname(ref_pal[rev(rownames(true_disp))]),
  stringsAsFactors = FALSE
)

p4 <- ggplot(dist_df, aes(component, reference, fill = closeness)) +
  geom_tile(color = "white", linewidth = 0.9) +
  geom_text(aes(label = label, color = txt_col), size = 4.0, fontface = "bold") +
  geom_text(
    data = reference_ann,
    aes(x = x, y = reference, label = label, color = col),
    inherit.aes = FALSE,
    hjust = 1,
    fontface = "bold",
    size = 4.0
  ) +
  facet_grid(. ~ method, scales = "free_x", space = "free_x") +
  scale_fill_gradient(
    low = "#F3F4F6",
    high = "#1D4ED8",
    guide = "none"
  ) +
  scale_color_identity() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey84", fill = NA, linewidth = 0.7),
    strip.text = element_text(face = "bold", size = 12),
    plot.margin = margin(6, 6, 6, 52)
  )

fig <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Synthetic IFU components versus PCA and NMF eigenspectra",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 20)
    )
  )

ggsave(
  filename = "site/images/synthetic-basis-mosaic.png",
  plot = fig,
  width = 13,
  height = 10,
  dpi = 320,
  bg = "white"
)

# This preserves the published visual language while we prepare the equivalent
# MaNGA comparison figure for PCA, vanilla NMF, and spatial NMF.
