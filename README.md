# SpaxNMF

`SpaxNMF` is an R package for IFU cube decomposition with regularized
non-negative matrix factorization. The current package centers on reshaping
cubes, fitting a baseline NMF model with smooth spectral components, and
inspecting the recovered spatial maps and spectra.

The model is

$$
X \approx A S
$$

where `X` is a spaxel-by-wavelength matrix, `A` contains spatial abundances,
and `S` contains component spectra.

Each spaxel spectrum is represented as a non-negative combination of shared
component spectra. The weights vary across the field, while the component
spectra are global to the fit.

## Installation

```r
# install.packages("remotes")
remotes::install_github("RafaelSdeSouza/SpaxNMF")
```

`torch` must be installed and configured in the local R environment.

## Basic usage

```r
library(SpaxNMF)

demo <- simulate_ifu_cube(nx = 16, ny = 16, n_wave = 100)

fit <- spax_nmf(
  demo$matrix,
  k = 3,
  lambda_smooth = 0.02,
  niter = 300,
  lr = 0.03
)

print(fit)
summary(fit)
```

## Accessors

The fitted object supports a compact interface inspired by `prcomp` and `NMF`.

```r
# abundance matrix (spaxels x components)
fit$spatial

# component spectra (components x wavelength)
fit$spectra

# NMF-style accessors
basis(fit)        # wavelength x components
coef(fit)         # components x spaxels

# fitted values and residuals
fitted(fit)
residuals(fit, x = demo$matrix)

# in-sample or new-data prediction
predict(fit)
predict(fit, newdata = demo$cube, type = "spatial")

# metadata carried by the fit
cube_metadata(fit)
```

## Visualization

```r
plot(fit, type = "spectra", wavelength = demo$wavelength)
plot(fit, type = "maps", nx = demo$nx, ny = demo$ny)
plot(fit, type = "loss")
plot_reconstruction(fit, demo$matrix, n = 4, wavelength = demo$wavelength)
```

## IFU workflow

```r
library(FITSio)

X <- readFITS("cube.fits")
Mat <- cube_to_matrix(X)
Mat[!is.finite(Mat)] <- 0
Mat <- pmax(Mat, 0)

fit <- spax_nmf(Mat, k = 5, lr = 0.01, niter = 5000)

maps <- predict(fit, type = "spatial")
cube_hat <- predict(
  fit,
  type = "cube",
  nx = dim(X$imDat)[1],
  ny = dim(X$imDat)[2]
)

# recover stored FITS-side metadata if present
meta <- cube_metadata(cube_hat)
```

When `cube_to_matrix()` receives a FITS-like list object, it now carries
non-image entries such as headers and other metadata through the matrix, the
fitted object, and reconstructed cubes.

## Roadmap

The package is being cleaned up around one main use case: IFU demonstrations
that compare what we learn from PCA, vanilla NMF, and a future spatially aware
NMF model on public MaNGA cubes. The synthetic example stays in the package so
the website and tests remain fast and reproducible, while the real-data demos
will move toward a small set of curated MaNGA examples.

## Documentation

Package articles are provided in the `vignettes/` directory.

Website: [https://rafaelsdesouza.github.io/SpaxNMF/](https://rafaelsdesouza.github.io/SpaxNMF/)
