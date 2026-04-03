# SpaxNMF

`SpaxNMF` is an R package for IFU cube decomposition with regularized
non-negative matrix factorization. It provides tools for reshaping IFU cubes,
fitting spatially regularized NMF models, and inspecting recovered
component-weight maps and spectra.

The fitted model factorizes the data matrix as

$$
X \approx A S
$$

where `X` is a spaxel-by-wavelength matrix, `A` contains non-negative
component weights, and `S` contains non-negative component spectra.

The objective currently implemented in `spax_nmf()` is

$$
\mathcal{L}(A, S)=
\|X-AS\|_F^2
+\lambda_S \|DS\|_F^2
+\sum_c \lambda_{A,c}\sum_{i,j} w_{ij}(A_{ic}-A_{jc})^2
+\sum_c \lambda_{C,c}\|A_{:,c}\|_1
$$

with Gaussian spatial couplings

$$
w_{ij} = \exp\left(-\frac{d_{ij}^2}{2\sigma^2}\right).
$$

In words:

- `lambda_smooth` keeps component spectra smooth in wavelength.
- `lambda_spatial` couples nearby spaxels so selected component-weight maps are
  spatially coherent.
- `lambda_sparse` can be applied to chosen components to favor compact,
  GC-like or point-like structure.

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
# component-weight matrix (spaxels x components)
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

## Applications

SpaxNMF is intended for methodological studies and case-study analyses of IFU
data, including matched comparisons among PCA, vanilla NMF, and spatially
regularized NMF on curated survey cubes.

## Documentation

Package articles are provided in the `vignettes/` directory.

Website: [https://rafaelsdesouza.github.io/SpaxNMF/](https://rafaelsdesouza.github.io/SpaxNMF/)
