\name{spectral_unmix}
\alias{spectral_unmix}
\alias{print.spectral_unmix}
\title{Regularized NMF for IFU Data}
\usage{
spectral_unmix(
  x,
  k = 3,
  lambda_smooth = 0.01,
  lr = 0.02,
  niter = 2000,
  center = FALSE,
  scale = FALSE,
  cuda = FALSE,
  verbose = FALSE,
  metadata = NULL
)

\method{print}{spectral_unmix}(x, ...)
}
\arguments{
\item{x}{Numeric matrix-like object with spectra in rows and wavelength
channels in columns.}

\item{k}{Number of components to estimate.}

\item{lambda_smooth}{Non-negative weight for the first-difference smoothness
penalty applied to the recovered spectra.}

\item{lr}{Learning rate passed to \code{torch::optim_adam()}.}

\item{niter}{Number of optimization iterations.}

\item{center}{Logical; if \code{TRUE}, center wavelength channels before
fitting.}

\item{scale}{Logical; if \code{TRUE}, scale wavelength channels before fitting.}

\item{cuda}{Logical; if \code{TRUE}, request CUDA execution.}

\item{verbose}{Logical; if \code{TRUE}, emit progress every 100 iterations.}

\item{metadata}{Optional metadata list to carry through the fitted object,
reconstructed cubes, and prediction outputs.}

\item{...}{Unused.}
}
\value{
A list with class \code{"spectral_unmix"} containing:
\itemize{
\item \code{spatial}: abundance matrix with one column per component.
\item \code{abundance}: alias of \code{spatial}.
\item \code{spectra}: component spectra with one row per component.
\item \code{basis}: component spectra arranged as wavelength-by-component.
\item \code{coef}: component weights arranged as component-by-spaxel.
\item \code{reconstruction}: reconstructed matrix \eqn{AS}.
\item \code{fitted}: alias of \code{reconstruction}.
\item \code{loss}: objective history over the optimization.
\item \code{metadata}: optional metadata carried from the input matrix or
provided explicitly.
\item \code{center}: centering values used during preprocessing, or
\code{FALSE}.
\item \code{scale}: scaling values used during preprocessing, or
\code{FALSE}.
\item \code{call}: matched function call.
}
}
\description{
Fits a non-negative matrix factorization model with an optional smoothness
penalty on the recovered spectra.
}
\details{
The factorization follows the linear mixing model
\deqn{X \approx A S}
where rows of \eqn{X} correspond to spectra from individual spatial elements,
\eqn{A} contains non-negative abundance weights, and \eqn{S} contains
non-negative spectral components.

This implementation is aimed at integral-field spectroscopy data after
reshaping the spatial dimensions into a two-dimensional matrix of spaxels by
wavelength. The objective function is
\deqn{\|X - AS\|^2 + \lambda \|D S\|^2}
where the second term penalizes rapid channel-to-channel variations in each
component spectrum.
}
\references{
Lee, D. D., and Seung, H. S. (1999). Learning the parts of objects by
non-negative matrix factorization. \emph{Nature}, 401, 788-791.
}
\examples{
cube <- simulate_ifu_cube(nx = 12, ny = 12, n_wave = 80)
\dontrun{
fit <- spectral_unmix(cube$matrix, k = 3, niter = 150, lr = 0.03)
print(fit)
component_map(fit, nx = cube$nx, ny = cube$ny, component = 1)
}
}
