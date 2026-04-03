\name{plot.spax_nmf}
\alias{plot.spax_nmf}
\alias{plot_reconstruction}
\title{Visualization Helpers for Spectral Unmixing Results}
\usage{
\method{plot}{spax_nmf}(
  x,
  type = c("spectra", "maps", "loss"),
  nx = NULL,
  ny = NULL,
  components = NULL,
  wavelength = NULL,
  ...
)

plot_reconstruction(fit, x, pixels = NULL, n = 6, wavelength = NULL, ...)
}
\arguments{
\item{x}{For \code{plot.spax_nmf()}, a result from
\code{spax_nmf()}. For \code{plot_reconstruction()}, the numeric matrix
used as input to the fit.}

\item{type}{Plot type: \code{"spectra"}, \code{"maps"}, or \code{"loss"}.}

\item{nx}{Number of x-axis pixels when plotting maps.}

\item{ny}{Number of y-axis pixels when plotting maps.}

\item{components}{Optional vector of component indices to display.}

\item{wavelength}{Optional wavelength vector.}

\item{...}{Additional graphical arguments passed to base plotting functions.}

\item{fit}{Result from \code{spax_nmf()}.}

\item{pixels}{Optional vector of spaxel indices to plot.}

\item{n}{Number of random spaxels to plot when \code{pixels} is omitted.}
}
\value{
These functions invisibly return the fitted object.
}
\description{
Convenience plotting functions for component spectra, component maps, loss
curves, and data-versus-reconstruction diagnostics.
}
\examples{
\dontrun{
demo <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
fit <- spax_nmf(demo$matrix, k = 3, niter = 50)

plot(fit, type = "spectra", wavelength = demo$wavelength)
plot(fit, type = "maps", nx = demo$nx, ny = demo$ny)
plot(fit, type = "loss")
plot_reconstruction(fit, demo$matrix, n = 4, wavelength = demo$wavelength)
}
}
