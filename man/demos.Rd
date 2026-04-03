\name{simulate_ifu_cube}
\alias{simulate_ifu_cube}
\alias{demo_ifu_unmix}
\title{Synthetic IFU Demos}
\usage{
simulate_ifu_cube(nx = 20, ny = 20, n_wave = 120, noise = 0.01)

demo_ifu_unmix(
  k = 3,
  nx = 18,
  ny = 18,
  n_wave = 120,
  niter = 300,
  plot = interactive()
)
}
\arguments{
\item{nx}{Number of x-axis pixels.}

\item{ny}{Number of y-axis pixels.}

\item{n_wave}{Number of wavelength channels.}

\item{noise}{Standard deviation of additive Gaussian noise.}

\item{k}{Number of components.}

\item{niter}{Number of optimization iterations.}

\item{plot}{Logical; if \code{TRUE}, display a compact base-R summary plot.}
}
\value{
\code{simulate_ifu_cube()} returns a list with the simulated \code{cube}, its
matrix representation, the true component spectra, abundance maps, wavelength
grid, and spatial dimensions.

\code{demo_ifu_unmix()} returns a list with the simulated data and fitted model.
}
\description{
Helpers for building reproducible IFU examples and package documentation.
}
\examples{
demo <- simulate_ifu_cube()
str(demo, max.level = 1)
\dontrun{
result <- demo_ifu_unmix(plot = TRUE, niter = 150)
}
}
