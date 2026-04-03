\name{basis}
\alias{basis}
\alias{basis.spax_nmf}
\alias{coef.spax_nmf}
\alias{fitted.spax_nmf}
\alias{predict.spax_nmf}
\alias{residuals.spax_nmf}
\alias{summary.spax_nmf}
\alias{print.summary.spax_nmf}
\title{Standard Model Methods for SpaxNMF Fits}
\usage{
basis(object, ...)

\method{coef}{spax_nmf}(object, ...)

\method{fitted}{spax_nmf}(object, nx = NULL, ny = NULL, ...)

\method{predict}{spax_nmf}(
  object,
  newdata = NULL,
  type = c("reconstruction", "spatial", "cube"),
  nx = NULL,
  ny = NULL,
  lr = 0.05,
  niter = 500,
  cuda = FALSE,
  verbose = FALSE,
  ...
)

\method{residuals}{spax_nmf}(object, x, nx = NULL, ny = NULL, ...)

\method{summary}{spax_nmf}(object, ...)

\method{print}{summary.spax_nmf}(x, ...)
}
\arguments{
\item{object}{Result from \code{spax_nmf()}.}

\item{...}{Unused.}

\item{nx}{Optional x-axis size for cube output.}

\item{ny}{Optional y-axis size for cube output.}

\item{newdata}{Optional new matrix or cube.}

\item{type}{Prediction target: \code{"reconstruction"}, \code{"spatial"}, or
\code{"cube"}.}

\item{lr}{Learning rate used when estimating abundance weights for new data.}

\item{niter}{Number of optimization iterations for new data.}

\item{cuda}{Logical; if \code{TRUE}, request CUDA execution.}

\item{verbose}{Logical; if \code{TRUE}, emit progress every 100 iterations.}

\item{x}{For \code{residuals()}, the numeric matrix originally used for the
fit. For \code{print.summary.spax_nmf()}, the summary object.}
}
\value{
\code{basis()} returns a wavelength-by-component matrix.

\code{coef()} returns a component-by-spaxel matrix.

\code{fitted()} returns the reconstruction matrix, or a cube when
\code{nx} and \code{ny} are supplied.

\code{predict()} returns either abundance weights, a reconstruction matrix,
or a reconstructed cube.

\code{residuals()} returns residuals as a matrix or cube.

\code{summary()} returns a list with class \code{"summary.spax_nmf"}.
}
\description{
Convenience methods that make \code{spax_nmf()} objects behave more like
standard decomposition objects from \code{prcomp()} and NMF-style interfaces.
}
\examples{
fake <- list(
  spatial = matrix(runif(12), 4, 3),
  spectra = matrix(runif(15), 3, 5),
  reconstruction = matrix(runif(20), 4, 5),
  loss = c(3, 2, 1)
)
class(fake) <- "spax_nmf"

basis(fake)
coef(fake)
fitted(fake)
summary(fake)

x <- matrix(runif(20), 4, 5)
residuals(fake, x = x)
}
