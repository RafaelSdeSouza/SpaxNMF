\name{cube_to_matrix}
\alias{cube_to_matrix}
\alias{matrix_to_cube}
\alias{component_map}
\alias{component_spectrum}
\alias{component_reconstruction}
\alias{cube_metadata}
\title{Helpers for Cube Reshaping and Map Extraction}
\usage{
cube_to_matrix(cube)

matrix_to_cube(x, nx, ny, metadata = NULL)

component_map(fit, nx, ny, component = 1)

component_spectrum(fit, component = 1)

component_reconstruction(fit, component = 1, nx = NULL, ny = NULL)

cube_metadata(object)
}
\arguments{
\item{cube}{Numeric 3D array, or a list-like object containing a cube in
\code{$imDat} or \code{$cube}.}

\item{x}{Numeric matrix with one row per spatial pixel.}

\item{nx}{Number of x-axis pixels.}

\item{ny}{Number of y-axis pixels.}

\item{metadata}{Optional metadata list to attach to the returned cube.}

\item{fit}{Result from \code{spax_nmf()}.}

\item{component}{Component index to extract.}

\item{object}{Matrix, cube, or fitted object.}
}
\value{
\code{cube_to_matrix()} returns a numeric matrix with \code{nx * ny} rows and
\code{n_wave} columns.

\code{matrix_to_cube()} returns a numeric 3D array.

\code{component_map()} returns a numeric matrix of size \code{nx x ny}.

\code{component_spectrum()} returns a numeric vector.

\code{component_reconstruction()} returns a numeric matrix, or a numeric 3D
array when \code{nx} and \code{ny} are supplied.

\code{cube_metadata()} returns a metadata list, or \code{NULL} when no metadata
is present.
}
\description{
Utility functions for moving between cube and matrix representations and for
extracting component spectra, abundance maps, metadata, and single-component
reconstructions from a fitted model.
}
\examples{
demo <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
x <- cube_to_matrix(demo$cube)
rebuilt <- matrix_to_cube(x, demo$nx, demo$ny)
dim(rebuilt)
\dontrun{
fit <- spax_nmf(x, k = 3, niter = 50)
img <- component_map(fit, demo$nx, demo$ny, 2)
spec2 <- component_spectrum(fit, 2)
comp2 <- component_reconstruction(fit, 2, nx = demo$nx, ny = demo$ny)
cube_metadata(comp2)
}
}
