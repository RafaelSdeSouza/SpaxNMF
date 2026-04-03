\name{SpaxNMF}
\alias{SpaxNMF}
\alias{SpaxNMF-package}
\docType{package}
\title{SpaxNMF: Regularized Non-negative Matrix Factorization for IFU Cubes}
\description{
Tools for reshaping IFU cubes, fitting smooth non-negative factorization
models, and building reproducible examples for package documentation and
comparison experiments.
}
\details{
The package is designed around a simple workflow:
\enumerate{
\item simulate or ingest a spectral cube;
\item reshape it with \code{cube_to_matrix()};
\item fit a low-rank unmixing model with \code{spectral_unmix()};
\item inspect component spectra and abundance maps.
}
}
