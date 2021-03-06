% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fast_functions.R
\name{fast_inv.FM_approx}
\alias{fast_inv.FM_approx}
\title{Fast Matrix Algebra for BMC}
\usage{
\method{fast_inv}{FM_approx}(FM, phi = NULL, tau = NULL, ...)
}
\arguments{
\item{FM}{an object of class "FM_approx"}

\item{phi}{current value of phi (can be omitted if phi is stored in FM object)}

\item{tau}{current value of tau (can be omitted if tau is stored in FM object)}

\item{...}{appropriate arguments for the corresponding class.}
}
\value{
Returns an approximation of the inverse of the covariance matrix
}
\description{
Function for (approximately) computing the inverse of a BMC Covariance matrix.
Near-quadratic runtime leads to significant time-savings when covariance matrix is moderate to large.
}
\examples{
n <- 100
R <- corr_matrix(seq(0,1,length.out=n), beta=0.2)
FM <- fast_process(R, method='approx', control=list(tol=1e-3))
#Compute determinant (approximately)
fast_det(FM, phi=1.5, tau=0.2)
#Compute inverse (approximately)
fast_inv(FM, phi=1.5, tau=0.2)
}
