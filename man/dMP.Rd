% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MP.R
\name{dMP}
\alias{dMP}
\title{Density of the Moment Penalization Prior}
\usage{
dMP(x, w1 = 1, w2 = 1, log = FALSE, normalized = FALSE, ...)
}
\arguments{
\item{x}{a vector of length p}

\item{w1}{normalized penalty associated with second moment. Default is 1}

\item{log}{logical. Should density be returned on a log scale?}

\item{normalized}{logical. Should density be normalized (default is FALSE)}

\item{...}{additional parameters passed to get_constMP (if norm=TRUE)}
}
\value{
returns the density of the MP(w1, w2) prior
}
\description{
Returns the unnormalized density of the MP prior with parameters w1 and w2
}
\examples{
x <- rnorm(10)
dMP(x, w1=5, w2=2)
X <- matrix(rnorm(10*30), nrow=30)
apply(X, 1, dMP, w1=5, w2=2)
}
