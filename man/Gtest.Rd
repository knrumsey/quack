% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gtest.R
\name{Gtest}
\alias{Gtest}
\title{G Test}
\usage{
Gtest(arr, df = NA)
}
\arguments{
\item{arr}{the true scalar value}

\item{df}{If not specified, test uses the standard mutual independence assumption: prod(dims) - sum(dims-1) - 1.}
}
\description{
Performs G-test for an arbitrary D-dimensional array
}
\examples{
data_2d <- cbind(rpois(3, lambda=10), rpois(3, lambda=10), rpois(3, lambda=10))
fisher.test(data_2d, alternative="greater")
chisq.test(data_2d)
Gtest(data_2d)

data_3d <- array(rpois(2*3*4, lambda=10), dim=c(2,3,4))
Gtest(data)
}
\references{
Fill this in.
}
