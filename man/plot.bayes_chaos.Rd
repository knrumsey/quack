% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chaos.R
\name{plot.bayes_chaos}
\alias{plot.bayes_chaos}
\title{Plot Method for class bayes_chaos}
\usage{
\method{plot}{bayes_chaos}(x, ...)
}
\arguments{
\item{x}{An object returned by the \code{bayes_chaos()} function.}

\item{...}{additional arguments passed to \code{plot}}
}
\description{
See \code{bayes_chaos()} for details.
}
\details{
Plot function for bayes_chaos object.
}
\examples{
X <- lhs::maximinLHS(500, 3)
y <- apply(X, 1, duqling::dms_simple) + rnorm(500, 0, 0.1)
fit <- bayes_chaos(X, y)
plot(fit)
}
\references{
Shao, Q., Younes, A., Fahs, M., & Mara, T. A. (2017). Bayesian sparse polynomial chaos expansion for global sensitivity analysis. Computer Methods in Applied Mechanics and Engineering, 318, 474-496.
}
