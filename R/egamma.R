#' @name egamma
#' @rdname egamma
#'
#' @title Extended Gamma Distribution
#'
#' The EGamma(alpha, beta, kappa) is best described by its CDF which has the form F(x) = wF_1(x) + (1-w)F_2(x)
#' where w = beta^-alpha/(beta^-alpha - (beta + kappa)^-alpha) > 1. F_1 is the CDF of a gamma(alpha, beta) distribution
#' and F_2 is the CDF of a gamma(alpha, beta + kappa) distribution.
#'
#' @param x a vector of length p
#' @param q a vector of quantiles
#' @param p a vector of probabilities
#' @param n number of observations. If `length(n) > 1`, the length is taken to be the number required.
#' @param k the moment desired
#' @param alpha shape of the gamma distribution
#' @param beta rate of the gamma distribution
#' @param kappa extension parameter (positive)
#' @param log logical. Should density/probability be returned on a log scale?
#' @param lower.tail logical; if TRUE (default) probabilities are P(X <= x), otherwise P(X > x).
#' @param ... optional additional arguments passed to `optim()`
#' @details `degamma` returns the density, `pegamma` returns the distribution function, `qegamma` returns the quantile function, `regamma` generates random samples (using inverse transform) and `megamma` returns the moments of the distribution.
NULL

#' @rdname egamma
#' @examples
#' n <- 100 #Number of observations
#' X <- regamma(n, 3, 1, 2)
#' hist(X, breaks=20, freq=FALSE)
#' curve(degamma(x, 3, 1, 2), add=TRUE)
#' @export
degamma <- function(x, alpha, beta, kappa, log=FALSE){
  w <- beta^-alpha/(beta^-alpha - (beta + kappa)^-alpha)
  d <- w*dgamma(x, alpha, beta) + (1-w)*dgamma(x, alpha, beta + kappa)
  if(log){
    return(log(d))
  }else{
    return(d)
  }
}

#' @rdname egamma
pegamma <- function(q, alpha, beta, kappa, log=FALSE, lower.tail=TRUE){
  w <- beta^-alpha/(beta^-alpha - (beta + kappa)^-alpha)
  p <- w*pgamma(q, alpha, beta) + (1-w)*pgamma(q, alpha, beta + kappa)
  p <- p + (1 - 2*p)*(1-lower.tail)
  if(log){
    return(log(p))
  }else{
    return(p)
  }
}

#' @rdname egamma
qegamma <- function(p, alpha, beta, kappa, ...){
  w <- beta^-alpha/(beta^-alpha - (beta + kappa)^-alpha)
  l <- pmin(qgamma(p, alpha, beta), 0)
  u <- qgamma(1 - (1-p)/w, alpha, beta)
  f2opt <- function(zz, pp) (pegamma(zz, alpha, beta, kappa) - pp)^2
  q <- rep(NA, length(p))
  for(i in seq_along(p)){
    if(p[i] <= 0){
      q[i] <- 0
    }else{
      if(p[i] > 1){
        q[i] <- 1
      }else{
        q[i] <- optim((l[i]+u[i])/2, fn=f2opt,
                      method="Brent", lower=l[i], upper = u[i],
                      pp = p[i], ...)$par
      }
    }

  }
  return(q)
}

#' @rdname egamma
regamma <- function(n, alpha, beta, kappa, ...){
  if(length(n) > 1) n <- length(n)
  u <- runif(n)
  x <- qegamma(u, alpha, beta, kappa, ...)
  return(x)
}

#' @rdname egamma
megamma <- function(k=1, alpha, beta, kappa, ...){
  w <- beta^-alpha/(beta^-alpha - (beta + kappa)^-alpha)
  term1 <- gamma(alpha + k)/gamma(alpha)/beta^k
  term2 <- w*(1-((w-1)/w)^(k/alpha + 1))
  return(term1*term2)
}


