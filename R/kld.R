#' KL Divergence based on samples from P and Q
#'
#' Returns the KLD estimator of Perez-Cruz (2008)
#'
#' @param x a vector of samples from P
#' @param y a vector of samples from Q
#' @param type See details
#' @param lambda parameter (between 0 and 1) for the Population Stability Index (ignored unless type = 4)
#' @param n_mc number of monte carlo samples used when type = 3 or 4
#' @details Set type = 1 for D(P||Q), type = 2 for symmetrized "Jeffreys Divergence", type = 3 for symmetric "Jensen-Shannon Distance", type = 4 for Population Stability Index (with parameter lambda between 0 and 1). When type = 3 or 4, n_mc must be set for approximation. When type = 4, lambda must be specified. lambda = 1/2 is equivalent to type = 3.
#' @examples
#' x <- rnorm(100)
#' y <- rgamma(200, 3, 2) - 1.5
#' kld(x, y)
#' @export
kld <- function(x, y, type=1, lambda=1/2, n_mc=1){
  if(type == 1){
    return(compute_kld(x, y))
  }
  if(type == 2){
    D1 <- compute_kld(x, y)
    D2 <- compute_kld(y, x)
    return((D1+D2)/2)
  }
  if(type == 3){
    lambda <- 1/2
  }
  if(type >= 3){
    nx <- length(x)
    ny <- length(y)
    N <- min(nx, ny)
    D <- rep(NA, n_mc)
    z <- rep(NA, N)
    for(i in 1:n_mc){
      z0 <- sample(c(x, y), size=N,
                  replace=FALSE,
                  prob=c(rep(lambda/nx, nx),
                          rep((1-lambda)/ny, ny)))
      #'I don't know why this is needed,
      #'but estimator gets weird without it.
      #'I guess it's because the samples have
      #'to be independent, and the way i'm doing
      #'the mixture the samples are dependent
      #'(e.g., z dependent on x)
      z <- z0 + rnorm(N, 0, 2*diff(range(z0))/N)
      D1  <- compute_kld(x, z)
      D2  <- compute_kld(y, z)
      lambda*D1 + (1-lambda)*D2
      D[i] <- lambda*D1 + (1-lambda)*D2
    }
    return(D)
  }
  stop("Not a valid type")
  return(FALSE)
}

compute_kld <- function(x, y){
  x <- sort(x)
  y <- sort(y)
  rxy <- range(c(x, y))
  eps <- min(diff(x))/2
  x <- sort(union(x, rxy + eps*c(-1, 1)))
  y <- sort(union(y, rxy + eps*c(-1, 1)))

  #xx should be univariate, z is a sorted vector
  pz <- function(xx, z){
    if(xx < min(z)) return(0)
    if(xx > max(z)) return(1)
    ind_lb <- which.min(abs(xx-z) + 2*max(diff(z))*(xx < z))
    ind_ub <- which.min(abs(xx-z) + 2*max(diff(z))*(xx > z))

    z_lb <- z[ind_lb]
    z_ub <- z[ind_ub]

    U_lb <- (z_lb > z) + 0.5*(z_lb == z)
    U_ub <- (z_ub > z) + 0.5*(z_ub == z)

    if(ind_lb != ind_ub){
      Fz <- mean(U_lb) + mean(U_ub - U_lb)/(z_ub - z_lb)*(xx - z_lb)
    }else{
      Fz <- mean(U_lb)
    }
    return(Fz)
  }

  D <- 0
  n <- length(x)
  Dx <- unlist(lapply(1:n, function(i, x, eps) log(pz(x[i], x) - pz(x[i]-eps, x)), x=x, eps=eps))
  Dy <- unlist(lapply(1:n, function(i, x, y, eps) log(pz(x[i], y) - pz(x[i]-eps, y)), x=x, y=y, eps=eps))
  D <- mean(Dx - Dy) - 1
  return(D)
}

x <- rexp(1000)
y <- rnorm(1000, 3, 2)
z <- rexp(1000)*(B <- rbinom(1000, 1, 1/2)) + (1-B)*rnorm(1000, 3, 2)
1/2*kld(x, z) + 1/2*kld(y, z)

kld(x, y, type=3)


