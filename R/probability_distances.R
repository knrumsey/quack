#' Lp Distance Between CDFs for Multivariate Samples
#'
#' @param x a vector of samples from P
#' @param y a vector of samples from Q
#' @param p norm order (p=Inf returns KS distance, p=1 returns Kantorovich distance)
#' @param n_pts number of evaluation points IN EACH DIMENSION. Run time is O(n_pts^ncol(x)).
#' @details Approximates distance between CDFs. Note that the performance rapidly deteriorates for large ncol(x)
#' @examples
#' x <- matrix(rnorm(300), ncol=3)
#' y <- matrix(rgamma(180, 3, 2) - 1.5, ncol=3)
#' cdf_dist_lp(x, y, p=1, n_pts=10)
#' cdf_dist_lp(x, y, p=2, n_pts=10)
#' cdf_dist_lp(x, y, p=Inf, n_pts=10)
#'
#' @export
cdf_dist_lp <- function(x, y, p=1, n_pts=30){
  if(!is.matrix(x) | !is.matrix(y)){
    x <- matrix(x, ncol=1)
    y <- matrix(y, ncol=1)
  }
  nx <- nrow(x)
  ny <- nrow(y)
  q <- ncol(x)
  if(q != ncol(y)){
    stop("x and y must be matrices with same number of columns.")
  }
  bounds <- matrix(NA, nrow=2, ncol=q)
  for(i in 1:q){
    bounds[,i] <- range(c(x[,i], y[,i]))
  }

  indx <- generate_cartesian_product(n_pts, q)
  ks <- -Inf
  res <- 0
  for(i in 1:nrow(indx)){
    flagx <- rep(TRUE, nx)
    flagy <- rep(TRUE, ny)
    for(j in 1:q){
      ii <- indx[i,j]
      z <- bounds[1,j] + ((ii - 1) / (n_pts - 1)) * diff(bounds[,j])
      flagx <- flagx & (x[,j] <= z)
      flagy <- flagy & (y[,j] <= z)
    }
    Fx <- sum(flagx)/nx
    Fy <- sum(flagy)/ny
    curr <- abs(Fx-Fy)
    res <- res + curr^p/n_pts^q
    ks <- max(ks, curr)
  }

  if(p == Inf){
    return(ks)
  }
  return(res^(1/p))
}

#' KS Distance Between CDFs for Multivariate Samples
#'
#'
#' @param x a vector of samples from P
#' @param y a vector of samples from Q
#' @param n_pts number of evaluation points IN EACH DIMENSION. Run time is O(n_pts^ncol(x)).
#' @details Multivariate KS distance
#' @examples
#' x <- matrix(rnorm(300), ncol=3)
#' y <- matrix(rgamma(180, 3, 2) - 1.5, ncol=3)
#' ks_dist(x, y)
#'
#' @export
ks_dist <- function(x, y, n_pts=30){
  cdf_dist_lp(x, y, p=Inf, n_pts=n_pts)
}


#' Kantorovich Distance Between CDFs for Multivariate Samples
#'
#'
#' @param x a vector of samples from P
#' @param y a vector of samples from Q
#' @param n_pts number of evaluation points IN EACH DIMENSION. Run time is O(n_pts^ncol(x)).
#' @details Multivariate KS distance
#' @examples
#' x <- matrix(rnorm(300), ncol=3)
#' y <- matrix(rgamma(180, 3, 2) - 1.5, ncol=3)
#' kantorovich_dist(x, y)
#'
#' @export
kantorovich_dist <- function(x, y, n_pts=30){
  cdf_dist_lp(x, y, p=1, n_pts=n_pts)
}

generate_cartesian_product <- function(n, p) {
  # Create a list of p sequences from 1 to n
  args_list <- replicate(p, 1:n, simplify = FALSE)

  # Use do.call to pass the list to expand.grid
  result <- do.call(expand.grid, args_list)

  return(result)
}


