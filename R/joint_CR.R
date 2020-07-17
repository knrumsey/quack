#' Mahalanobis Distance
#'
#' Computes the mahalanobis distance between two vectors given a precision matrix S
#'
#' @param x a vector
#' @param y an optional vector. The zero vector by default
#' @param S a precision matrix of dimension equal to length(x)
#' @examples x <- rnorm(5); y <- rnorm(5); S <- GPfit::corr_matrix((1:5)/6, beta=0.2)
#' d_mahal(x, y, S)
d_mahal <- function(x, y=NULL, S=diag(rep(1,length(x)))){
  if(is.null(y)) y <- rep(0, length(x))
  n <- length(x)
  v <- matrix(x-y, nrow=1)
  return(as.numeric(v%*%S%*%t(v)))
}

#' Joint Confidence Regions using Mahalanobis Distance
#'
#' Finds the minimal convex hull, with respect to Mahalanobis distance, which contains (1-alpha)x100% of the points.
#'
#' @param X a matrix - each row represents a sample
#' @param alpha the confidence level is 1-alpha. Defaults to 0.042 (because 0.05 is so arbitrary)
#' @param two_dim logical. If TRUE, then the convex hull surrounding the interior points in the first two dimensions is returned.
#'  If FALSE, then we return the entire set of interior points (for arbitrary dimension).
#' @return the samples laying on the convex hull which contains (1-alpha)100% of the samples.
#' @details Note that a set of N-dimensional points is always returned. When two_dim=TRUE,
#' this set of points represents the points on the convex hull of the JCR - but only in the first two dimensions.
#' This is a limitation of the built-in chull() function. Future versions may extend the functionality.
#' For now, one can operate in higher dimensions by setting two_dim = FALSE to return the full set of
#' (1-alpha)*N points which are in the interior of the JCR.
#' @examples
#' x <- rnorm(1000); y <- x + rnorm(1000, 0, .5)
#' X <- cbind(x, y)
#' CR <- joint_CR(X)
#' plot(X)
#' polygon(CR, border='blue', lwd=2)
#'
#' #A point inside of JCR
#' inside_CR(c(0,0), CR)
#' #A point on the convex hull of JCR
#' lambda <- 0.3
#' inside_CR(lambda*CR[3,] + (1-lambda)*CR[4,], CR)
#' #A point outside of JCR
#' inside_CR(c(-3, 7), CR)
#' @export
joint_CR <- function(X, alpha=0.042, two_dim=TRUE){
  d <- apply(X, 1, d_mahal,
             y=apply(X, 2, mean),
             S=solve(cov(X)))
  interior <- X[which(d <= quantile(d, 1-alpha)),]
  if(two_dim){
    return(interior[chull(interior),])
  }else{
    return(interior)
  }
}

#' Is Inside Credible Region
#'
#' Checks to see whether a 2-D query point is inside, on or outside the convex hull
#' representing an empirical credible region
#'
#' @param q a query point
#' @param CR the vertices of the convex hull returned by joint_CR (with two_dim = TRUE)
#' @return Returns +1, 0, or -1 for inside JCR, on JCR or outside JCR respectively.
#' @details Only implemented (currently) for 2 dimensions.
#' @export
inside_CR <- function(q, CR){
  if(dim(CR)[2] != 2) warning('Only meant for 2-D points')
  N <- dim(CR)[1]
  ref <- apply(CR, 2, mean)
  res <- 1
  for(i in 1:N){
    res <- min(res, get_position(q, CR[i,], CR[(i%%N)+1,]))
  }
  return(res)
}

get_position <- function(q, p0, p1, tol=1e-6){
  #Center at 0
  q <- q-p0
  p1 <- p1-p0
  return(-sign(round(p1[1]*q[2] - p1[2]*q[1], digits=9)))
}


