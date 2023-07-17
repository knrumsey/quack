#' Maximal Couplings
#'
#' A function to sample from the maximal coupling of two distributions (using Johnson (1998))
#'
#' @param n number of requested samples
#' @param d1 density of distribution 1
#' @param d2 density of distribution 2
#' @param g1 generator for distribution 1
#' @param g2 generator for distribution 2
#' @return An nx2 matrix of draws where the ith column has (marginal) distribution i. The Pr(X[,1] == X[,2]) is maximized.
#' @details Uses the gamma-coupling algorithm of Johnson 1998
#' @export
#' @examples
#' d1 <- function(x) dnorm(x)
#' d2 <- function(x) dunif(x, -3, 3)
#' g1 <- function()  rnorm(1)
#' g2 <- function()  runif(1, -3, 3)
#' maximal_coupling(1000, d1, d2, g1, g2)
maximal_coupling <- function(n, d1, d2, g1=NULL, g2=NULL){
  X <- matrix(NA, nrow=n, ncol=2)
  for(i in 1:n){
    x <- g1()
    w <- runif(1, 0, d1(x))
    if(w <= d2(x)){
      X[i,] <- c(x, x)
    }else{
      flag <- TRUE
      while(flag){
        y <- g2()
        w <- runif(1, 0, d2(y))
        if(w > d1(y)){
          X[i,] <- c(x, y)
          flag <- FALSE
        }
      }
    }
  }
  return(X)
}



