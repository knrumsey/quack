#' Customized "smart" (robust) LHS function
#'
#' A default lhs function to use (to avoid accidental crashing of maximinLHS)
#'
#' @param n the number of rows or samples
#' @param k the number of columns or parameters/variables
#' @param maximin_n the number of rows or samples in the original maximin design
#' @param ... additional arguments passed to randomLHS (n > 50k) or maximinLHS (n <= 50k)
#' @return a latin hypercube design
#' @details Tries to use a maxmin LHS design when possible (n <= maximin_n = 1000).
#' If n is greater than 100,000, we simply return a random LHS.
#' In between these cases, we build an initial maximin design of size maximin_n and then augment this design with augmentLHS.
#' @examples
#' x <- myLHS(1000, 8)
#' y <- apply(x, 1, quack::ff_borehole)
#' @export
smartLHS <- function(n, k, maximin_n=1e3, ...){
  if(n > 1e5){
    return(randomLHS(n, k, ...))
  }
  if(n <= maximin_n){
    return(maximinLHS(n, k, ...))
  }
  n2 <- n - maximin_n
  x <- maximinLHS(maximin_n, k, ...)
  return(augmentLHS(x, n2))
}
