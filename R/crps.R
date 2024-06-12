#' Estimate CRPS(F, y) given samples from F
#'
#' Returns the an estimator between the INT and PWM estimators of Zamo
#'
#' @param y the true scalar value
#' @param x a vector of samples from F
#' @param w parameter (between 0 and 1). When w=0, the estimator is unbiased ("fair" and "PWM" from Zamo & Naveau 2017). When w=1, the estimator has lower variance ("INT" and "NRG" estimator from Zamo & Naveau 2017).
#' @references Zamo, M., & Naveau, P. (2018). Estimation of the continuous ranked probability score with limited information and applications to ensemble weather forecasts. Mathematical Geosciences, 50(2), 209-234.
#' @examples
#' F1 <- function(x) pnorm(x)
#' I1 <- function(x, y) as.numeric(x >= y)
#' Q <- function(x, y) (F1(x) - I1(x, y))^2
#' y <- 0.1
#' x <- rnorm(100)
#' integrate(Q, lower=-10, upper=10, y=y)$value
#' crps(y, x, w=0) #default w
#' crps(y, x, w=0.5)
#' crps(y, x, w=1)

#' @export
crps <- function(y, x, w=0){
  M <- length(x)
  term1 <- mean(abs(x-y))
  if(M <= 6500){
    term2 <- sum(outer(x, x, function(a, b) abs(a-b))) # Fastest way for small M
  }else{
    idx <- unlist(lapply(2:M, function(i) 1:i))
    term2 <- 2*sum(abs(x[rep(2:M, 2:M)] - x[idx]))     # Faster for big M
  }
  res <- term1 - (1 - w/M)*term2/(2*M*(M-1))
  return(res)
}
