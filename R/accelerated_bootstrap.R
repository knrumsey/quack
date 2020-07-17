#' Estimation of the "Accelration Constant"
#'
#' Estimates the accleration constant using a simple jackknife estimator (Efron 1987).
#'
#' @param x the data vector
#' @param est_theta function to estimate parameter
#' @examples
#' est_accle(rnorm(30), function(z) quantile(z, 0.75))
#' @export
est_accel <- function(x, est_theta, ...){
  n <- length(x)
  tmp_indx <- 1:n
  theta_hat <- est_theta(x, ...)
  I <- rep(NA, n)
  for(i in 1:n){
    tmp_indx <- (1:n)[-i]
    theta_j <- est_theta(x[tmp_indx], ...)
    I[i] <- (n-1)*(theta_hat - theta_j)
  }
  return((sum(I^3)/sum(I^2)^1.5)/6)
}

#' Accelerated Bootstrap
#'
#' Flexible code for constructing a confidence interval using accelrated bootstrap
#'
#' @param x the data vector
#' @param est_theta function to estimate parameter
#' @param B number of bootstrap samples
#' @param alpha significance level. Default is 0.042 (because 0.05 is arbitrary)
#' @param a acceleration constant. See details
#' @param ... additional parameters passed to est_theta
#' @details By default, the acceleration constant is a=NULL which leads to estimation of the constant
#' using the est_accel() function. This can be skipped by specifying a particular value for a. Note that
#' a=0 corresponds to the usual "quantile" bootstrap.
#'
#' The code uses a special variable called 'tmp_indx' which can be exploited to handle more complicated
#' cases using the dynGet() function. See the example below for more details.
#' @examples
#' #Simulate data
#' x <- rnorm(50)
#' y <- x + rnorm(50)
#' #Simple example
#' boot_accel(x, function(z) quantile(z, 0.75) , B=1e4)
#'
#' #Bootstrap for correlation using dynGet()
#' my_cor <- function(x, y, indx=dynGet('tmp_indx')){
#'    cor(x, y[indx])
#' }
#' boot_accel(x, my_cor, y=y)
#' @export
boot_accel <- function(x, est_theta, B=1e4, alpha=0.042, a=NULL, ...){
  n <- length(x)
  tmp_indx <- 1:n
  theta_hat <- est_theta(x, ...)
  if(is.null(a)){
    a <- est_accel(x, est_theta, ...)
  }

  #Begin Bootstrap
  theta_boot <- rep(NA, B)
  for(b in 1:B){
    tmp_indx <- sample(n, length(x), replace=TRUE)
    theta_boot[b] <- est_theta(x[tmp_indx], ...)
  }

  #Return accel-boot interval
  z0 <- qnorm(mean(theta_boot <= theta_hat))
  zu <- qnorm(c(alpha/2, 1-alpha/2))
  u_adj <- pnorm(z0 + (z0+zu)/(1-a*(z0+zu)))
  return(quantile(theta_boot, u_adj))
}







