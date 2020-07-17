#' Fast Matrix Algebra for BMC
#'
#' Function for obtaining the determinant in linear time - after pre-processing - for covariance matrices of suitable structure
#'
#' @param FM an object created by the fast_process() function
#' @param phi current value of phi (can be omitted if phi is stored in FM object)
#' @param tau current value of tau (can be omitted if tau is stored in FM object)
#' @param log logical. Should log-determinant be returned? Default is TRUE.
#' @return Returns the (log) determinant of the covariance matrix.
#' @examples
#' n <- 100
#' R <- GPfit::corr_matrix(seq(0,1,length.out=n), beta=0.2)
#' FM <- fast_process(R, method='approx', control=list(tol=1e-9))
#' #Compute determinant (approximately)
#' fast_det(FM, phi=1.5, tau=0.2)
#' #Compute inverse (approximately)
#' fast_inv(FM, phi=1.5, tau=0.2)
#' @export
fast_det <- function(FM, phi=NULL, tau=NULL, log=TRUE){
  phi <- ifelse(length(phi) > 0, phi, FM$phi)
  tau <- ifelse(length(tau) > 0, phi, FM$tau)
  if(is.null(phi) || is.null(tau)){
    stop('Could not find phi and/or tau')
  }
  out <- sum(log(phi*FM$r + tau))
  if(log){
    return(out)
  }else{
    return(exp(out))
  }
}

#' Fast Matrix Algebra for BMC
#'
#' Function for computing the inverse of a BMC Covariance matrix. Leads to speedup of roughly 2 when matrix dimension is moderate or large.
#'
#' @param FM an object of class "FM_fast"
#' @param phi current value of phi (can be omitted if phi is stored in FM object)
#' @param tau current value of tau (can be omitted if tau is stored in FM object)
#' @return Returns the inverse of the covariance matrix
#' @export
fast_inv.FM_fast <- function(FM, phi=NULL, tau=NULL){
  phi <- ifelse(length(phi) > 0, phi, FM$phi)
  tau <- ifelse(length(tau) > 0, tau, FM$tau)
  if(is.null(phi) || is.null(tau)){
    stop('Could not find phi and/or tau')
  }
  W <- sqrt(tau/(phi*FM$r + tau))
  tmpcntr <- 0 #Hacky temp counter for speed
  QW <- apply(FM$Q, 2, FUN=function(z, W){tmpcntr <<- tmpcntr + 1; z*W[tmpcntr]}, W=W)
  out <- tcrossprod(QW)/tau
  return(out)
}


#' Fast Matrix Algebra for BMC
#'
#' Function for (approximately) computing the inverse of a BMC Covariance matrix.
#' Near-quadratic runtime leads to significant time-savings when covariance matrix is moderate to large.
#'
#' @param FM an object of class "FM_approx"
#' @param phi current value of phi (can be omitted if phi is stored in FM object)
#' @param tau current value of tau (can be omitted if tau is stored in FM object)
#' @return Returns an approximation of the inverse of the covariance matrix
#' @examples
#' n <- 100
#' R <- GPfit::corr_matrix(seq(0,1,length.out=n), beta=0.2)
#' FM <- fast_process(R, method='approx', control=list(tol=1e-9))
#' #Compute determinant (approximately)
#' fast_det(FM, phi=1.5, tau=0.2)
#' #Compute inverse (approximately)
#' fast_inv(FM, phi=1.5, tau=0.2)
#' @export
fast_inv.FM_approx <- function(FM, phi=NULL, tau=NULL){
  phi <- ifelse(length(phi) > 0, phi, FM$phi)
  tau <- ifelse(length(tau) > 0, tau, FM$tau)
  if(is.null(phi) || is.null(tau)){
    stop('Could not find phi and/or tau')
  }
  N0 <- FM$N0
  W <- sqrt(1/(phi*FM$r[1:N0] + tau))
  tmpcntr <- 0 #Hacky temp counter for speed
  QW0 <- apply(FM$Q[,1:N0], 2, FUN=function(z, W){tmpcntr <<- tmpcntr + 1; z*W[tmpcntr]}, W=W)
  out <- tcrossprod(QW0) + FM$Z/tau
  return(out)
}


#' Fast Matrix Algebra for BMC
#'
#' Function for (approximately) computing the inverse of a BMC Covariance matrix.
#' Near-quadratic runtime leads to significant time-savings when covariance matrix is moderate to large.
#'
#' @param FM an object of class "FM_ts"
#' @param phi current value of phi (can be omitted if phi is stored in FM object)
#' @return Returns an approximation of the inverse of the covariance matrix
#' @examples
#' n <- 100
#' R <- GPfit::corr_matrix(seq(0,1,length.out=n), beta=0.2)
#' FM <- fast_process(R, method='approx', control=list(tol=1e-9))
#' #Compute determinant (approximately)
#' fast_det(FM, phi=1.5, tau=0.2)
#' #Compute inverse (approximately)
#' fast_inv(FM, phi=1.5, tau=0.2)
#' @export
fast_inv.FM_ts <- function(FM, phi=NULL){
  phi <- ifelse(length(phi) > 0, phi, FM$phi)
  if(is.null(phi)){
    stop('Could not find phi and/or tau')
  }
  N <- length(FM$r)
  out <- matrix(0, nrow=N, ncol=N)
  for(p in 0:FM$P){
    out <- out + (-1)^p*(phi-FM$a)^p*FM$C[,,p+1]
  }
  return(out)
}
