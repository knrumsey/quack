#' Fast Matrix Algebra for BMC
#'
#' Functions for obtaining matrix inversions and determinants quickly in the context of BMC.
#' Processes matrices of the form Sigma = phi x R + tau x I
#'
#' @param R fixed correlation matrix
#' @param phi Optional. Indicates that phi is fixed throughout the analysis.
#' @param tau Optional. Indicates that tau is fixed throughout the analysis. Required when method == 'ts'
#' @param method Must be 'fast', 'approx', or 'ts'.
#' @param control a list of control parameters. Will be ignored when method == 'fast'.
#' Method 'approx' needs control$tol (1e-9) recomended, and method 'ts' needs an expansion value control$a and an expansion order control$P
#' @return returns an object of class "FM_method" containing the information needed for fast computation of matrix inverse and determinant
#' @examples
#' n <- 100
#' R <- corr_matrix(seq(0,1,length.out=n), beta=0.2)
#' FM <- fast_process(R, method='approx', control=list(tol=1e-9))
#' #Compute determinant (approximately)
#' fast_det(FM, phi=1.5, tau=0.2)
#' #Compute inverse (approximately)
#' fast_inv(FM, phi=1.5, tau=0.2)
#'
#' @export
fast_process <- function(R, phi=NULL, tau=NULL, method='fast', control=list()){
  if(method == 'fast'){
    eig <- eigen(R)
    out <- list(Q=eig$vectors, r=eig$values, R=R)
    if(length(phi) > 0){
      out$phi <- phi
    }
    if(length(tau) > 0){
      out$tau <- tau
    }
    class(out) <- 'FM_fast'
    return(out)
  }
  if(method == 'approx'){
    eig <- eigen(R)
    Q <- eig$vectors
    N <- nrow(Q)
    N0 <- which(eig$values < control$tol)[1]
    if(is.na(N0)) stop("There are no eigenvalues smaller than control$tol.")
    Z <- array(0, dim=c(N, N))
    ind <- (N0+1):N
    for(i in 1:N){
      for(j in i:N){
        Z[i,j] <- Z[j,i] <- sum(Q[i,ind]*Q[j,ind])
      }
    }
    out <- list(Q=Q, r=eig$values, tau=tau, R=R, Z=Z, N0=N0, N=N, tol=control$tol)
    if(length(phi) > 0){
      out$phi <- phi
    }
    if(length(tau) > 0){
      out$tau <- tau
    }
    class(out) <- 'FM_approx'
    return(out)
  }
  if(method == 'ts'){
    eig <- eigen(R)
    Q <- eig$vectors
    r <- eig$values
    N <- nrow(R)
    C <- array(0, dim=c(N, N, control$P+1))
    a <- control$a
    for(i in 1:N){
      for(j in i:N){
        for(p in 0:control$P){
          C[i,j,p+1] <- C[j,i,p+1] <- sum(Q[i,]*Q[j,]*r^p/(r*a+tau)^(p+1))
        }
      }
    }
    out <- list(Q=Q, r=r, R=R, P=control$P, a=a, C=C, tau=tau)
    if(length(phi) > 0){
      out$phi <- phi
    }
    if(length(tau) > 0){
      out$tau <- tau
    }
    class(out) <- 'FM_ts'
    return(out)
  }
  stop('method not recognized')
}


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
#' R <- corr_matrix(seq(0,1,length.out=n), beta=0.2)
#' FM <- fast_process(R, method='approx', control=list(tol=1e-9))
#' #Compute determinant (approximately)
#' fast_det(FM, phi=1.5, tau=0.2)
#' #Compute inverse (approximately)
#' fast_inv(FM, phi=1.5, tau=0.2)
#' @export
fast_det <- function(FM, phi=NULL, tau=NULL, log=TRUE){
  phi <- ifelse(length(phi) > 0, phi, FM$phi)
  tau <- ifelse(length(tau) > 0, tau, FM$tau)
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
#' @param ... appropriate arguments for the corresponding class.
#' @return Returns the inverse of the covariance matrix
#' @export
fast_inv.FM_fast <- function(FM, phi=NULL, tau=NULL, ...){
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
#' @param ... appropriate arguments for the corresponding class.
#' @return Returns an approximation of the inverse of the covariance matrix
#' @examples
#' n <- 100
#' R <- corr_matrix(seq(0,1,length.out=n), beta=0.2)
#' FM <- fast_process(R, method='approx', control=list(tol=1e-3))
#' #Compute determinant (approximately)
#' fast_det(FM, phi=1.5, tau=0.2)
#' #Compute inverse (approximately)
#' fast_inv(FM, phi=1.5, tau=0.2)
#' @export
fast_inv.FM_approx <- function(FM, phi=NULL, tau=NULL, ...){
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
#' @param ... appropriate arguments for the corresponding class.
#' @return Returns an approximation of the inverse of the covariance matrix
#' @examples
#' n <- 100
#' R <- corr_matrix(seq(0,1,length.out=n), beta=0.2)
#' FM <- fast_process(R, method='approx', control=list(tol=1e-9))
#' #Compute determinant (approximately)
#' fast_det(FM, phi=1.5, tau=0.2)
#' #Compute inverse (approximately)
#' fast_inv(FM, phi=1.5, tau=0.2)
#' @export
fast_inv.FM_ts <- function(FM, phi=NULL, ...){
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

#' Fast Matrix Algebra for BMC
#'
#' See fast_inv.class for details.
#'
#' @param FM an object of class "FM_fast", "FM_approx" or "FM_ts"
#' @param ... appropriate arguments for the corresponding class.
#' @return Returns (possibly an approximation of) the inverse of the covariance matrix
#'
#' @export
fast_inv <- function(FM, ...){
  UseMethod("fast_inv")
}


#' Correlation Matrix
#'
#' Genersates the power exponential or Matern correlation matrix for a set of n design points in d-dimensional space.
#' Default is Gaussian (power=2). See GPfit::corr_matrix for details.
#'
#' @param X nxd matrix of design points
#' @param beta correlation parameter. If 1/kappa is length scale then beta=log10(kappa).
#' @param corr list specifying the correlation function. See GPfit package for details.
#' @examples
#' R <- corr_matrix(seq(0, 1, length.out=10), beta=log10(1.0))
#'
#' @export
corr_matrix <- function (X, beta, corr = list(type = "exponential", power = 2))
{
  if (!is.matrix(X)) {
    X = as.matrix(X)
  }
  d = ncol(X)
  n = nrow(X)
  if (d != length(beta)) {
    stop("The dimensions of beta and X do not match. \n")
  }
  R <- matrix(0, nrow = n, ncol = n)
  rcoord <- cbind(rep(seq_len(n - 1L), times = rev(seq_len(n -
                                                             1L))), unlist(lapply(X = rev(seq_len(n - 1L)), FUN = function(nn,
                                                                                                                           nm) seq_len(nn) + nm - nn, nm = n)))
  absdiff <- abs(X[rcoord[, 2L], ] - X[rcoord[, 1L], ])
  switch(corr$type, exponential = {
    power <- corr$power
    if (is.null(power)) {
      power <- 1.95
    }
    absdiff <- absdiff^power
    Beta <- matrix(beta, nrow = (length(absdiff)/d), ncol = d,
                   byrow = TRUE)
    Rtemp <- (10^Beta) * absdiff
    Rtemp <- rowSums(Rtemp)
    R[rcoord] <- Rtemp
    R <- R + t(R)
    R <- exp(-R)
    R
  }, matern = {
    nu <- corr$nu
    if (is.null(nu)) {
      nu <- 2.5
    }
    Beta <- 10^beta
    Beta <- matrix(Beta, ncol = d, nrow = (length(absdiff)/d),
                   byrow = TRUE)
    absdiff <- 2 * sqrt(nu) * absdiff * (Beta)
    pos <- which(absdiff == 0)
    Rtemp <- 1/(gamma(nu) * 2^(nu - 1)) * (absdiff)^nu *
      besselK(x = absdiff, nu = nu)
    Rtemp[pos] <- 1
    Rtemp <- apply(X = Rtemp, MARGIN = 1L, FUN = prod)
    R[rcoord] <- Rtemp
    R <- R + t(R)
    diag(R) <- 1
    R
  }, stop("corr type must be 'exponential' or 'matern'"))
}


