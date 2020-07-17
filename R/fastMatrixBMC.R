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
#' @return returns an object of class "FM_<method>" containing the information needed for fast computation of matrix inverse and determinant
#' @examples
#' n <- 100
#' R <- GPfit::corr_matrix(seq(0,1,length.out=n), beta=0.2)
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




