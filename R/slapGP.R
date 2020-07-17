#' Sequence of Local Approximate Gaussian Processes
#'
#' This function is a modification of the LA-GP framework of Gramacy and Apley
#' designed for cases where parallel predictions are not possible (i.e. MCMC).
#' The slapGP framework offers users a time-accuracy tradeoff based on the rho parameter.
#'
#' @param Xnew the location at which prediction is requested
#' @param X a matrix of training locations (1 row for each training instance)
#' @param Y a vector of training responses (length(y) == nrow(X))
#' @param rho parameter controlling time-accuracy tradeoff (default = 0.95)
#' @param hubs a list of current prediction hubs
#' @param scale logical. Do we want the scale parameter to be returned for predictions? If TRUE,
#'            the matrix K^{-1} will be stored for each hub.
#' @param iso logical. Is correlation function isotropic? (Currently not supported)
#' @param n local neighborhood size
#' @param start number of starting points for neighborhood (between 6 and n inclusive)
#' @param ... optional arguments to be passed to laGP()
#' @return a univariate prediction and an updated list of hubs. Also returns scale parameter if scale=TRUE
#' @examples
#'
#' hubs <- list()
#' X <- matrix(runif(100), nrow=50, ncol=2)
#' Y <- apply(X, 1, prod)
#' preds <- rep(NA, 1000)
#' for(i in 1:1000){
#'    Xnew <- runif(2)
#'    emulator <- slapGP(Xnew, X, Y, hubs=hubs)
#'    hubs <- emulator$hubs
#'    preds[i] <- emulator$pred
#' }
#' @export
slapGP <- function(Xnew, X, Y, rho=0.95, hubs=list(),
                   scale=F, iso=TRUE, n=NA, start=NA, ...){
  #Get parameters
  N <- length(Y)
  if(is.na(n)){
    n <- ceiling(sqrt(N))
  }
  if(is.na(start)){
    start <- floor(max(6, 0.1*n))
  }

  #Find nearest hub (if possible)
  if(length(hubs) > 0){
    #Matrix of hub locations
    H <- matrix(unlist(lapply(hubs, function(h) h$coord)), nrow=length(hubs), byrow=TRUE)
    #Find nearest hub using kd-tree
    NH <- RANN::nn2(H, matrix(Xnew, nrow=1), k=1)
    curr_hub <- hubs[[NH$nn.idx]]
    #Is hub close enough to make a prediction?
    if(NH$nn.dists < curr_hub$epsilon){
      #Make prediction and return hubs unmodified
      kvec <- rep(NA, n)
      for(i in 1:n){
        ii <- curr_hub$neigh[i]
        kvec[i] <- exp(-curr_hub$kappa*sum((X[ii,]-Xnew)^2))
      }
      #Prepare object for return
      out <- list()
      out$pred <- matrix(kvec, nrow=1)%*%curr_hub$a
      out$pred_hub <- NH$nn.idx
      if(scale){
        out$scale <- curr_hub$phi/n*(1 - matrix(kvec, nrow=1)%*%curr_hub$Kinv%*%matrix(kvec, ncol=1))
      }
      out$hubs <- hubs
      return(out)
    }
  }
  #Make a new hub (and use it for prediction)
  h <- laGP::laGP(Xref=matrix(Xnew, ncol=ncol(X)), X=X, Z=Y,
            start=start, end=n, ...)
  #Prepare object for return
  out <- list(pred=h$mean, scale=h$s2)

  #Make a new hub list
  new_hub <- list(coord=Xnew, neigh=h$Xi, kappa=1/h$mle$d,
                  nugget=1/h$mle$g)
  new_hub$epsilon <- sqrt(-log(rho)/new_hub$kappa)
  K <- diag(rep(1, n))
  for(i in 2:n){
    for(j in 1:(n-1)){
      ii <- new_hub$neigh[i]
      jj <- new_hub$neigh[j]
      K[i,j] <- K[j, i] <- exp(-new_hub$kappa*sum((X[ii,]-X[jj,])^2))
    }
  }
  Kinv <- solve(K+diag(rep(1e-9, n)))
  new_hub$a <- Kinv%*%matrix(Y[new_hub$neigh], nrow=n)
  if(scale){
    new_hub$Kinv <- Kinv
    new_hub$phi <- matrix(Y[new_hub$neigh], ncol=n)%*%(new_hub$a)
  }
  #Add new hub to hubs
  i <- length(hubs) + 1
  hubs[[i]] <- new_hub
  out$hubs <- hubs
  return(out)
}


