#' Localized Ensemble of Approximate Gaussian Processes
#'
#' This function is a modification of the LA-GP framework of Gramacy and Apley
#' designed for cases where parallel predictions are not possible (i.e. MCMC).
#' The leapGP offers a quadratic training algorithm which leads to fast predictions.
#'
#' @param X a matrix of training locations (1 row for each training instance)
#' @param Y a vector of training responses (length(y) == nrow(X))
#' @param H the number of prediction hubs desired. Defaults to ceiling(sqrt(length(Y))).
#' @param scale logical. Do we want the scale parameter to be returned for predictions? If TRUE,
#'            the matrix K^{-1} will be stored for each hub.
#' @param iso boolean. Is correlation function isotropic? (Current version ignores this argument)
#' @param n local neighborhood size
#' @param start number of starting points for neighborhood (between 6 and n inclusive)
#' @param frac logical. If TRUE, information is returned about the fraction of training points which are in their own prediction neighborhoods.
#' @param verbose  logical. Deault is FALSE
#' @param justdoit logical. Force leapGP to run using specified parameters (may be incredibly time consuming).
#' @param ... optional arguments to be passed to laGP()
#' @return a univariate prediction and an updated list of hubs. Also returns scale parameter if scale=TRUE
#' @examples
#' Xnew <- matrix(runif(100), nrow=50, ncol=2)
#' X <- matrix(runif(100), nrow=50, ncol=2)
#' Y <- apply(X, 1, prod)
#' preds1 <- leapGP(Xnew, X, Y)
#'
#' #Or equivalently
#' leap <- leapGP_build(X, Y)
#' preds2 <- rep(NA, 100)
#' for(i in 1:100){
#'    preds2[i] <- leapGP_predict(leap, Xnew[i,])
#' }
#'
#' #Or used with slapGP
#' leap <- leapGP_synch(leap, rho=0.95)
#' hubs <- leap$hubs
#' preds3 <- rep(NA, 100)
#' for(i in 1:100){
#'    emulator <- slapGP(Xnew[i,], X, Y, hubs=hubs)
#'    preds3[m] <- emulator$pred
#'    hubs <- emulator$hubs
#' }
#' @export
leapGP_build <- function(X, Y, H=NA,
                         scale=F, iso=T, n=NA, start=NA, frac=TRUE,
                         verbose=FALSE, justdoit=FALSE, ...){
  #Get parameters
  N <- length(Y)
  if(is.na(n)){
    n <- ceiling(sqrt(N))
  }
  if(n > 300 | justdoit){
    stop('n is very large - consider choosing a smaller value. Set justdoit=TRUE to override.')
  }

  if(is.na(H)){
    H <- ceiling(sqrt(N))
  }
  if(is.na(start)){
    start <- floor(max(6, 0.1*n))
  }

  if(H < N){
    if(N < 5000 | justdoit){
      #Use PAM to find medoids
      coord_ind <- cluster::pam(X, H, pamonce=5)$id.med
    }else{
      warning('N is very large, a subset of 5000 points was used. Set justdoit=TRUE to override.')
      rand_ind <- sample(N,5001)
      coord_ind <- rand_ind[cluster::pam(X[rand_ind,], H, pamonce=5)$id.med]
    }
  }else{
    coord_ind <- 1:N
  }

  #Start building hubs
  hubs <- list()
  for(hh in 1:H){
    h <- laGP::laGP(Xref=matrix(X[coord_ind[hh], ], ncol=ncol(X)),
              X=X, Z=Y,
              start=start, end=n, ...)
    new_hub <- list(coord_id=coord_ind[hh], neigh=h$Xi, kappa=1/h$mle$d,
                    nugget=1/h$mle$g)
    #Calculate psi (called a in code)
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
    hubs[[hh]] <- new_hub
    if(verbose){
      if((hh %% ((H - H%%10)/10)) == 0){
        cat('a dime towards a dollar\n')
      }
    }
  }
  out <- list(hubs=hubs, X=X, Y=Y)
  class(out) <- 'leapGP'
  #Calculate frac_J if needed
  if(frac){
    Jset <- hubs[[1]]$neigh
    for(h in 2:length(hubs)){
      Jset <- union(Jset, hubs[[h]]$neigh)
    }
    frac <- length(Jset)/N
    out$frac_J <- frac
  }
  return(out)
}


#' Localized Ensemble of Approximate Gaussian Processes
#'
#' This function synchronizes an object of class "leapGP", so that is also an object of "slapGP"
#'
#' @param leapGP an object of class "leapGP"
#' @param rho parameter controlling time-accuracy tradeoff (default = 0.95)
#' @return an object of class "leapGP" AND of class "slapGP"
#' @examples
#' Xnew <- matrix(runif(100), nrow=50, ncol=2)
#' X <- matrix(runif(100), nrow=50, ncol=2)
#' Y <- apply(X, 1, prod)
#' preds1 <- leapGP(Xnew, X, Y)
#'
#' #Or equivalently
#' leap <- leapGP_build(X, Y)
#' preds2 <- rep(NA, 100)
#' for(i in 1:100){
#'    preds2[i] <- leapGP_predict(leap, Xnew[i,])
#' }
#'
#' #Or used with slapGP
#' leap <- leapGP_synch(leap, rho=0.95)
#' hubs <- leap$hubs
#' preds3 <- rep(NA, 100)
#' for(i in 1:100){
#'    emulator <- slapGP(Xnew[i,], X, Y, hubs=hubs)
#'    preds3[m] <- emulator$pred
#'    hubs <- emulator$hubs
#' }
#' @export
leapGP_synch <- function(leapGP, rho=0.95){
  hubs <- leapGP$hubs; X <- leapGP$X
  for(h in 1:length(hubs)){
    hubs[[h]]$epsilon <- sqrt(log(1/rho)/hubs[[h]]$kappa)
    hubs[[h]]$coord <- X[hubs[[h]]$coord_id, ]
  }
  leapGP$hubs <- hubs
  class(leapGP) <- c('leapGP', 'slapGP')
  return(leapGP)
}



#' Localized Ensemble of Approximate Gaussian Processes
#'
#' A function for prediction with an object of class leapGP.
#'
#' @param leapGP an object of class "leapGP"
#' @param Xnew the location at which prediction is requested
#' @param scale logical. Do we want the scale parameter to be returned for predictions? If TRUE,
#'            the matrix K^{-1} will be stored for each hub.
#' @param iso logical. Is correlation function isotropic? (Currently not supported)
#' @return an object of class "leapGP" AND of class "slapGP"
#' @examples
#' Xnew <- matrix(runif(100), nrow=50, ncol=2)
#' X <- matrix(runif(100), nrow=50, ncol=2)
#' Y <- apply(X, 1, prod)
#' preds1 <- leapGP(Xnew, X, Y)
#'
#' #Or equivalently
#' leap <- leapGP_build(X, Y)
#' preds2 <- rep(NA, 100)
#' for(i in 1:100){
#'    preds2[i] <- leapGP_predict(leap, Xnew[i,])
#' }
#'
#' #Or used with slapGP
#' leap <- leapGP_synch(leap, rho=0.95)
#' hubs <- leap$hubs
#' preds3 <- rep(NA, 100)
#' for(i in 1:100){
#'    emulator <- slapGP(Xnew[i,], X, Y, hubs=hubs)
#'    preds3[m] <- emulator$pred
#'    hubs <- emulator$hubs
#' }
#' @export
leapGP_predict <- function(leapGP, Xnew,
                           scale=F, iso=T, ...){
  #Get parameters
  X <- leapGP$X; Y <- leapGP$Y; hubs <- leapGP$hubs
  N <- length(Y); n <- length(hubs[[1]]$neigh)
  if(ncol(X) != length(Xnew)){
    stop('wrong dimensions for Xnew')
  }
  #Assume single prediction requested
  H <- X[unlist(lapply(hubs, function(h) h$coord_id)),]
  #Find nearest hub using kd-tree
  NH <- RANN::nn2(H, matrix(Xnew, nrow=1), k=1)
  curr_hub <- hubs[[NH$nn.idx]]
  #Make prediction and return hubs unmodified
  kvec <- rep(NA, n)
  for(i in 1:n){
    ii <- curr_hub$neigh[i]
    kvec[i] <- exp(-curr_hub$kappa*sum((X[ii,]-Xnew)^2))
  }
  pred <- as.numeric(matrix(kvec, nrow=1)%*%curr_hub$a)
  if(scale){
    scale <- curr_hub$phi/n*(1 - matrix(kvec, nrow=1)%*%curr_hub$Kinv%*%matrix(kvec, ncol=1))
    out <- c(pred, scale)
    names(out) <- c('pred', 'scale')
    return(c(pred, scale))
  }
  return(pred)
}



#' Localized Ensemble of Approximate Gaussian Processes
#'
#' A wrapper for combining training and prediction for leapGP
#'
#' @param X a matrix of training locations (1 row for each training instance)
#' @param Y a vector of training responses (length(y) == nrow(X))
#' @param H the number of prediction hubs desired. Defaults to ceiling(sqrt(length(Y))).
#' @param scale logical. Do we want the scale parameter to be returned for predictions? If TRUE,
#'            the matrix K^{-1} will be stored for each hub.
#' @param iso boolean. Is correlation function isotropic? (Current version ignores this argument)
#' @param n local neighborhood size
#' @param frac logical. If TRUE, information is returned about the fraction of training points which are in their own prediction neighborhoods.
#' @param verbose  logical. Deault is FALSE
#' @param justdoit logical. Force leapGP to run using specified parameters (may be incredibly time consuming).
#' @param ... optional arguments to be passed to laGP()
#' @return a univariate prediction and an updated list of hubs. Also returns scale parameter if scale=TRUE
#' @examples
#' Xnew <- matrix(runif(100), nrow=50, ncol=2)
#' X <- matrix(runif(100), nrow=50, ncol=2)
#' Y <- apply(X, 1, prod)
#' preds1 <- leapGP(Xnew, X, Y)
#'
#' #Or equivalently
#' leap <- leapGP_build(X, Y)
#' preds2 <- rep(NA, 100)
#' for(i in 1:100){
#'    preds2[i] <- leapGP_predict(leap, Xnew[i,])
#' }
#'
#' #Or used with slapGP
#' leap <- leapGP_synch(leap, rho=0.95)
#' hubs <- leap$hubs
#' preds3 <- rep(NA, 100)
#' for(i in 1:100){
#'    emulator <- slapGP(Xnew[i,], X, Y, hubs=hubs)
#'    preds3[m] <- emulator$pred
#'    hubs <- emulator$hubs
#' }
#' @export
leapGP <- function(Xnew, X, Y, H=NA,
                   scale=F, n=NA, frac=FALSE, start=NA,
                   verbose=TRUE, justdoit=FALSE, ...){
  if(verbose){ cat('Training...\n') }
  leap <- leapGP_build(X, Y, H=H, scale=scale, n=n, start=start, frac=frac, ...)
  if(verbose){ cat('completed. Starting predictions.') }
  if(!is.matrix(Xnew)){
    res <- leapGP_predict(leap, Xnew, scale=scale, ...)
    return(res)
  }
  M <- nrow(Xnew) #Number of predictions requested
  preds <- list()
  for(m in 1:M){
    preds[[m]] <- leapGP_predict(leap, Xnew[m,], scale=scale, ...)
  }
  if(scale){
    return(preds)
  }
  return(unlist(preds))
}




