#' EM Algorithm to fit Finite Mixtures of Multivariate Normal Distributions
#'
#' Returns a list of mixture distribution for various number of mixture components
#'
#' @param X data, a nxp matrix
#' @param K a vector (or scalar for a single fit) of
#' @param kappa Regularization parameter. Each Sigma is shrank towards Sx weight 1/nrow(X)^kappa
#' @param epsilon stopping criterion for EM algorithm
#' @param max_iter maximum number of iterations of EM algorithm
#' @param verbose should information on status be printed
#' @param print_iter how often should information being printed? Ignored when verbose=FALSE.
#' @param crit a subset of c("aic", "bic", "hqic") corresponding to Aikike's, Bayesian and Hanaan-Quinn information criteria.
#' @param minibatch WARNING: This doesn't seem to work right now. For stochastic udpating. epsilon should be set to -Inf.
#' @return a list with components corresponding to the fit using the number of components specified by `K`
#' @details A regularized approach based loosely on Chi & Lange (2014).
#'
#' Chi, Eric C., and Kenneth Lange. "Stable estimation of a covariance matrix guided by nuclear norm penalties." Computational statistics & data analysis 80 (2014): 117-128.
#' @examples
#' n1 <- 300
#' mu1 <- c(0, 0)
#' sigma1 <- matrix(c(1, 0.5, 0.5, 1), nrow=2, byrow=TRUE)
#' X1 <- mvtnorm::rmvnorm(n1, mu1, sigma1)
#'
#' n2 <- 500
#' mu2 <- c(1, 1.5)
#' sigma2 <- matrix(c(1, -0.7, -0.7, 1), nrow=2, byrow=TRUE)
#' X2 <- mvtnorm::rmvnorm(n2, mu2, sigma2)
#'
#' X <- rbind(X1, X2)
#' plot(X)
#'
#' #Fit mixture models for 1 to 3 compoenents
#' obj <- mvn_mix(X, K=1:3, epsilon=1e-3, verbose=TRUE, crit=c("aic", "bic", "hqic"))
#'
#' plot(1:3, obj$aic, type='l', col="orange")
#' lines(1:3, obj$bic, col="dodgerblue")
#' lines(1:3, obj$hqic, col="firebrick")
#' legend("bottomleft", c("aic", "bic", "hqic"), lwd=2, col=c("orange", "dodgerblue", "firebrick"))
#' @export
mvn_mix <- function(X, K=1:3, kappa=0.75, epsilon=1e-5, max_iter=2000, verbose=TRUE, print_iter = 10, crit=c("aic", "bic", "hqic"),
                    minibatch=nrow(X)){
  if(length(dim(X)) != 2){
    stop("X must be a matrix")
  }
  if(length(epsilon) == 1){
    epsilon <- rep(epsilon, length(K))
  }
  if(length(epsilon) != length(K)){
    stop("epsilon should be either scalar or a vector of length = length(K)")
  }
  if(minibatch < nrow(X)){
    if(max(epsilon) > -Inf){
      warning("Stochastic EM algorithm with minibatching ignores epsilon. Set epsilon = -Inf to ignore this message")
      epsilon <- rep(-Inf, length(K))
    }
    warning("Minibatch version of EM does not seem to work at this point. ")
  }

  out <- list()
  cnt <- 1
  for(k in K){
    if(verbose){
      cat("Number of components:", k)
    }
    N_comp <- k
    # Initialize model
    SX <- cov(X)
    a_n <- function(n, kap=0.75) n^(-kap)
    pi <- (tmp <- runif(N_comp))/sum(tmp)
    mu <- sigma <- list()
    for(j in 1:N_comp){
      mu[[j]]    <- X[sample(nrow(X), 1),]
      sigma[[j]] <- as.matrix(rWishart(1, df=ncol(X)+2, SX)[,,1])
    }
    n <- nrow(X)
    p <- N_comp
    d <- nrow(SX)
    pi_ij <- matrix(0, nrow=n, ncol=p)

    iter <- 1
    tol <- Inf
    log_lik <- rep(NA, max_iter)
    eps0 <- epsilon[which(K==k)]

    while(tol > eps0 & iter < max_iter){
      sub_samp <- sample(n, minibatch, replace=FALSE)
      #E-STEP
      for(i in sub_samp){
        for(j in 1:p){
          xi <- X[i,]
          pi_ij[i,j] <- pi[j] * dmvnorm(as.numeric(xi), as.numeric(mu[[j]]), sigma[[j]])
        }
      }
      for(i in sub_samp){
        pi_ij[i,] <- pi_ij[i,]/sum(pi_ij[i,])
      }

      # M-STEP
      Sj <- list()
      for(j in 1:p){
        tmp <- matrix(0, nrow=d, ncol=d)
        for(i in sub_samp){
          xi <- X[i,]
          tmp <- tmp + pi_ij[i,j]*tcrossprod(xi - mu[[j]])
        }
        Sj[[j]] <- tmp
      }
      #Update component probabilities
      pi <- apply(pi_ij[sub_samp,,drop=F], 2, mean)

      #Update covariance matrices
      for(j in 1:p){
        sigma[[j]] <- as.matrix((2*a_n(n, kappa)*SX + Sj[[j]])/(2*a_n(n, kappa) + n*pi[j]))
      }

      #Update mean vectors
      for(j in 1:p){
        tmp <- 0
        for(i in sub_samp){
          xi <- X[i,]
          tmp <- tmp + pi_ij[i,j]*xi
        }
        mu[[j]] <- tmp/(minibatch*pi[j])
      }

      #Compute log likelihood
      log_lik[iter] <- log(n) - log(minibatch)
      for(i in sub_samp){
        tmp <- 0
        xi <- X[i,]
        for(j in 1:p){
          tmp <- tmp + pi[j]*mvtnorm::dmvnorm(xi, mu[[j]], sigma[[j]])
        }
        log_lik[iter] <- log_lik[iter] + log(tmp)
      }
      if((iter %% print_iter) == 0){
        if(verbose)
          cat("\nIteration:", iter, ", tol:", tol, ", log-lik:", log_lik[iter])
      }
      if(iter > 1){
        tol <- log_lik[iter] - log_lik[iter - 1]
      }
      iter <- iter + 1
    }
    log_lik <- log_lik[1:(iter - 1)]

    out[[cnt]] <- list(k=k,
                       log_lik = log_lik,
                       pi = pi,
                       mu = mu,
                       sigma = sigma)
    cnt <- cnt + 1
    if(verbose){
      cat("\n\n")
    }
  }

  if(!is.null(crit)){
    nvar <- ncol(X)
    Kmax <- length(out)
    bic <- ll <- rep(NA, Kmax)
    aic <- hqic <- rep(NA, Kmax)
    kvec <- rep(NA, Kmax)

    kk <- 1
    for(k in K){
      kvec[kk] <- k
      if(k == 1){
        bic[kk] <- -2*sum(log(dmvnorm(X, apply(X, 2, mean), cov(X)))) +
          (k-1 + nvar*k + nvar*(nvar+1)/2*k)*log(nrow(X))
        aic[kk] <- -2*sum(log(dmvnorm(X, apply(X, 2, mean), cov(X)))) +
          (k-1 + nvar*k + nvar*(nvar+1)/2*k)*2
        hqic[kk] <- -2*sum(log(dmvnorm(X, apply(X, 2, mean), cov(X)))) +
          (k-1 + nvar*k + nvar*(nvar+1)/2*k)*2*log(log(n))
        ll[kk] <- sum(log(dmvnorm(X, apply(X, 2, mean), cov(X))))
      }else{
        tmp <- out[[kk]]$log_lik
        bic[kk] <- -2*max(tmp) + (k-1 + nvar*k + nvar*(nvar+1)/2*k)*log(nrow(X))
        aic[kk] <- -2*max(tmp) + (k-1 + nvar*k + nvar*(nvar+1)/2*k)*2
        hqic[kk]<- -2*max(tmp) + (k-1 + nvar*k + nvar*(nvar+1)/2*k)*2*log(log(nrow(X)))
        ll[kk] <- max(tmp)
      }
      kk <- kk + 1
    }

    out2 <- list(fits=out)
    out2$k <- kvec
    if("aic" %in% crit){
      out2$aic <- aic
    }
    if("bic" %in% crit){
      out2$bic <- bic
    }
    if("hqic" %in% crit){
      out2$hqic <- hqic
    }
    out2$ll <- ll
    return(out2)
  }
  return(out)
}

#' Generating samples from mvn_mix object
#'
#' @param n number of samples
#' @param obj an element of the list generated by mvn_mix()
#'
#' Chi, Eric C., and Kenneth Lange. "Stable estimation of a covariance matrix guided by nuclear norm penalties." Computational statistics & data analysis 80 (2014): 117-128.
#' @examples
#' x <- rgamma(100, 3, 1.2)
#' fit <- mvn_mix(x, K=1:3)
#' obj <- fit[[3]]
#' xx <- rmvn_mix(1000, obj)
#' hist(xx)
#' mvtnorm::dmvnorm(c(2, 3, 4), obj)
#' @export
rmvn_mix <- function(n, obj){
  p <- length(obj$mu[[1]])
  L <- length(obj$pi)
  indx <- sample(L, n, replace=TRUE, obj$pi)
  res <- matrix(NA, nrow=n, ncol=p)
  cnt <- 0
  for(i in 1:L){
    j <- which(indx == i)
    if(length(j) > 0){
      cnt <- max(cnt) + 1:length(j)
      res[cnt,] <- rmvnorm(length(j),
                         mean=obj$mu[[i]],
                         as.matrix(obj$sigma[[i]]))
    }
  }
  # Shuffle indices
  res <- res[sample(n,n,FALSE),]
  return(res)
}


#' Density of mvn_mix object
#'
#' @param x vector of locations to evaluate density
#' @param obj an element of the list generated by mvn_mix()
#'
#' Chi, Eric C., and Kenneth Lange. "Stable estimation of a covariance matrix guided by nuclear norm penalties." Computational statistics & data analysis 80 (2014): 117-128.
#' @examples
#' x <- rgamma(100, 3, 1.2)
#' fit <- mvn_mix(x, K=1:3)
#' obj <- fit[[3]]
#' xx <- rmvn_mix(1000, obj)
#' hist(xx)
#' dmvnorm(c(2, 3, 4), obj)
#' @export
dmvn_mix <- function(x, obj, log=FALSE){
  n <- length(x)
  p <- length(obj$mu[[1]])
  L <- length(obj$pi)
  res <- matrix(0, nrow=n, ncol=p)
  for(i in 1:n){
    for(j in 1:L){
      res[i] <- res[i] + obj$pi[j]*mvtnorm::dmvnorm(x[i],
                                           obj$mu[[j]],
                                           as.matrix(obj$sigma[[j]]))
    }
  }
  return(res)
}







