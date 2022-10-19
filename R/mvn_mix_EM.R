#' EM Algorithm to fit Finite Mixtures of Multivariate Normal Distributions
#'
#' Returns a list of mixture distribution for various number of mixture components
#'
#' @param X data, a nxp matrix
#' @param K a vector (or scalar for a single fit) of
#' @param kappa Regularization parameter. Each Sigma is shrank towards Sx at iteration t with weight 1/t^kappa
#' @param epsilon stopping criterion for EM algorithm
#' @param max_iter maximum number of iterations of EM algorithm
#' @param verbose should information on status be printed
#' @param print_iter how often should information being printed? Ignored when verbose=FALSE.
#' @param crit a subset of c("aic", "bic", "hqic") corresponding to Aikike's, Bayesian and Hanaan-Quinn information criteria.
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
mvn_mix <- function(X, K=1:3, kappa=0.75, epsilon=1e-5, max_iter=2000, verbose=FALSE, print_iter = 10, crit=NULL,
                    minibatch=nrow(X)){
  if(length(epsilon) == 1){
    epsilon <- rep(epsilon, length(K))
  }
  if(length(epsilon) != length(K)){
    stop("epsilon should be either scalar or a vector of length = length(K)")
  }
  if(minibatch < nrow(X)){
    if(max(epsilon) > 0){
      warning("Stochastic EM algorithm with minibatching ignores epsilon. Set epsilon = 0 to ignore this message")
      epsilon <- rep(0, length(K))
    }
  }

  out <- list()
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
      sigma[[j]] <- rWishart(1, df=ncol(X)+2, SX)[,,1]
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
      for(i in 1:n){
        pi_ij[i,] <- pi_ij[i,]/sum(pi_ij[i,])
      }

      # M-STEP
      Sj <- list()
      for(j in 1:p){
        tmp <- matrix(0, nrow=d, ncol=d)
        for(i in ind){
          xi <- X[i,]
          tmp <- tmp + pi_ij[i,j]*tcrossprod(xi - mu[[j]])
        }
        Sj[[j]] <- tmp
      }
      #Update component probabilities
      pi <- apply(pi_ij, 2, mean)

      #Update covariance matrices
      for(j in 1:p){
        sigma[[j]] <- (2*a_n(n)*SX + Sj[[j]])/(2*a_n(n) + n*pi[j])
      }

      #Update mean vectors
      for(j in 1:p){
        tmp <- 0
        for(i in ind){
          xi <- X[i,]
          tmp <- tmp + pi_ij[i,j]*xi
        }
        mu[[j]] <- tmp/(n*pi[j])
      }

      #Compute log likelihood
      log_lik[iter] <- 0
      for(i in ind){
        tmp <- 0
        xi <- X[i,]
        for(j in 1:p){
          tmp <- tmp + pi[j]*dmvnorm(xi, mu[[j]], sigma[[j]])
        }
        log_lik[iter] <- log_lik[iter] + log(tmp) + log(n) - log(minibatch)
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

    out[[k]] <- list(log_lik = log_lik,
                     pi = pi,
                     mu = mu,
                     sigma = sigma)
    if(verbose){
      cat("\n\n")
    }
  }

  if(!is.null(crit)){
    nvar <- ncol(X)
    Kmax <- length(out)
    bic <- ll <- rep(NA, Kmax)
    aic <- hqic <- rep(NA, Kmax)
    k <- 1
    bic[1] <- -2*sum(log(dmvnorm(X, apply(X, 2, mean), cov(X)))) +
      (k-1 + nvar*k + nvar*(nvar+1)/2*k)*log(nrow(X))
    aic[1] <- -2*sum(log(dmvnorm(X, apply(X, 2, mean), cov(X)))) +
      (k-1 + nvar*k + nvar*(nvar+1)/2*k)*2
    hqic[1] <- -2*sum(log(dmvnorm(X, apply(X, 2, mean), cov(X)))) +
      (k-1 + nvar*k + nvar*(nvar+1)/2*k)*2*log(log(n))
    ll[1] <- sum(log(dmvnorm(X, apply(X, 2, mean), cov(X))))


    for(k in 2:Kmax){
      tmp <- out[[k]]$log_lik
      bic[k] <- -2*max(tmp) + (k-1 + nvar*k + nvar*(nvar+1)/2*k)*log(nrow(X))
      aic[k] <- -2*max(tmp) + (k-1 + nvar*k + nvar*(nvar+1)/2*k)*2
      hqic[k]<- -2*max(tmp) + (k-1 + nvar*k + nvar*(nvar+1)/2*k)*2*log(log(nrow(X)))
      ll[k] <- max(tmp)
    }

    out2 <- list(fits=out)
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



