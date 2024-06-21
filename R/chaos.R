#' Bayesian Sparse Polynomial Chaos Expansion
#'
#' The emulation approach of Shao et al. (2017)
#'
#' @param X A dataframe or matrix of predictors scaled to be between 0 and 1
#' @param y a reponse vector
#' @param max_degree_init The maximum polynomial degree (for initial model)
#' @param max_order_init The maximum order of interaction (for initial model)
#' @param max_degree The maximum degree allowed (bounds the computation time).
#' @param lambda A parameter fed to the KIC calculations. Larger values will lead to sparser models.
#' @param verbose Logical. Should progress information be printed?
#' @details Implements the Bayesian sparse PCE method described by Shao et al. (2017).
#' @references Shao, Q., Younes, A., Fahs, M., & Mara, T. A. (2017). Bayesian sparse polynomial chaos expansion for global sensitivity analysis. Computer Methods in Applied Mechanics and Engineering, 318, 474-496.
#' @examples
#' X <- lhs::maximinLHS(500, 3)
#' y <- apply(X, 1, duqling::dms_simple) + rnorm(500, 0, 0.1)
#' fit <- bayes_chaos(X, y)
#' plot(fit)
#' @export
bayes_chaos <- function(X, y, max_degree_init=2, max_order_init=1, max_degree=20, lambda=1, verbose=TRUE){
  n <- nrow(X)
  p <- ncol(X)
  md <- max_degree_init
  mo <- max_order_init

  mu_y <- mean(y)
  sig_y <- sd(y)
  y <- (y - mu_y)/sig_y

  res <- bayes_chaos_wrapper(X, y, n, p, md, mo, mu_y, sig_y, max_degree, lambda, verbose)
  return(res)
}

bayes_chaos_wrapper <- function(X, y, n, p, md, mo, mu_y, sig_y, realmd, lambda, verbose){
  # Create a list of p sequences from 1 to n
  if(verbose){
    cat("Starting model with max degree = ", md, " and max order = ", mo, "\n", sep="")
  }
  if(verbose) cat("\tComputing initial phi matrix\n")
  A_set <- generate_A(p, md, mo)
  A_deg <- apply(A_set, 1, sum)
  A_ord <- apply(A_set, 1, function(aa) sum(aa > 0))
  N_alpha <- nrow(A_set)

  phi <- matrix(NA, nrow=n, ncol=N_alpha)
  rr <- rep(NA, N_alpha)
  for(i in 1:N_alpha){
    curr <- rep(1, n)
    for(j in 1:p){
      curr <- curr * ss_legendre_poly(X[,j], A_set[i,j])
    }
    phi[,i] <- curr
    rr[i] <- cor(curr, y)
  }
  ord <- rev(order(rr^2))
  A_set <- A_set[ord,]
  A_deg <- A_deg[ord]
  A_ord <- A_ord[ord]
  phi <- phi[,ord]
  rr <- rr[ord]

  # Get conditional covariances
  if(verbose) cat("\tComputing partial correlation coefficients\n")
  for(i in 2:N_alpha){
    if(N_alpha > 50 && ((i %% round(N_alpha/10)) == 0)){
      if(verbose) cat("\t\t completed ", i, " out of ", N_alpha, ",   ", "\n",sep="")
    }
    eps_y <- lm(y ~ phi[,1:(i-1)])$residuals
    eps_p <- lm(phi[,i] ~ phi[,1:(i-1)])$residuals
    rr[i] <- cor(eps_y, eps_p)
  }
  ord <- rev(order(rr^2))
  A_set <- A_set[ord,]
  A_deg <- A_deg[ord]
  A_ord <- A_ord[ord]
  phi <- phi[,ord]
  rr <- rr[ord]

  # Get KIC for various models
  if(verbose) cat("\tRanking models based on KIC\n")

  # Add coefficient column in
  phi <- cbind(rep(1, n), phi)
  A_set <- rbind(rep(0, p), A_set)

  KIC <- rep(NA, N_alpha+1)
  KIC[1] <- Inf
  Caa <- diag(c(1, (A_deg + A_ord -1)*A_ord^2))
  Caa_inv <- diag(1/diag(Caa))
  best <- list(k=0, KIC=Inf, map=NULL)
  for(k in 2:(N_alpha+1)){
    map_k <- estimate_map(y, phi[,1:k,drop=FALSE], Caa_inv[1:k, 1:k, drop=FALSE])
    yhat_k <- phi[,1:k]%*%map_k$a
    KIC[k] <- -2*sum(dnorm(y, yhat_k, map_k$sig, log=TRUE)) -
              2*sum(dnorm(map_k$a, 0, sqrt(diag(Caa)), log=TRUE)) -
              lambda*(k+1)*log(2*pi) + log(det(map_k$Caa_inv))
    if(KIC[k] <= best$KIC){
      best$k   <- k
      best$KIC <- KIC[k]
      best$map <- map_k
    }
  }
  if(md <= realmd - 2){
    if(max(A_deg[1:best$k]) + 1 >= md){
      res <- bayes_chaos_wrapper(X, y, n, p, md+2, mo+1, mu_y, sig_y, realmd, lambda, verbose)
      return(res)
    }
  }
  #else
  obj <- list(map=best$map, phi=phi[,1:best$k,drop=FALSE],
              vars=A_set[1:best$k,,drop=FALSE],
              mu_y=mu_y, sigma_y=sig_y, KIC=best$KIC, X=X, y=y*sig_y + mu_y)
  class(obj) <- "bayes_chaos"
  return(obj)
}

#' Predict Method for class bayes_chaos
#'
#' See \code{bayes_chaos()} for details.
#'
#' @param object An object returned by the \code{bayes_chaos()} function.
#' @param newdata A dataframe of the same dimension as the training data.
#' @param samples How many posterior samples should be taken at each test point? If 0 or FALSE, then the MAP estimate is returned.
#' @details Predict function for bayes_chaos object.
#' @references Shao, Q., Younes, A., Fahs, M., & Mara, T. A. (2017). Bayesian sparse polynomial chaos expansion for global sensitivity analysis. Computer Methods in Applied Mechanics and Engineering, 318, 474-496.
#' @examples
#' X <- lhs::maximinLHS(500, 3)
#' y <- apply(X, 1, duqling::dms_simple) + rnorm(500, 0, 0.1)
#' fit <- bayes_chaos(X, y)
#' pred <- predict(fit, X)
#'
#' @export
predict.bayes_chaos <- function(object, newdata=NULL, samples=FALSE){
  if(is.null(newdata)){
    newdata <- object$X
  }
  XX <- newdata
  n <- nrow(XX)
  p <- ncol(XX)
  N_alpha <- nrow(object$vars)
  phi <- matrix(NA, nrow=n, ncol=N_alpha)
  for(i in 1:N_alpha){
    curr <- rep(1, n)
    for(j in 1:p){
      curr <- curr * ss_legendre_poly(XX[,j], object$vars[i,j])
    }
    phi[,i] <- curr
  }

  if(samples){
    pred <- matrix(NA, nrow=samples, ncol=n)
    for(i in 1:samples){
      a_sigma <- sqrt(1/diag(object$map$Caa_inv))
      a_hat <- object$map$a
      coeff <- rnorm(a_hat, a_hat, a_sigma)
      noise <- sqrt(1/rgamma(1, (n+2)/2, scale=2/(n*object$map$sig^2)))
      y_hat <- phi%*%coeff + rnorm(n, 0, noise)
      pred[i,] <- y_hat
    }
  }else{
    pred <- phi%*%object$map$a
  }
  pred <- object$mu_y + object$sigma_y * pred
  return(pred)
}

#' Plot Method for class bayes_chaos
#'
#' See \code{bayes_chaos()} for details.
#'
#' @param x An object returned by the \code{bayes_chaos()} function.
#' @param ... additional arguments passed to \code{plot}
#' @details Plot function for bayes_chaos object.
#' @references Shao, Q., Younes, A., Fahs, M., & Mara, T. A. (2017). Bayesian sparse polynomial chaos expansion for global sensitivity analysis. Computer Methods in Applied Mechanics and Engineering, 318, 474-496.
#' @examples
#' X <- lhs::maximinLHS(500, 3)
#' y <- apply(X, 1, duqling::dms_simple) + rnorm(500, 0, 0.1)
#' fit <- bayes_chaos(X, y)
#' plot(fit)
#' @export
plot.bayes_chaos <- function(x, ...){
  pred <- predict(x, x$X, samples=1000)
  yhat <- colMeans(pred)
  plot(x$y, yhat, ...)
  abline(0, 1, lwd=2, col='orange')

  ci <- apply(pred, 2, function(yy) quantile(yy, c(0.025, 0.975)))
  for(i in 1:ncol(ci)){
    segments(x$y[i], ci[1,i], x$y[i], ci[2,i])
  }
}






