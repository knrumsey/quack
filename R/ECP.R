#' Emulation of the Conditional Posterior (General)
#'
#' This function models the conditional posterior distribution of alpha conditional on a
#' set of nuisance parameters gamma. Draws from the modularization posterior can be obtained
#' by pairing this function with ECP_sample. This function is designed to give maximal flexibility.
#' Consider using ECP_norm for a more straightforward implementation and ECP_multi for multiple alpha
#'
#' @param lpost The log posterior. First argument must be a vector of parameters for which
#' full Bayesian inference is desired. Second argument is a vector of nuisance parameters
#' for which modularization is desired. Second argument must be named "gamma".
#' @param init a set of initial values for alpha. Can also be a function which generates starting values
#' based on gamma
#' @param L the budget.
#' @param MCMC a named list of arguments for the Metro_Hastings function.
#' @param f_alpha an optional function for transforming the physical parameters. Helpful for reducing
#' to a single dimension.
#' @param estimate_psi function whose argument is the result of f_alpha(al), where al is a
#' posterior sample from the conditional posterior. Should return a vector of r parameter estimates.
#' Defaults to univariate normal.
#' @param gamma_generate can be either (i) an L x q matrix of nuisance parameters, (ii) a function
#' (with no arguments) that generates a single instance of gamma (length q) from the prior or (iii)
#' a named list with two components (mu and sigma) which are q-length vectors of prior mean/sd's.
#' Please ensure that class(gamma_generate) is either 'matrix', 'function' or 'list'.
#' @param ... additional arguments passed to lpost
#' @examples
#' set.seed(42)
#' a <- 0; g <- 0.1
#' n1 <- 10; n2 <- 90
#' y <- rnorm(n1+n2, a + c(rep(0, n1), rep(g, n2)), 1)
#' lpost <- function(a, gamma, y, n1, n2){
#'   res <- sum(dnorm(y, a + c(rep(0, n1), rep(gamma, n2)), 1, log=TRUE))
#'   res <- res + dnorm(a, 0, 1, log=TRUE) + dnorm(gamma, 0, 1, log=TRUE)
#'   return(res)
#' }
#' ecp <- ECP(lpost, L=30, init=0, gamma_generate=list(mu=0, sigma=1), y=y, n1=n1, n2=n2)
#' hist(ECP_sample(ecp, M=1000, gamma_generate=list(mu=0, sigma=1)), breaks=30)
#' ECP_profile_plot(ecp)
#' @export
ECP <- function(lpost, L=30, init=NULL,
                MCMC=list(pars=NULL, iterations=5000, burn_in=1000, thin=1,
                          prop_sigma=NULL, adapt_par=c(100,20,0.5,0.75)),
                f_alpha=function(x) x[1], estimate_psi=function(x) c(mean(x), sd(x)),
                gamma_generate=NULL, verbose=TRUE, ...){
  #Process gamma_generate
  if(is.function(gamma_generate)){
    gamma0 <- gamma_generate()
    q <- length(gamma0)
    gamma <- matrix(NA, nrow=L, ncol=q)
    gamma[1,] <- gamma0
    for(ell in 2:L){
      gamma[ell,] <- gamma_generate()
    }
  }
  if(is.matrix(gamma_generate)){
    gamma <- gamma_generate
    q <- ncol(gamma)
    if(nrow(gamma) != L){
      warning('Ignroing L: Budget L differs from nrow(gamma)')
    }
  }
  if(is.list(gamma_generate)){
    q <- length(gamma_generate[[1]])
    if(L < 1000){
      gamma <- lhs::maximinLHS(L, q)
    }else{
      gamma <- lhs::randomLHS(L, q)
    }
    for(i in 1:q){
      gamma[,i] <- qnorm(gamma[,i], gamma_generate$mu[i], gamma_generate$sigma[i])
    }
  }

  print(gamma)
  for(ell in 1:L){
    if(is.function(init)){
      alpha0 <- init(gamma[ell,])
    }else{
      alpha0 <- init
    }
    chains <- MHadaptive::Metro_Hastings(lpost, pars=alpha0, prop_sigma=MCMC$prop_sigma,
                             iterations = MCMC$iterations, burn_in = MCMC$burn_in,
                             adapt_par = MCMC$adapt_par, quiet=TRUE, gamma=gamma[ell,],
                             ...)$trace
    chains <- matrix(chains, nrow=MCMC$iterations+1-MCMC$burn_in)
    chains <- matrix(chains[seq(1, nrow(chains), by=MCMC$thin),], nrow=MCMC$iterations+1-MCMC$burn_in)
    alpha <- apply(chains, 1, f_alpha)
    psi_curr <- estimate_psi(alpha)
    if(!exists("tmp_PSI")){
      r <- length(psi_curr)
      tmp_PSI <- matrix(NA, nrow=L, ncol=r)
    }
    tmp_PSI[ell,] <- psi_curr
    if(verbose){
      if((ell %% ceiling(L/10)) == 0){
        print('10% of Budget Spent')
      }
    }
  }

  GP_list <- mlegp::mlegp(X=gamma, Z=tmp_PSI, verbose=0)
  rm(tmp_PSI)
  return(GP_list)
}



#' Sampling from the Modularization Posterior with ECP
#'
#' This function produces draws from the Modularization Posterior of alpha (with respect to gamma).
#'
#' @param ecp An object returned from the ECP() function.
#' @param M number of samples requested
#' @param gamma_generate can be either (i) an L x q matrix of nuisance parameters, (ii) a function
#' (with no arguments) that generates a single instance of gamma (length q) from the prior or (iii)
#' a named list with two components (mu and sigma) which are q-length vectors of prior mean/sd's.
#' Please ensure that class(gamma_generate) is either 'matrix', 'function' or 'list'.
#' @param alpha_generate a function to generate random samples of alpha conditional on psi. Should
#' correspond to the estimate_psi function used in ECP() function. Takes one argument: psi and returns
#' a single alpha draw. Defaults to univariate normal.
#' @examples
#' set.seed(42)
#' a <- 0; g <- 0.1
#' n1 <- 10; n2 <- 90
#' y <- rnorm(n1+n2, a + c(rep(0, n1), rep(g, n2)), 1)
#' lpost <- function(a, gamma, y, n1, n2){
#'   res <- sum(dnorm(y, a + c(rep(0, n1), rep(gamma, n2)), 1, log=TRUE))
#'   res <- res + dnorm(a, 0, 1, log=TRUE) + dnorm(gamma, 0, 1, log=TRUE)
#'   return(res)
#' }
#' ecp <- ECP(lpost, L=30, init=0, gamma_generate=list(mu=0, sigma=1), y=y, n1=n1, n2=n2)
#' hist(ECP_sample(ecp, M=1000, gamma_generate=list(mu=0, sigma=1)), breaks=30)
#' ECP_profile_plot(ecp)
#' @export
ECP_sample <- function(ecp, M=1, gamma_generate=NULL, alpha_generate=NULL){
  if(is.null(alpha_generate)){
    alpha_generate <- function(psi){rnorm(1, psi[1], psi[2])}
  }
  #Process gamma_generate
  if(is.function(gamma_generate)){
    gamma0 <- gamma_generate()
    q <- length(gamma0)
    gamma <- matrix(NA, nrow=M, ncol=q)
    gamma[1,] <- gamma0
    for(m in 2:M){
      gamma[m,] <- gamma_generate()
    }
  }
  if(is.matrix(gamma_generate)){
    gamma <- gamma_generate
    q <- ncol(gamma)
    if(nrow(gamma) != M){
      warning('Ignoring M: Input M differs from nrow(gamma)')
    }
  }
  if(is.list(gamma_generate)){
    q <- length(gamma_generate[[1]])
    if(M < 1000){
      gamma <- lhs::maximinLHS(M, q)
    }else{
      gamma <- lhs::randomLHS(M, q)
    }
    for(i in 1:q){
      gamma[,i] <- qnorm(gamma[,i], gamma_generate$mu[i], gamma_generate$sigma[i])
    }
  }

  r <- ecp$numGPs
  out <- matrix(NA, nrow=M, ncol=ecp$numDim)
  for(m in 1:M){
    psi_curr <- rep(NA, r)
    for(j in 1:r){
      psi_curr[j] <- as.numeric(predict(ecp[[j]], gamma[m,]))
    }
    out[m,] <- alpha_generate(psi_curr)
  }
  return(out)
}

#' Profile Plot for ECP Algorithm
#'
#' Plots quantiles of the posterior for alpha conditional on gamma.
#'
#' @param ecp an object generated from the ecp function (requires normality)
#' @param bounds defines the bounds of gamma (from, to, by).
#' Should be a q x 3 matrix when there are multiple gamma.
#' @param fix when q > 1, what values should the remaining modularization parameters be fixed to?
#' @param quantiles the quantiles of the conditional distribution to plot
#' @param ... additional arguments passed to plot() and lines()
#' @examples
#' set.seed(42)
#' a <- 0; g <- 0.1
#' n1 <- 10; n2 <- 90
#' y <- rnorm(n1+n2, a + c(rep(0, n1), rep(g, n2)), 1)
#' lpost <- function(a, gamma, y, n1, n2){
#'   res <- sum(dnorm(y, a + c(rep(0, n1), rep(gamma, n2)), 1, log=TRUE))
#'   res <- res + dnorm(a, 0, 1, log=TRUE) + dnorm(gamma, 0, 1, log=TRUE)
#'   return(res)
#' }
#' ecp <- ECP(lpost, L=30, init=0, gamma_generate=list(mu=0, sigma=1), y=y, n1=n1, n2=n2)
#' hist(ECP_sample(ecp, M=1000, gamma_generate=list(mu=0, sigma=1)), breaks=30)
#' ECP_profile_plot(ecp)
#' @export
ECP_profile_plot <- function(ecp, bounds=c(0,1,.01), fix=NULL, quantiles=c(0.95, 0.5, 0.05), ...){
  q <- length(ecp$params)
  r <- ecp$numGPs
  if(r != 2) stop('only normal distribution is currently supported')
  if(q==1){
    gamma_vec <- seq(bounds[1], bounds[2], by=bounds[3])
    alpha_vec <- matrix(NA, nrow=length(gamma_vec), ncol=length(quantiles))
    for(i in 1:length(gamma_vec)){
      psi_curr <- rep(NA, r)
      for(j in 1:r){
        psi_curr[j] <- as.numeric(predict(ecp[[j]], gamma_vec[i]))
      }
      alpha_vec[i,] <- qnorm(quantiles, psi_curr[1], psi_curr[2])
    }
    plot(NULL, xlim=bounds[1:2], ylim=range(alpha_vec), xlab='gamma', ylab='alpha', ...)
    for(i in 1:length(quantiles)){
      lines(gamma_vec, alpha_vec[,i], ...)
    }
  }
  if(q > 1){
    if(length(fix) != q){stop('fix must be specified when using multiple gamma')}
    for(cnt in 1:q){
      gamma_vec <- seq(bounds[cnt,1], bounds[cnt,2], by=bounds[cnt,3])
      alpha_vec <- matrix(NA, nrow=length(gamma_vec), ncol=length(quantiles))
      for(i in 1:length(gamma_vec)){
        gamma_curr <- fix
        gamma_curr[cnt] <- gamma_vec[i]
        psi_curr <- rep(NA, r)
        for(j in 1:r){
          psi_curr[j] <- as.numeric(predict(ecp[[j]], gamma_curr))
        }
        alpha_vec[i,] <- qnorm(quantiles, psi_curr[1], psi_curr[2])
      }
      plot(NULL, xlim=bounds[1:2], ylim=range(alpha_vec), xlab=paste0('gamma', cnt), ylab='alpha', ...)
      for(i in 1:length(quantiles)){
        lines(gamma_vec, alpha_vec[,i], ...)
      }
    }
  }
}





