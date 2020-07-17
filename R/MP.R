#' Density of the Moment Penalization Prior
#'
#' Returns the unnormalized density of the MP prior with parameters w1 and w2
#'
#' @param x a vector of length p
#' @param w1 normalized penalty associated with first moment. Default is 1
#' @param w1 normalized penalty associated with second moment. Default is 1
#' @param log logical. Should density be returned on a log scale?
#' @param normalized logical. Should density be normalized (default is FALSE)
#' @param ... additional parameters passed to get_constMP (if norm=TRUE)
#' @return returns the density of the MP(w1, w2) prior
#' @examples
#' x <- rnorm(10)
#' dMP(x, w1=5, w2=2)
#' X <- matrix(rnorm(10*30), nrow=30)
#' apply(X, 1, dMP, w1=5, w2=2)
#' @export
dMP <- function(x, w1=1, w2=1, log=FALSE, normalized=FALSE, ...){
  p <- length(x); M <- mean(x); V <- var(x)
  ld <- (-p*w1/2*M^2 - (p-1)*w2/4*(V-1)^2) + log(min(w1,w2) > 0)
  if(normalized){
    ld <- ld + log(get_constMP(p, w1, w2, ...))
  }
  ifelse(log, return(ld), return(exp(ld)))
}

#' Normalizing Constant for the Moment Penalization Prior
#'
#' Returns the normalizing constant for the Moment Penalization Prior with parameters p, w1 and w2
#' @param p the dimension of the vector having an MP prior
#' @param w1 normalized penalty associated with first moment. Default is 1
#' @param w1 normalized penalty associated with second moment. Default is 1
#' @param MC number of MC iterations. Default is 1e7, but larger values are recommended when speed is not important.
#' @return returns the normalizing constant of the MP prior
#' @examples get_constMP(7, 1, 3, MC=1e9)
#' @export
get_constMP <- function(p, w1=1, w2=1, MC=1e7){
  w_1 <- max(w1, w2); w_2 <- min(w1, w2)
  a <- 2-sqrt(3); b <- 2+sqrt(3)
  U <- rchisq(MC, p-1)
  res <- (2*base::pi/w1)^(p/2)*mean(exp(-(p-1)*w1/4*(w2*(U/((p-1)*w1)-1)^2-2/(p-1)*U)))
  return(1/res)
}


#' Random Sampling from the Moment Penalization Prior
#'
#' Uses an optimized rejection sampling framework to sample from the MP prior
#'
#' @param n number of samples requested
#' @param p dimension of random vector
#' @param w1 normalized penalty associated with first moment. Default is 1
#' @param w1 normalized penalty associated with second moment. Default is 1
#' @param acceptance logical. If TRUE, a list is returned with an attribute giving the acceptance rate.
#' @param ... additional parameters passed to get_constMP (if norm=TRUE)
#' @return returns the density of the MP(w1, w2) prior
#' @examples X <- rMP(100, 7, 1, 3)
#' @export
rMP <- function(n, p, w1, w2, acceptance=FALSE){
  res <- matrix(NA, nrow=n, ncol=p)
  s <- optim(max(1, w2)/2,
             fn=function(s) s^(-p/2)*exp((p-1)*s*(2*w2+s)/(4*w2)),
             lower=0, upper=max(1, w2), method='Brent')$par
  C <- log((2*pi/s)^(p/2)*exp((p-1)*s*(2*w2+s)/(4*w2))) + 1e-3
  cnt <- 1
  cnt2 <- 0
  while(cnt <= n){
    cand <- rnorm(p, 0, sqrt(1/s))
    acc <- dMP(cand, w1, w2, log=T) - C - sum(dnorm(cand, 0, sqrt(1/s), log=T))
    if(log(runif(1)) < acc){
      res[cnt,] <- cand
      cnt <- cnt + 1
    }
    cnt2 <- cnt2 + 1
  }
  if(acceptance){
    res <- list(samples=res, acceptance=n/cnt2)
  }
  return(res)
}


#' Probability of Prior Coherency
#'
#' Computes the probability of prior coherency using Monte Carlo
#'
#' @param M a vector of posterior M values (mean of a parameter set)
#' @param V a vector of posterior V values (variance of a parameter set)
#' @param p dimension of the parameter set
#' @param type indicates what should be returned. See details
#' @param MC number of Monte Carlo samples.
#' @details type = 1 leads to PPC calculation based on posterior mean of M and V (as described in paper)
#'
#' type = 2 is the posterior mean of PPC
#'
#' type = 3 returns posterior samples of PPC
#'
#' type = 4 returns additional information (used for plotting)
#' @examples
#' gamma <- cbind(rnorm(1000, 1.2, 0.2),
#'                rnorm(1000, -0.4, 0.3),
#'                rnorm(1000, -0.6, 0.15))
#' M <- apply(gamma, 1, mean)
#' V <- apply(gamma, 1, var)
#' ppc(M, V, 3)
#' diagnose_ppc(M, V, p, control=list(levels=c(0.5, 0.9, 0.99), fill_col='purple'))
#' @export
ppc <- function(M, V, p, type=1, MC=1e4){
  if(length(M) != length(V)){
    stop('M and V must have same length')
  }

  #Perform Monte Carlo Simulations
  Mi <- rnorm(MC, 0, 1/sqrt(p))
  Vi <- rchisq(MC, p-1)/(p-1)
  Di <- dnorm(Mi, 0, 1/sqrt(p))*dchisq((p-1)*Vi, p-1)

  if(type == 1){
    D0 <- dnorm(mean(M), 0, 1/sqrt(p))*dchisq((p-1)*mean(V), p-1)
    return(mean(Di < D0))
  }

  if(type == 2){
    D0 <- dnorm(M, 0, 1/sqrt(p))*dchisq((p-1)*V, p-1)
    return(mean(vapply(D0, function(z, Di) mean(Di < z), Di=Di, FUN.VALUE=1)))
  }

  if(type >= 3){
    D0 <- dnorm(M, 0, 1/sqrt(p))*dchisq((p-1)*V, p-1)
    ppc <- vapply(D0, function(z, Di) mean(Di < z), Di=Di, FUN.VALUE=1)
    if(type == 3){ return(ppc) }
    return(list(ppc=ppc, Di=Di))
  }
}

#' Diagnostic Plot for Probability of Prior Coherency
#'
#' Diagnostic plot for detecting overfitting using probability of prior coherency
#'
#' @param M a vector of posterior M values (mean of a parameter set)
#' @param V a vector of posterior V values (variance of a parameter set)
#' @param p dimension of the parameter set
#' @param MC number of Monte Carlo samples.
#' @param control a list of control parameters for plotting. See details.
#' @param ... additional parameters passed to plot() function.
#' @details control is a named list consisting of any of the following.
#'
#' levels - the desired contour levels for plotting (default c(0.5, 0.9, 0.95, 0.99))
#'
#' h - step-size for plotting contours (default 0.005)
#'
#' tol - tolerance for plotting contours (default 1e-7)
#'
#' prob - a probability bound for setting xlim and ylim (default 0.9999)
#'
#' fill_col - the color of the contours regions (default 'dodgerblue')
#'
#' point_col - the color of the points (default 'orange')
#'
#' cex - 0.2 the size of the points representing PPC posterior samples (default 0.5)
#'
#' leg_cex - the size of the legend (default 1.4)
#'
#' @examples X <- rMP(100, 7, 1, 3)
#' @export
diagnose_ppc <- function(M, V, p, MC=1e4, control=NULL, ...){
  #Process control arguments
  default_control <- list(h=.005, tol=1e-7, prob=.9999, levels=c(.5, .9, .95, .99), fill_col='dodgerblue', point_col='orange', cex=0.5, leg_cex=1.4)
  if(is.null(control)){
    control <- default_control
  }else{
    if(class(control) != 'list'){
      stop('Control must be a named list')
    }
    not_provided <- setdiff(names(default_control), names(control))
    control[not_provided] <- default_control[not_provided]
  }

  ppc <- ppc(M, V, p, type=4, MC=MC)

  #Get contours
  xbounds <- c(min(mean(M), -qnorm(control$prob, 0, 1/sqrt(p))), max(mean(M), qnorm(control$prob, 0, 1/sqrt(p))))
  ybounds <- c(0, max(mean(V), qchisq(control$prob, p-1)/(p-1)))
  f <- function(x, y, p){
    dnorm(x, 0, 1/sqrt(p))*dchisq((p-1)*y, p-1)
  }
  q <- quantile(ppc$Di, 1-control$levels)
  contours <- list()
  for(i in 1:length(control$levels)){
    q_curr <- q[i]
    x_points <- y_points <- NULL
    for(yy in seq(ybounds[1], ybounds[2], by=control$h)){
      f_curr <- function(x){
        (f(x, yy, p) - q_curr)^2 - control$tol
      }
      roots <- rootSolve::uniroot.all(f_curr, lower=xbounds[1], upper=xbounds[2])
      if(length(roots) > 0){
        x_points <- c(x_points, roots)
        y_points <- c(y_points, rep(yy, length(roots)))
      }
    }
    contours[[i]] <- list(x=x_points, y=y_points)
  }

  #Make Plot
  plot(NULL, xlim=xbounds, ylim=ybounds, xlab='Mean', ylab='Variance', ...)
  for(i in 1:length(control$levels)){
    temp <- contours[[i]]
    ord <- chull(temp$x, temp$y)
    polygon(temp$x[ord], temp$y[ord], border='lightgray', lty=3, lwd=1, col=adjustcolor(control$fill_col, alpha.f=1/(length(control$levels)+1)))
  }
  nl <- length(control$levels)
  myfill <- rep(NA, nl)
  for(i in 1:nl){
    myfill[i] <- adjustcolor(control$fill_col, alpha.f= (nl+1-i)/(nl+1))
  }
  legend('topleft', c(expression(paste('p'[pc], ' levels')), 1 - control$levels),  border=c(NA, rep('lightgray', nl)), bty='n',
         fill=c(NA, myfill), cex=control$leg_cex)
  points(M, V, pch=16, cex=control$cex, col=adjustcolor(control$point_col, alpha.f=0.5))
  points(mean(M), mean(V), pch=21, cex=2, bg=control$point_col)
  points(0, (p-2)/(p-1), pch='+', cex=1.5)
  text(0.65, 3, expression(paste('p'[pc],' =')), cex=2)
  text(.99, 3.1, round(mean(ppc$ppc), 4), cex=2)
}

#' Random Sampling from Z-Regularization Prior
#'
#' Draws samples from Z-regularization prior. Alternative to MP prior with w1=w2=Inf (when sigma_R = 0).
#'
#' @param n number of samples requested
#' @param p dimension of random vector
#' @param sigma_R relaxation parameter (default 0)
#' @return returns the density of the Z regularization prior
#' @examples X <- rZreg(100, 7)
#' @export
rZreg <- function(n, p, sigma_R=0){
  Z <- matrix(rnorm(n*p), nrow=n)
  Zeta <- rnorm(n, 0, sigma_R)
  Zout <- matrix(NA, nrow=n, ncol=p)
  for(i in 1:n){
    Mz <- sum(Z[i,] + Zeta[i])/p
    Sz <- sqrt(sum((Z[i,] - Mz)^2)/(p-1-sigma_R/p))
    Zout[i,] <- (Z[i,] - Mz)/Sz
  }
  return(Zout)
}













