#' KL Divergence estimated via KNN
#'
#' Returns the KLD estimator of Wang et al. (2009). Useful in high dimensions. For univariate samples, use \code{kld()} instead.
#'
#' @param x an n by p matrix of samples from P
#' @param y an m by p matrix of samples from Q
#' @param automate logical; If TRUE, eps is set using the data-driven approach of Eq 26-29 (using k=1, by default).
#' @param k number of nearest neighbors to use.
#' @param eps neighborhood radius to use (only k or epsilon should be specified). Uses generalized estimator (Eq 17 and Eq 25).
#' @param whiten logical; Should x and y be whitened using Eq 32-33 first?
#' @param type See details
#' @param lambda parameter (between 0 and 1) for the Population Stability Index (ignored unless type = 4)
#' @param n_mc number of monte carlo samples used when type = 3 or 4
#' @details Follows Wang et al. 2009.
#' @examples
#' x <- matrix(rnorm(2*500), ncol=2)
#' x[,2] <- 0.6*x[,2] + rnorm(2*500, 0, 0.4)
#' y <- matrix(rnorm(2*450, 0, 0.9), ncol=2)
#'
#' kld_knn(x, y, TRUE, k=3)
#' kld_knn(x, y, TRUE, k=3, type=2)
#' kld_knn(x, y, TRUE, k=3, whiten=TRUE, type=2)
#' kld_knn(x, y, FALSE, k=5)
#' kld_knn(x, y, FALSE, eps = 1.5)
#'
#' @export
kld_knn <- function(x, y, automate=TRUE, k=NULL, eps=NULL, whiten=FALSE, type=1, lambda=0.5, n_mc=1){
  if(type == 1){
    return(compute_kld_knn(x, y, automate, k, eps, whiten))
  }
  if(type == 2){
    D1 <- compute_kld_knn(x, y, automate, k, eps, whiten)
    D2 <- compute_kld_knn(y, x, automate, k, eps, whiten)
    return((D1+D2)/2)
  }
  if(type == 3){
    if(lambda != 1/2){
      warning("lambda must = 1/2 when type = 3. Set type = 4 if you want to adjust lambda.")
    }
    lambda <- 1/2
  }
  if(type >= 3){
    nx <- nrow(x)
    ny <- nrow(y)
    p <- ncol(x)
    if(ncol(y) != p){
      stop("x and y must have same number of columns")
    }
    N <- min(nx, ny)
    D <- rep(NA, n_mc)
    z <- matrix(NA, nrow=N, ncol=p)
    for(i in 1:n_mc){
      ind <- sample(seq(1, nx+ny), size=N,
                   replace=FALSE,
                   prob=c(rep(lambda/nx, nx),
                          rep((1-lambda)/ny, ny)))
      z0 <- rbind(x, y)[ind,]
      # I don't know why this is needed,
      # but estimator gets weird without it.
      # I guess it's because the samples have
      # to be independent, and the way i'm doing
      # the mixture the samples are dependent
      # (e.g., z dependent on x)
      z <- z0 + rnorm(N, 0, 2*diff(range(z0))/N)
      D1  <- compute_kld_knn(x, z, automate, k, eps, whiten)
      D2  <- compute_kld_knn(y, z, automate, k, eps, whiten)
      lambda*D1 + (1-lambda)*D2
      D[i] <- lambda*D1 + (1-lambda)*D2
    }
    return(D)
  }
  stop("Not a valid type")
  return(FALSE)
}


compute_kld_knn <- function(x, y, automate, k, eps, whiten){
  n <- nrow(x)
  m <- nrow(y)
  p <- ncol(x)

  # Equation 30-33
  if(whiten){
    C <- cov(rbind(x, y))
    L <- solve(chol(C))
    mu <- colMeans(rbind(x, y))

    x <- (x - rep(mu, each=n))%*%L
    y <- (y - rep(mu, each=m))%*%L
  }

  if(automate){
    method = 26 # Equation 26
  }else{
    if(!is.null(k)){
      method = 5 # Equation 5
    }else{
      if(!is.null(eps)){
        method = 25 # Equation 25
      }else{
        stop("When automate is FALSE, k or eps must be specified")
      }
    }
  }

  # Fixed k
  if(method == 5){
    nn_x <- RANN::nn2(x, x, k=k+1) # + 1 to exclude self
    nn_y <- RANN::nn2(y, x, k=k)

    rho <- nn_x$nn.dists[,k+1]
    nu  <- nn_y$nn.dists[,k]

    res <- p * sum(log(nu) - log(rho)) / n + log(m) - log(n-1)
    return(res)
  }

  # Fixed epsilon
  if(method == 25){
    nn_x <- RANN::nn2(x, x, searchtype="radius", radius=eps, k=n)
    nn_y <- RANN::nn2(y, x, searchtype="radius", radius=eps, k=m)

    li <- apply(nn_x$nn.idx, 1,
                function(jj){
                  first_zero <- which(jj == 0)[1]
                  if(is.na(first_zero)){
                    res <- length(jj)
                  }else{
                    res <- first_zero - 1 - 1 # Extra -1 to exclude self
                  }
                })

    ki <- apply(nn_y$nn.idx, 1,
                 function(jj){
                   first_zero <- which(jj == 0)[1]
                   if(is.na(first_zero)){
                     res <- length(jj)
                   }else{
                     res <- first_zero - 1
                   }
                 })

    if(any(c(li, ki) == 0)){
      warning("no neighbors for at least one point inside radius. Expact erratic estimates.")
      term <- mean(digamma(li + 0.001) - digamma(ki + 0.001))
    }else{
      term <- mean(digamma(li) - digamma(ki))
    }
    res <- term + log(m) - log(n-1)
    return(res)
  }

  if(method == 26){
    if(is.null(k)) k <- 1
    nn_x <- RANN::nn2(x, x, k=k+1) # + 1 to exclude self
    nn_y <- RANN::nn2(y, x, k=k)
    # Get epsilon
    eps_vec <- pmax(nn_x$nn.dists[,k+1], nn_y$nn.dists[,k])

    # Make another pass using eps_vec as radius
    li <- ki <- rho <- nu <- rep(NA, n)
    for(i in 1:n){
      xx <- matrix(x[i,], nrow=1)

      # Do x vs x
      nn_x_v2 <- nn2(x, matrix(xx, nrow=1), searchtype="radius", radius=eps_vec[i] * 1.001, k=n)
      jj <- nn_x_v2$nn.idx
      first_zero <- which(jj == 0)[1]
      li[i] <- first_zero - 1 - 1 # extra -1 to exclude self
      rho[i] <- nn_x_v2$nn.dists[,first_zero-1]

      # Do x vs y
      nn_y_v2 <- nn2(y, matrix(xx, nrow=1), searchtype="radius", radius=eps_vec[i] * 1.001, k=m)
      jj <- nn_y_v2$nn.idx
      first_zero <- which(jj == 0)[1]
      ki[i] <- first_zero - 1
      nu[i] <- nn_y_v2$nn.dists[,first_zero-1]
    }

    # Compute estimator eq 29
    res <- p * sum(log(nu) - log(rho)) / n
    res <- res + mean(digamma(li) - digamma(ki))
    res <- res + log(m) - log(n-1)
    return(res)
  }
}
