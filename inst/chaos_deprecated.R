
bayes_chaos_old <- function(X, y, max_degree_init=2, max_order_init=1, max_degree=6, lambda=1, rho=0.02, verbose=TRUE){
  n <- nrow(X)
  p <- ncol(X)
  md <- max_degree_init
  mo <- max_order_init

  mu_y <- mean(y)
  sig_y <- sd(y)
  y <- (y - mu_y)/sig_y

  res <- bayes_chaos_wrapper_old(X, y, n, p, md, mo, mu_y, sig_y, max_degree, lambda, rho, verbose)
  return(res)
}

bayes_chaos_wrapper_old <- function(X, y, n, p, md, mo, mu_y, sig_y, realmd, lambda, rho, verbose){
  # Create a list of p sequences from 1 to n
  if(verbose){
    cat("Starting model with max degree = ", md, " and max order = ", mo, "\n", sep="")
  }
  if(verbose) cat("\tComputing initial phi matrix\n")
  A_set <- generate_A(p, md, mo)
  A_deg <- apply(A_set, 1, sum)
  A_ord <- apply(A_set, 1, function(aa) sum(aa > 0))
  N_alpha <- nrow(A_set)
  if(verbose) cat("\tFound ", N_alpha, " combinations.\n", sep="")
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

  # Get partial correlation coefficients
  if(verbose) cat("\tComputing partial correlation coefficients\n")
  rho_vec <- rep(rho, N_alpha)
  if(N_alpha > n || N_alpha > 2000){
    # Get preliminaries
    ZZ <- cbind(y, phi)
    p_alpha <- N_alpha + 1
    R <- cor(ZZ)
    rho_max <- max(abs(R-diag(rep(1, N_alpha + 1))))
    S_full <- R * tcrossprod(apply(ZZ, 2, sd))
    J_star <- get_J_star(c(n, N_alpha))

    if(J_star < p_alpha){
      if(verbose) cat("\t\t Estimating the full precision matrix with glasso\n")
      P_full <- glasso::glasso(S_full, rho_max*rho)$wi
      P_full <- (P_full + t(P_full))/2
      sparsity <- (sum(P_full == 0)-N_alpha)/(N_alpha^2 - N_alpha)
      if(verbose){
        cat("\t\t ", round(100*sparsity, 3), "% of the off-diagonal entries are zero\n", sep="")
        cat("\t\t Computing inverse of D matrix (Schur complement)")
      }

      # Get Schur Complement
      A_curr <- P_full[1:J_star, 1:J_star, drop=FALSE]
      B_curr <- P_full[1:J_star, (J_star+1):p_alpha, drop=FALSE]
      #C_curr <- P_full[(J_star+1):p_alpha, 1:J_star, drop=FALSE] # t(B_curr)
      D_curr <- P_full[(J_star+1):p_alpha, (J_star+1):p_alpha, drop=FALSE]
      Di_curr <- solve(D_curr)
    }
    for(j in 1:N_alpha){
      jj <- j+1
      if(jj <= J_star){
        # Get with glasso directly
        P_j <- glasso::glasso(S_full[1:jj, 1:jj], rho_max*rho_vec[j])$wi
      }else{
        if(j == N_alpha){
          P_j <- P_full
          rm(P_full)
        }else{
          # Partition-inverse equations (Barnett 1979)
          vv <- 1/Di_curr[1,1]
          gg <- Di_curr[2:(p_alpha - jj + 1), 1,drop=FALSE]
          Di_curr <- Di_curr[2:(p_alpha - jj + 1), 2:(p_alpha - jj + 1), drop=FALSE] +
            tcrossprod(gg)*vv
          A_curr <- cbind(rbind(A_curr, P_full[jj,1:(jj-1),drop=FALSE]), P_full[1:jj, jj,drop=FALSE])
          B_curr <- rbind(B_curr[,-1,drop=FALSE], P_full[jj,(jj+1):p_alpha,drop=FALSE])

          # Schur complement
          P_j <- A_curr - B_curr %*% Di_curr %*% t(B_curr)
        }
      }
      diag(P_j) <- pmax(diag(P_j), 1e-9)
      if(is.nan(sqrt(P_j[1,1]*P_j[jj,jj]))) browser()
      rr[j] <- -P_j[1,jj]/sqrt(P_j[1,1]*P_j[jj,jj])
    }
    # End big N_alpha case
  }else{
    if(verbose) cat("\t\tFitting linear models: 0/", N_alpha, ", ", sep="")
    for(i in 2:N_alpha){
      if(N_alpha > 20 && ((i %% round(N_alpha/5)) == 0)){
        if(verbose) cat(i, "/", N_alpha, ", ",sep="")
      }
      eps_y <- lm(y ~ phi[,1:(i-1)])$residuals
      eps_p <- lm(phi[,i] ~ phi[,1:(i-1)])$residuals
      rr[i] <- cor(eps_y, eps_p)
    }
  }
  if(any(abs(rr) > 1)){
    rr <- rr/max(abs(rr))
  }
  ord <- rev(order(rr^2))
  A_set <- A_set[ord,]
  A_deg <- A_deg[ord]
  A_ord <- A_ord[ord]
  phi <- phi[,ord]
  rr <- rr[ord]

  # Get KIC for various models
  if(verbose) cat("\n\tRanking models based on KIC\n")
  plot(rr^2)
  # Add coefficient column in
  phi <- cbind(rep(1, n), phi)
  A_set <- rbind(rep(0, p), A_set)
  K_trunc <- max(which(rr^2 >= 0.005)) + 1
  if(verbose) cat("\t\t Throwing out bases with low r^2. ", K_trunc, " remain.\n", sep="")
  abline(v=K_trunc)
  browser()

  KIC <- rep(NA, K_trunc)
  KIC[1] <- Inf
  Caa <- diag(c((md+mo-1)*mo^2, (A_deg + A_ord -1)*A_ord^2))
  Caa_inv <- diag(1/diag(Caa))
  best <- list(k=0, KIC=Inf, map=NULL)
  for(k in 2:K_trunc){
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
    if(max(A_deg[1:(best$k-1)]) + 1 >= md){
      res <- bayes_chaos_wrapper_old(X, y, n, p, md+2, mo+1, mu_y, sig_y, realmd, lambda, rho, verbose)
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

get_J_star <- function(np){
  n <- np[1]
  p <- np[2]
  f <- function(x, n, p){
    (1/4)*x^4 + (n/3-5/6)*x^3 + (n/2+4*p - 1/4)*x^2 + ((n-1)/6 + p*(1-4*p))*x
  }

  coefficients <- c((n-1)/6 + p*(1-4*p),
                    (n + 8*p - 1/2),
                    (n-5/2),
                    1)
  roots <- Re(polyroot(coefficients))
  roots <- roots[which(roots > 0)]
  roots <- roots[which(roots < p)]
  J_star <- max(2, min(p, ceiling(roots[which.min(f(roots, n, p))])))
  return(J_star)
}

# This function was found by
#   1. Fitting LASSO model to J_star ~ (n,p) to select coefficients
#   2. Fitting lm() model with selected coefficients
#   3. Rounding to nice numbers and verification that it works reasonably well.
get_J_star_hueristic <- function(np){
  n <- np[1]
  p <- np[2]
  #res <- 25.6 + n*.0091 + p*.1028 - 7.8*log(n) + 3.7*log(p) # From LASSO
  #res <- 43.5 - 11.4*log(n) - n/257 + 7*log(p) + p/16       # From lm
  res <- 50 - 11*log(n) - n/250 + 7*log(p) + p/16
  min(max(1, ceiling(res)), p)
}
