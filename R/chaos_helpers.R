legendre_poly <- function(x, j){
  if(j == 0){
    return(rep(1, length(x)))
  }
  if(j == 1){
    return(x)
  }
  n <- j - 1
  res <- ((2*n+1)*x*legendre_poly(x, n) - n*legendre_poly(x, n-1))/(n+1)
  return(res)
}

ss_legendre_poly <- function(x, j){
  sqrt(2*j+1)*legendre_poly(2*x-1, j)
}

generate_A <- function(p, d, q) {
  # Initialize the result matrix
  res <- NULL

  # Recursive function
  generate <- function(current_set, current_sum, current_nonzero, index) {
    if (current_sum > d || current_nonzero > q) return()
    if (index > p) {
      if(current_nonzero > 0){
        res <<- rbind(res, current_set)
      }
      return()
    }

    for (value in 0:(d - current_sum)) {
      new_set <- current_set
      new_set[index] <- value
      generate(new_set,
               current_sum + value,
               current_nonzero + as.numeric(value > 0),
               index + 1)
    }
  }

  # Start the generation process
  generate(integer(p), 0, 0, 1)
  return(res)
}


estimate_map <- function(y, phi, Cinv, sig0=1, tol=1e-3, max_iter=1000){
  n <- length(y)
  p <- ncol(phi)
  sig_curr <- sig0
  Cinv <- Cinv
  flag <- TRUE
  cnt <- 1
  while(flag){
    Cinv_curr <-  crossprod(phi)/sig_curr^2 + Cinv
    tmp <- tryCatch({
      solve(Cinv_curr)
    }, error=function(e){ NULL })
    if(is.null(tmp)){
      Cinv_curr <- as.matrix(Matrix::nearPD(Cinv_curr)$mat)
    }
    C_curr <- solve(Cinv_curr)
    a_curr <- (tcrossprod(C_curr, phi)%*%y)/sig_curr^2
    yhat   <- phi%*%a_curr
    sig_new <- sqrt(mean((y-yhat)^2))
    delta <- abs(sig_new-sig_curr)
    sig_curr <- sig_new

    if(delta < tol){
      flag <- FALSE
    }
    if(cnt >= max_iter){
      flag <- FALSE
    }
    cnt <- cnt + 1
  }
  out <- list(a=a_curr, sig=sig_curr, Caa_inv=Cinv_curr)
  return(out)
}
