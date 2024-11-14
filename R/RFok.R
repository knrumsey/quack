#' Conformal Random Forest with Out-of-bag KNN
#'
#' A conformal RF based on Johansson et al. (2014). Rather than using a traditional calibration set, conformal intervals are constructed using out-of-bag samples for each tree.
#'
#' @param X a data frame or a matrix of predictors
#' @param y a response vector
#' @param k Number of nearest neighbors to use for out-of-bag error estimaters
#' @param beta Sensitivity parameter
#' @return An object with class "rfok"
#' @references Johansson, U., Boström, H., Löfström, T., & Linusson, H. (2014). Regression conformal prediction with random forests. Machine learning, 97, 155-176.
#' @examples
#'
#' X <- matrix(runif(150), nrow=50, ncol=3)
#' y <- apply(X, 1, duqling::ishigami)
#' fit <- rfok(X, y)
#'
#' Xnew <- matrix(runif(150), nrow=50, ncol=3)
#' predict(fit, Xnew)
#' @export
rfok <- function(X, y, k=5, beta=sd(y)/30, ...){
  n <- nrow(X)

  # Fit RF model
  fit <- randomForest::randomForest(X, y, keep.inbag=TRUE, ...)

  # Get predictions
  preds <- predict(fit, newdata=X, predict.all=TRUE)
  yhat <- preds$aggregate
  preds <- preds$individual

  # Get out of bag predictions
  out_of_bag_sets <- apply(fit$inbag, 1, function(xx) which(xx==0))
  yhat_oob <- rep(NA, n)
  for(i in 1:n){
    yhat_oob[i] <- mean(preds[i, out_of_bag_sets[[i]]])
  }

  # Get k nearest neighbors for each point
  neighbors <- RANN::nn2(X, X, k=min(k, nrow(X)))$nn.idx

  # Estimate the mus
  mu <- rep(NA, n)
  for(i in 1:n){
    ind <- neighbors[i,]
    mu[i] <- mean(abs(y[ind] - yhat_oob[ind]))
  }

  # Calculate non-conformity scores
  alpha <- abs(y - yhat_oob)/(mu + beta)
  #alpha_s <- quantile(alpha, conf)

  # Return object
  object <- list(fit=fit, k=k, beta=beta, alpha=alpha, oob_error=yhat_oob-y, X=X)
  class(object) <- "rfok"
  return(object)
}


#' Predict method for rfok
#'
#' A conformal RF based on Johansson et al. (2014). Rather than using a traditional calibration set, conformal intervals are constructed using out-of-bag samples for each tree.
#'
#' @param object Object returned by \code{rfok()}
#' @param newdata a data frame or matrix of new data
#' @param conf a vector of desired confidence intervals
#' @param samples Number of samples from the predictive distribution.
#' @return An object with class "rfok"
#' @examples
#'
#' X <- matrix(runif(150), nrow=50, ncol=3)
#' y <- apply(X, 1, duqling::ishigami)
#' fit <- rfok(X, y)
#'
#' Xnew <- matrix(runif(150), nrow=50, ncol=3)
#' predict(fit, Xnew)
#' @export
predict.rfok <- function(object, newdata, conf=0.95, samples=NULL, ...){
  pred <- predict(object$fit, newdata)
  n <- length(pred)

  if(is.null(conf) & is.null(samples)){
    return(pred)
  }

  # Get nearest neighbors
  neighbors <- RANN::nn2(object$X, newdata, k=min(object$k, nrow(object$X)))$nn.idx

  # Estimate mus
  mu <- rep(NA, n)
  for(i in 1:n){
    ind <- neighbors[i,]
    mu[i] <- mean(abs(object$oob_error[ind]))
  }

  if(is.null(samples)){
    # Generate predictions and confidence intervals
    moe <- matrix(NA, nrow=length(pred), ncol=length(conf))
    colnames(moe) <- gsub("\\.", "_", paste(conf))
    for(i in seq_along(conf)){
      moe[,i] <- quantile(object$alpha, conf[i])*(mu + object$beta)
    }
    out <- list(pred=pred, moe=moe)
    return(out)
  }else{
    # Generate predictive samples
    preds <- matrix(NA, nrow=samples, ncol=n)
    for(i in 1:n){
      alpha_delta <- sample(object$alpha, samples, replace=TRUE)
      preds[,i] <- pred[i] + alpha_delta*(mu[i] + object$beta)
    }
  }
  return(preds)
}
