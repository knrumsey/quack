#' G Test
#'
#' Performs G-test for an arbitrary D-dimensional array
#'
#' @param arr the true scalar value
#' @param df If not specified, test uses the standard mutual independence assumption: prod(dims) - sum(dims-1) - 1.
#' @references Fill this in.
#' @examples
#' data_2d <- cbind(rpois(3, lambda=10), rpois(3, lambda=10), rpois(3, lambda=10))
#' fisher.test(data_2d, alternative="greater")
#' chisq.test(data_2d)
#' Gtest(data_2d)
#'
#' data_3d <- array(rpois(2*3*4, lambda=10), dim=c(2,3,4))
#' Gtest(data)

#' @export
Gtest <- function(arr, df=NA){
  dims <- dim(arr)
  D <- length(dims)
  N <- sum(arr)
  if(is.na(df)){
    df <- prod(dims) - sum(dims-1) - 1
    if(df <= 0){
      stop("degrees of freedom must be positive")
    }
    if(N <= df){
      warning("N is less than df")
    }
  }

  # Get marginal totals
  marginals <- lapply(1:D, function(d) apply(arr, d, sum))
  total_cells <- prod(dims)
  expected <- array(NA, dim=dims)

  #Loop over cells
  for(i in 1:total_cells){
    indices <- arrayInd(i, .dim=dims)

    #Calculate expected counts
    numerator <- 1
    for(d in 1:D){
     numerator <- numerator * marginals[[d]][indices[d]]
    }
    E_ind <-  numerator / N^(D-1)
    expected <- do.call("[<-", c(list(expected), as.list(indices), list(value = E_ind))) # Multi-indexing
  }

  if(min(expected) <= 0){
    warning("cannot have expected counts be <= 0. Check for level with no observations.")
  }

  # Compute G-statistico
  observed <- arr
  mask <- observed > 0
  Gstatistic <- 2 * sum(observed[mask] * log(observed[mask] / expected[mask]))

  # Step 4: Compute p-value
  pvalue <- pchisq(Gstatistic, df = df, lower.tail = FALSE)

  # Return object
  out <- list(G=Gstatistic, pval=pvalue,
              observed=observed,
              expected=expected)
  class(out) <- "Gtest"
  return(out)
}

#' @export
print.Gtest <- function(x, ...){
  print(paste0("G = ", x$Gstatistic))
  print(paste0("p = ", x$pval))
}
