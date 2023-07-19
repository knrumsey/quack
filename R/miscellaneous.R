#' Customized "smart" (robust) LHS function
#'
#' A default lhs function to use (to avoid accidental crashing of maximinLHS)
#'
#' @param n the number of rows or samples
#' @param k the number of columns or parameters/variables
#' @param maximin_n the number of rows or samples in the original maximin design
#' @param ... additional arguments passed to randomLHS (n > 50k) or maximinLHS (n <= 50k)
#' @return a latin hypercube design
#' @details Tries to use a maxmin LHS design when possible (n <= maximin_n = 1000).
#' If n is greater than 100,000, we simply return a random LHS.
#' In between these cases, we build an initial maximin design of size maximin_n and then augment this design with augmentLHS.
#' @examples
#' x <- myLHS(1000, 8)
#' y <- apply(x, 1, quack::ff_borehole)
#' @export
smartLHS <- function(n, k, maximin_n=1e3, ...){
  if(n > 1e5){
    return(randomLHS(n, k, ...))
  }
  if(n <= maximin_n){
    return(maximinLHS(n, k, ...))
  }
  n2 <- n - maximin_n
  x <- maximinLHS(maximin_n, k, ...)
  return(augmentLHS(x, n2))
}

#' Kernel Density Estimate with Beta Kernels
#'
#' A kernel density estimate which will always integrate to 1 on the interval (0, 1)
#'
#' @param x the data from which the estimate is to be computed.
#' @param v a smoothness parameter (smaller is smoother)
#' @param n the number of points in (0, 1) on which the density is evaluated
#' @return a list with components x and y
#' @details KDE with beta kernels centered at each point in x. v is the "sample size" parameter of the beta.
#' @examples
#' x = rbeta(150, 3, 2)
#' hist(x, freq=F, breaks=30, xlim=c(0,1))
#' lines(density01(x, v=50), lwd=2, col='firebrick')
#' lines(density01(x, v=100), lwd=2, col='orange')
#' lines(density01(x, v=1000), lwd=2, col='dodgerblue')
#' @export
density01 <- function(x, v=100, n=512){
  nx <- length(x)
  xx <- seq(0, 1, length.out=n)
  dens <- rep(0, n)
  for(i in 1:nx){
    v0 <- v
    a <- b <- 0
    while(min(a, b) < 1){
      a = x[i]*v0
      b = (1-x[i])*v0
      v0 <- v0*1.05
    }
    dens <- dens + dbeta(xx, a, b)/nx
  }
  out <- list(x=xx, y=dens)
  class(out) <- "density"
  return(out)
}


#' A summary function to make things readily available for ggplot
#'
#' Requires plyr pacakge
#'
#' @param data the data frame
#' @param measurevar the response variable (in quotes)
#' @param groupvars a vector of grouping variables (in quotes)
#' @param na.rm should na values be removed
#' @param conf.interval coverage probability desired for t-intervals
#' @return a data.frame of mean,sd,ci,n grouped by groupvars
#' @details A function for easy handling of data frame to be used for ggplot2
#' @examples
#' set.seed(111)
#' df <- data.frame(foo = sample(c("a", "b", "c"), 100, replace=TRUE),
#'                  fie = rpois(100, 0.5))
#' df$y <- rnorm(100, df$fie + 1, 0.1)*(df$foo == "a") +
#'         rnorm(100, df$fie + 2, 0.1)*(df$foo == "b") +
#'         rnorm(100, df$fie + 3, 0.1)*(df$foo == "c")
#'
#' summ <- ggsummary(df, measurevar = "y", groupvars = c("foo", "fie"))
#' #ggplot(summ, aes(x=fie, y=mean, fill=foo)) +
#' #     geom_bar(position=position_dodge2(padding=0.1), stat="identity") +
#' #     geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci),
#' #   width=.05,
#' #     position=position_dodge(.9))
#'
#' @export
ggsummary <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {


  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )

  # Rename the "mean" column
  #datac <- plyr::rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}
