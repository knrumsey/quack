#' @name test_functions
#' @rdname test_functions
#'
#' @title Some general testing functions for UQ
#'
#' @param x a vector in inputs. type get_input_length("function_name") to get the minimum requirement for length(x)
#' @param scale01 default TRUE. Expects inputs on [0,1] space. When FALSE, inputs are expected in their natural ranges
#' @details List of functions and details
#' \enumerate{
#'    \item borehole
#'    \item [lim](https://www.sfu.ca/~ssurjano/limetal02pol.html)
#'    \item [gramacy_lee](https://www.sfu.ca/~ssurjano/grlee08.html)
#'    \item synth_3level
#'    }
#'
NULL

#' @rdname test_functions
#' @examples
#' n <- 10
#' x <- matrix(runif(3*n), nrow=n)
#' y <- apply(x, 1, f_lim)
#' @export
ff_borehole <- function(x, scale01=TRUE){
  if(scale01){
    lb <- c(0.05, 100, 63070, 990, 63.1, 700, 1220, 9855)
    ub <- c(0.15, 40000, 115600, 1110, 116, 820, 1680, 12045)
    x[1:8] <- x[1:8]*(ub-lb) + lb
  }
  xx <- x
  rw <- xx[1]
  r  <- xx[2]
  Tu <- xx[3]
  Hu <- xx[4]
  Tl <- xx[5]
  Hl <- xx[6]
  L  <- xx[7]
  Kw <- xx[8]

  frac1 <- 2 * pi * Tu * (Hu-Hl)

  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)

  y <- frac1 / frac2
  return(y)
}

#' @rdname test_functions
#' @export
ff_lim <- function(x){
  xx <- x
  (18 +
    5*xx[1] -
    35*xx[2] +
    5*xx[1]*xx[2] +
    38*xx[2]^2 -
    15*xx[1]^3 -
    5*xx[1]*xx[2]^2 -
    11*xx[2]^4 +
    xx[1]^3*xx[2]^2)/2
}

#' @rdname test_functions
#' @export
ff_gl2008 <- function(x, scale01=TRUE){
  if(scale01){
    xx <- x*8 -2
  }else{
    xx <- x
  }
  x[1]*exp(-x[1]^2 - x[2]^2)
}

#' @rdname test_functions
#' @export
ff_ripples <- function(x){
  cc <- c(-0.1379894,
          -0.3483005,
          1.6280129,
          -1.4515784,
          0.9552673,
          1.9053080,
          1.5198398,
          -1.5328629,
          0.1851551,
          -1.4395361)
  res <- 0
  for(i in 1:5){
    res <- res + sin((2*i+1)*pi*(cc[2*i-1]*x[1]+cc[2*i]*x[2]))
  }
  return(res)
}


#' @rdname test_functions
#' @export
ff_synth_3level <- function(x){
  0.55*ff_lim(x) + 5*ff_gl2008(x) + 0.05*ff_ripples(x)
}

#' @rdname test_functions
#' @export
get_input_length <- function(function_name){
  if(function_name=="ff_borehole") return(8)
  if(function_name=="ff_ripples") return(2)
  if(function_name=="ff_gl2008") return(2)
  if(function_name=="ff_lim") return(2)
  if(function_name=="ff_synth_3level") return(2)
  stop("function not recognized")
}






