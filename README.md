
<!-- README.md is generated from README.Rmd. Please edit that file -->

# quack - R tools for UQ and calibration

<div class="figure">

<img src="inst/logos/QUACK.png" alt="This logo was designed by Imagine AI Art Studio" width="40%" />
<p class="caption">
This logo was designed by Imagine AI Art Studio
</p>

</div>

### Description

The `quack` package contains several tools for UQ and Bayesian model
calibration, with an emphasis on inference for physical parameters in
complex physical systems. Many of the tools found in this package are
implementations of the ideas found
[here](https://digitalrepository.unm.edu/cgi/viewcontent.cgi?article=1165&=&context=math_etds&=&sei-redir=1&referer=https%253A%252F%252Fscholar.google.com%252Fscholar%253Fhl%253Den%2526as_sdt%253D0%25252C32%2526q%253DKellin%252BRumsey%252Bphysical%252Bparameters%2526btnG%253D#search=%22Kellin%20Rumsey%20physical%20parameters%22)

### Note (7/17/23)

The R package `MHadaptive` has been removed from CRAN. As of 7/17/2023,
in order to install the `quack` package, users will first need to
install `MHadaptive` from the Archive as

``` r
#install.packages("devtools")
link <- "https://cran.r-project.org/src/contrib/Archive/MHadaptive/MHadaptive_1.1-8.tar.gz"
devtools::install_url(link)
```

## Installation

To install the `quack` package, type

``` r
# install.packages("devtools")
devtools::install_github("knrumsey/quack")
```

## Tools

See
[manual](https://github.com/knrumsey/quack/blob/master/quack_0.0.0.9000.pdf)
for details.

- [**Moment
  penalization**](https://epubs.siam.org/doi/pdf/10.1137/19M1283707?casa_token=javZzkQnG4oAAAAA:NqF-i_Wpuz5I8j0IHYK-j-q4QzoJr04ohxO1PBMHLwKE640bTAD1MsalHKtFBu1-VTsuu4sR):
  functions for the moment penalization (MP) prior, which are useful for
  improving identifiability and reducing bias when a physical system
  contains a set of nuisance parameters which can be viewed as iid
  samples from a specified distribution. Also contains methods for
  computing the probability of prior coherency, a useful diagnostic
  tool.
- **Emulating the conditional posterior**: A fast approach for
  approximating the [cut
  distribution](https://idp.springer.com/authorize/casa?redirect_uri=https://link.springer.com/article/10.1007/s11222-014-9503-z&casa_token=K78vRWO3qrQAAAAA:12SChC4sxODP5fp0NYdGSQ-1YybCxPndiT6JZb_NXdKtVlf3-eRa87JvfULTvZs1D30yMOJkWZfOpdQ)
  in low to moderate dimensions.
- [**Fast matrix algebra for
  BMC**](https://www.tandfonline.com/doi/pdf/10.1080/00949655.2020.1850729?casa_token=aa1qKL4wpeoAAAAA:R5KvrccY4aWIDYLi5OqQKYHCfV5JMK-nEKytBwjmAqwSgW8miqZMlfqs9yStkQDXqtW3PHg74j4):
  Near-quadratic time approximations to the inverse of the covariance
  matrix given a particular structure. Useful for model calibration when
  a large number of sequential inversions are needed.
- **Sequential local approximate GPs**: Implementations of `slapGP` and
  `leapGP` which are two “global-model” extensions of the
  [laGP](https://www.jstatsoft.org/index.php/jss/article/view/v072i01/1030)
  framework.
- **Accelerated bootstrap**: An implementation of the accelerated
  bootstrap.
- [**Joint credible
  regions**](https://stats.stackexchange.com/questions/361350/joint-credible-regions-from-mcmc-draws):
  Function to find elliptical joint credible regions given a sample of
  points $x \in \mathbb R^p$,

## References

Rumsey, Kellin N. Methods of uncertainty quantification for physical
parameters. Diss. The University of New Mexico, 2020.

Rumsey, Kellin, et al. “Dealing with measurement uncertainties as
nuisance parameters in Bayesian model calibration.” SIAM/ASA Journal on
Uncertainty Quantification 8.4 (2020): 1287-1309.

Plummer, Martyn. “Cuts in Bayesian graphical models.” Statistics and
Computing 25.1 (2015): 37-43.

Rumsey, Kellin N., and Gabriel Huerta. “Fast matrix algebra for Bayesian
model calibration.” Journal of Statistical Computation and Simulation
91.7 (2021): 1331-1341.

Gramacy, Robert B. “laGP: large-scale spatial modeling via local
approximate Gaussian processes in R.” Journal of Statistical Software 72
(2016): 1-46
