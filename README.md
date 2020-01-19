
<!-- README.md is generated from README.Rmd. Please edit that file -->

## The `BayesMFSurv` package

*Minnie M. Joo, Sergio Bejar, Nicolas Schmidt, Bumba
Mukherjee*

<!-- badges: start -->

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/BayesMFSurv)](https://cran.r-project.org/package=BayesMFSurv)
[![CRAN RStudio
mirrordownloads](https://cranlogs.r-pkg.org/badges/BayesMFSurv?color=brightgreen)](https://www.r-pkg.org/pkg/BayesMFSurv)
<!-- badges: end -->

### Description

Contains a split population survival estimator that models the
misclassification probability of failure versus right-censored events.
The split population survival estimator is described in Bagozzi et
al. (2019) <doi:10.1017/pan.2019.6>.

### Installation

``` r
# Install BayesMFSurv from CRAN
install.packages("BayesMFSurv")

# The development version from GitHub:
if (!require("remotes")) install.packages("remotes")
remotes::install_github("Nicolas-Schmidt/BayesMFSurv")
```

## Functions

| Function                      | Description                                                                                                                                                      |
| ----------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `mfsurv`                      | fits a parametric Bayesian MF model via Markov Chain Monte Carlo (MCMC) to estimate the misclassification in the first stage and the hazard in the second stage. |
| `mfsurv.stats`                | A function to calculate the deviance information criterion (DIC) for fitted model objects of class mfsurv                                                        |
| `mfsurv.summary`              | Returns a summary of a mfsurv object via `coda::summary.mcmc`                                                                                                    |
| `betas.post2`                 | log-posterior distribution of betas with pth element fixed as betas.p                                                                                            |
| `betas.slice.sampling2`       | slice sampling for betas                                                                                                                                         |
| `univ.betas.slice.sampling2`  | univariate slice sampling for betas.p                                                                                                                            |
| `gammas.post2`                | log-posterior distribution of gammas with pth element fixed as gammas.p                                                                                          |
| `gammas.slice.sampling2`      | slice sampling for gammas                                                                                                                                        |
| `univ.gammas.slice.sampling2` | univariate slice sampling for gammas.p                                                                                                                           |
| `lambda.post2`                | log-posterior distribution of lambda                                                                                                                             |
| `lambda.slice.sampling2`      | univariate slice sampling for lambda                                                                                                                             |

### Example

``` r
library(BayesMFSurv)

set.seed(95)
bgl <- Buhaugetal_2009_JCR
bgl <- subset(bgl, coupx == 0)
bgl <- na.omit(bgl)
Y   <- bgl$Y
X   <- as.matrix(cbind(1, bgl[,1:7]))
C   <- bgl$C
Z1  <- matrix(1, nrow = nrow(bgl))
Y0  <- bgl$Y0
model1 <- mfsurv(Y ~ X | C ~ Z1, Y0 = Y0,
                 N = 150,
                 burn = 50,
                 thin = 30,
                 w = c(0.1, .1, .1),
                 m = 10,
                 form = "Weibull",
                 na.action = 'na.omit')


mfsurv.summary(model1, parameter = c("betas"))
#> 
#> Iterations = 1:3
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 3 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                  Mean     SD Naive SE Time-series SE
#> X.intercept    1.3526 0.5162   0.2980         0.2980
#> X1             4.5132 1.9285   1.1134         1.1134
#> Xlndistx      -1.4516 0.7722   0.4458         0.4458
#> Xconfbord    -10.0149 2.2335   1.2895         1.2895
#> Xborddist     -1.1794 3.4491   1.9913         1.9913
#> Xfigcapdum     1.7497 1.0601   0.6121         0.6121
#> Xlgdp_onset   -5.8896 0.3726   0.2151         0.2151
#> Xsip2l_onset  -4.4748 0.8591   0.4960         0.4960
#> Xpcw           0.8724 1.4270   0.8239         0.8239
#> 
#> 2. Quantiles for each variable:
#> 
#>                  2.5%      25%     50%     75%   97.5%
#> X.intercept    0.9706   1.0602   1.160  1.5485  1.8985
#> X1             3.2518   3.4041   3.573  5.1524  6.5736
#> Xlndistx      -2.2769  -1.6799  -1.017 -1.0058 -0.9961
#> Xconfbord    -11.3148 -11.3044 -11.293 -9.3643 -7.6287
#> Xborddist     -4.4765  -2.8915  -1.130  0.5573  2.0762
#> Xfigcapdum     1.0414   1.1410   1.252  2.1094  2.8814
#> Xlgdp_onset   -6.2076  -6.0913  -5.962 -5.7241 -5.5098
#> Xsip2l_onset  -5.0866  -4.9646  -4.829 -4.1622 -3.5619
#> Xpcw          -0.5952   0.2549   1.199  1.6535  2.0621
```

#### Citation

To cite package`BayesMFSurv` in publications, please use:

``` r
citation(package = 'BayesMFSurv')
```
