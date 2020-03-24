
<!-- README.md is generated from README.Rmd. Please edit that file -->

## The `BayesMFSurv` package

*Minnie M. Joo, Nicolas Schmidt, Sergio Bejar, Bumba Mukherjee, Vineeta
Yadav*

<!-- badges: start -->

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/BayesMFSurv)](https://cran.r-project.org/package=BayesMFSurv)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN RStudio
mirrordownloads](https://cranlogs.r-pkg.org/badges/grand-total/BayesMFSurv?color=blue)](https://www.r-pkg.org/pkg/BayesMFSurv)
[![CRAN RStudio
mirrordownloads](https://cranlogs.r-pkg.org/badges/BayesMFSurv?color=blue)](https://www.r-pkg.org/pkg/BayesMFSurv)
<!-- badges: end -->

### Description

Contains a split population survival estimator that models the
misclassification probability of failure versus right-censored events.
The split population survival estimator is described in Bagozzi et
al. (2019) <doi:10.1017/pan.2019.6>.

### Installation

You can install the released version (`0.1.0`) of `BayesMFSurv` from
[CRAN](https://cran.r-project.org/) with:

``` r
install.packages("BayesMFSurv")
```

And the development version (`0.2.0`) from GitHub with:

``` r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("Nicolas-Schmidt/BayesMFSurv")
```

## Functions

| Function   | Description                                                                                                                                                                                                                                                                                 |
| ---------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `mfsurv`   | fits a parametric Bayesian MF model via Markov Chain Monte Carlo (MCMC) to estimate the probability of misclassification in the first stage and the hazard in the second stage. Slice sampling is employed to draw the posterior sample of the model’s split and survival stage parameters. |
| `mcmcsurv` | estimates a Bayesian Exponential or Weibull survival model via Markov Chain Monte Carlo (MCMC). Slice samplig is employed to draw the posterior sample of the model’s survival stage parameters.                                                                                            |
| `stats`    | a function to calculate the deviance information criterion (DIC) and the log-likelihood for fitted model objects of class mfsurv or mcmcsurv.                                                                                                                                               |
| `summary`  | returns a summary of a mfsurv or mcmcsurv object via `coda::summary.mcmc`.                                                                                                                                                                                                                  |

### Example

The data used to estimate the following examples come from Reenock,
Bernhard and Sobek (2007) -DOI: 10.111/j.1468-2478.2007.00469.x-. The
RBS (2007) dataset uses continuous-time event history techniques to code
episodes of democratic breakdown in all democracies from 1961 to 1995.
In addition, it provides data on a number of economic and political
variables.

| Variable     | Description                                         |
| ------------ | --------------------------------------------------- |
| **calinv**   | inverse of per capita daily caloric supply          |
| **lnlevel**  | natural log of economic development                 |
| **calileve** | inverse of per capita daily caloric supply\*lnlevel |
| **necon**    | economic performance                                |
| **presi**    | presidential regime                                 |
| **tag**      | effective number of parties                         |
| **rel**      | religious fractionalization                         |
| **ethn**     | ethnic fractionalization                            |
| **prevdem**  | numbers of previous democratic episodes             |
| **openc**    | trade openness                                      |
| **Y**        | years in current democratic episode                 |
| **Y0**       | years in current democratic episode (lagged)        |
| **C**        | breakdown of democratic episode                     |

#### Misclassified-Failure

`mfsurv` estimated the probability of misclassification failure in the
first (split) stage and hazard in the second (survival) stage.

`mfsurv` should be used when user suspects that data of survival cases
could be right-sensored (i.e. when there is a probability that failure
events are misclassified).

Example with N = 100000 is
[here](https://github.com/Nicolas-Schmidt/BayesMFSurv/tree/master/data-raw).

``` r

# Baseline Bayesian misclassified failure (MF) model. 
# Misclassification stage only includes the intercept while the survival stage 
# includes all covariates described above.  

library(BayesMFSurv)

set.seed(95)
data(RBS)
RBS <- na.omit(RBS)
Y   <- RBS$Y
X   <- as.matrix(cbind(1, RBS[,1:10]))
C   <- RBS$C
Z1  <- cbind(rep(1,nrow(RBS)))
Y0  <- RBS$Y0
model1 <- mfsurv(Y ~ X | C ~ Z1, Y0 = Y0,
                 N = 100,
                 burn = 50,
                 thin = 5,
                 w = c(0.5, 0.5, 0.5),
                 m = 20,
                 form = 'Weibull')


stats(model1)
#> $DIC
#> [1] 256.8576
#> 
#> $Loglik
#> [1] 478.8454

summary(model1, parameter = c("betas"))
#> 
#> Iterations = 1:10
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 10 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean    SD Naive SE Time-series SE
#> X.intercept  -5.0477 3.384   1.0700         1.0700
#> X1          -10.6678 3.797   1.2006         1.2006
#> Xcalinv      12.8464 7.331   2.3182         6.9599
#> Xlnlevel     -3.2283 3.595   1.1367         1.1367
#> Xcalileve     0.9831 2.667   0.8435         0.8369
#> Xnecon       -8.9798 4.207   1.3302         2.9719
#> Xpresi       -2.4573 3.032   0.9588         0.9588
#> Xtag         -3.7535 2.600   0.8222         0.8222
#> Xrel         -8.1554 3.174   1.0037         2.0935
#> Xethn         2.2446 3.124   0.9878         1.9027
#> Xprevdem     13.8614 2.904   0.9182         0.9182
#> Xopenc       -2.4537 1.367   0.4324         0.4324
#> 
#> 2. Quantiles for each variable:
#> 
#>                2.5%      25%    50%     75%   97.5%
#> X.intercept  -9.040  -7.8035 -5.141 -3.3208  0.6360
#> X1          -17.594 -12.4190 -9.431 -8.0545 -6.4100
#> Xcalinv       2.862   5.6273 15.445 19.0850 20.6433
#> Xlnlevel     -8.261  -5.2289 -3.487 -1.6446  3.1120
#> Xcalileve    -2.559  -1.2788  1.068  2.5555  5.1182
#> Xnecon      -15.669 -10.7354 -9.695 -5.7098 -3.5264
#> Xpresi       -8.147  -3.5865 -2.260 -0.5025  1.2235
#> Xtag         -8.032  -4.6135 -3.771 -1.3903 -0.7852
#> Xrel        -12.638 -10.7618 -7.601 -6.1422 -3.8367
#> Xethn        -1.735  -0.6924  2.470  4.3050  6.8935
#> Xprevdem      9.677  11.7561 13.544 16.2759 18.0732
#> Xopenc       -4.466  -3.6585 -1.993 -1.5582 -0.9043
```

#### Non Misclassified-Failure

`mcmcsurv` estimates a Bayesian equivalent of standard survival models
(i.e. Exponential or Weibull).

Example with N = 15000 can be found
[here](https://github.com/Nicolas-Schmidt/BayesMFSurv/tree/master/data-raw).

``` r
set.seed(95)
model2 <- mcmcsurv(Y = Y, Y0 = Y0, C =  C,  X = X, 
                   N = 100, 
                   burn = 50, 
                   thin = 5, 
                   w = c(0.5, 0.5, 0.5),
                   m = 20, 
                   form = 'Weibull')


stats(model2)
#> $DIC
#> [1] -976.9334
#> 
#> $Loglik
#> [1] 444.0505

summary(model2, parameter = c("betas"))
#> 
#> Iterations = 1:10
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 10 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>               Mean       SD Naive SE Time-series SE
#> X1        -0.49109 0.746752 0.236144       0.416615
#> Xcalinv    0.08710 2.488545 0.786947       0.786947
#> Xlnlevel  -0.60347 0.090054 0.028477       0.034769
#> Xcalileve -1.17459 1.770004 0.559724       0.426229
#> Xnecon    -2.65799 2.004900 0.634005       0.634005
#> Xpresi     0.45379 0.336383 0.106374       0.106374
#> Xtag       0.05672 0.085752 0.027117       0.027117
#> Xrel       1.02890 0.462122 0.146136       0.146136
#> Xethn      1.93289 0.569031 0.179943       0.432835
#> Xprevdem   0.08207 0.262127 0.082892       0.082892
#> Xopenc    -0.01040 0.005484 0.001734       0.001872
#> 
#> 2. Quantiles for each variable:
#> 
#>               2.5%      25%      50%       75%     97.5%
#> X1        -1.95461 -0.53040 -0.27203 -0.048772  0.269946
#> Xcalinv   -2.88002 -1.91528 -0.09226  0.964488  4.146867
#> Xlnlevel  -0.73460 -0.65810 -0.61260 -0.556346 -0.462481
#> Xcalileve -4.29853 -1.86115 -1.24021 -0.137511  1.121879
#> Xnecon    -6.10899 -4.17259 -1.91485 -1.281245 -0.347827
#> Xpresi     0.19266  0.23482  0.30154  0.629261  1.082617
#> Xtag      -0.07137  0.01431  0.05023  0.089369  0.204853
#> Xrel       0.48865  0.71068  0.96574  1.158598  1.863872
#> Xethn      1.10966  1.42295  2.13381  2.190843  2.733231
#> Xprevdem  -0.30593 -0.10253  0.08783  0.327342  0.387947
#> Xopenc    -0.01962 -0.01124 -0.01010 -0.007111 -0.003705
```

#### Citation

To cite package`BayesMFSurv` in publications, please use:

``` r
citation(package = 'BayesMFSurv')
```
