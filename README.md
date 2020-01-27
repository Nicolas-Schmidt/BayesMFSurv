
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

| Function                      | Description                                                                                                                                                                     |
| ----------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `mfsurv`                      | fits a parametric Bayesian MF model via Markov Chain Monte Carlo (MCMC) to estimate the probability of misclassification in the first stage and the hazard in the second stage. |
| `mfsurv.stats`                | a function to calculate the deviance information criterion (DIC) and the log-likelihood for fitted model objects of class mfsurv                                                |
| `mfsurv.summary`              | returns a summary of a mfsurv object via `coda::summary.mcmc`                                                                                                                   |
| `betas.post2`                 | log-posterior distribution of betas with pth element fixed as betas.p                                                                                                           |
| `betas.slice.sampling2`       | slice sampling for betas                                                                                                                                                        |
| `univ.betas.slice.sampling2`  | univariate slice sampling for betas.p                                                                                                                                           |
| `gammas.post2`                | log-posterior distribution of gammas with pth element fixed as gammas.p                                                                                                         |
| `gammas.slice.sampling2`      | slice sampling for gammas                                                                                                                                                       |
| `univ.gammas.slice.sampling2` | univariate slice sampling for gammas.p                                                                                                                                          |
| `lambda.post2`                | log-posterior distribution of lambda                                                                                                                                            |
| `lambda.slice.sampling2`      | univariate slice sampling for lambda                                                                                                                                            |
| `mcmcsurv`                    | estimates a Bayesian Exponential or Weibull survival model via Markov Chain Monte Carlo (MCMC)                                                                                  |

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

``` r
library(BayesMFSurv)

# Baseline Bayesian misclassified failure (MF) model. 
# Misclassification stage only includes the intercept while the survival stage 
# includes all covariates described above.  

library(BayesMFSurv)

set.seed(95)
RBS <- na.omit(RBS)
Y   <- RBS$Y
X   <- as.matrix(cbind(1, RBS[,1:10]))
C   <- RBS$C
Z1  <- cbind(rep(1,nrow(RBS)))
Y0  <- RBS$Y0
model1 <- mfsurv(Y ~ X | C ~ Z1, Y0 = Y0,
                 N = 100000,
                 burn = 10000,
                 thin = 100,
                 w = c(0.5, .5, .5),
                 m = 20,
                 form = "Weibull",
                 na.action = 'na.omit')
```

``` r
mfsurv.stats(model1)
#> $DIC
#> [1] -91438.64
#> 
#> $Loglik
#> [1] 1243.269

mfsurv.summary(model1, parameter = c("betas"))
#> 
#> Iterations = 1:900
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 900 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                 Mean     SD Naive SE Time-series SE
#> X.intercept  1.59992 2.5488  0.08496        0.11147
#> X1           1.56001 1.6683  0.05561        0.08020
#> Xcalinv      0.16505 2.5411  0.08470        0.19144
#> Xlnlevel    -1.12098 1.4762  0.04921        0.09873
#> Xcalileve   -0.03761 2.4731  0.08244        0.15828
#> Xnecon      -2.13090 2.4107  0.08036        0.07463
#> Xpresi       0.23041 1.9100  0.06367        0.05790
#> Xtag        -0.01230 0.7878  0.02626        0.03822
#> Xrel         1.45198 1.7404  0.05801        0.04774
#> Xethn        1.01890 1.3000  0.04333        0.03564
#> Xprevdem     0.13070 1.2525  0.04175        0.11232
#> Xopenc      -0.19741 1.7926  0.05975        0.18939
#> 
#> 2. Quantiles for each variable:
#> 
#>                 2.5%      25%      50%     75%     97.5%
#> X.intercept -1.13533  0.69903  1.52827  2.5172  4.865523
#> X1          -1.40643  0.62558  1.55816  2.4422  4.804547
#> Xcalinv     -3.85120 -1.13823  0.00846  1.2496  4.676094
#> Xlnlevel    -1.42372 -1.12625 -1.01540 -0.8972 -0.676814
#> Xcalileve   -3.85666 -1.07941  0.04700  1.2389  4.015589
#> Xnecon      -6.11097 -3.09541 -2.00216 -0.9263  1.233515
#> Xpresi      -0.38909  0.05153  0.26004  0.4705  0.945025
#> Xtag        -0.25624 -0.06322  0.03359  0.1123  0.255185
#> Xrel         0.07307  0.93358  1.38643  1.8464  2.834912
#> Xethn       -0.32831  0.54553  0.99991  1.4344  2.498845
#> Xprevdem    -0.66111 -0.15803  0.05060  0.2317  0.689183
#> Xopenc      -0.03202 -0.01575 -0.01156 -0.0081 -0.001084
```

#### Non Misclassified-Failure

``` r
model2 <- mcmcSurv(Y = Y, Y0 = Y0, C =  C,  X = X, 
                   N = 5000, 
                   burn = 500, 
                   thin = 5, 
                   w = c(0.5, 0.5, 0.5),
                   m = 10, 
                   form = "Weibull")
```

``` r
str(model2)
#> List of 2
#>  $ betas: num [1:900, 1:11] 3.22 3.04 2.85 3.07 3.12 ...
#>  $ rho  : num [1:900] 1.032 0.948 1.384 1.322 1.467 ...
mfsurv.summary(model2, parameter = c("betas"))
#> 
#> Iterations = 1:900
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 900 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>           Mean       SD  Naive SE Time-series SE
#>  [1,]  2.72471 1.324951 0.0441650      0.2703780
#>  [2,] -0.05396 2.124243 0.0708081      0.1220671
#>  [3,] -0.96244 0.167482 0.0055827      0.0400096
#>  [4,]  0.16176 2.254210 0.0751403      0.1362043
#>  [5,] -2.19964 1.753517 0.0584506      0.0843977
#>  [6,]  0.27741 0.293103 0.0097701      0.0109187
#>  [7,]  0.02540 0.116009 0.0038670      0.0073200
#>  [8,]  1.39552 0.632141 0.0210714      0.0361141
#>  [9,]  1.13872 0.604012 0.0201337      0.0450722
#> [10,]  0.04308 0.300894 0.0100298      0.0123460
#> [11,] -0.01116 0.005696 0.0001899      0.0002637
#> 
#> 2. Quantiles for each variable:
#> 
#>            2.5%      25%      50%      75%     97.5%
#> var1   0.318390  1.82712  2.64239  3.40875  5.866844
#> var2  -4.436232 -1.28413 -0.02135  1.13929  4.312856
#> var3  -1.357789 -1.05501 -0.94940 -0.85087 -0.670396
#> var4  -4.281074 -1.17035  0.13294  1.37366  5.164420
#> var5  -5.965880 -3.33416 -1.97787 -0.96312  0.721732
#> var6  -0.293328  0.09357  0.27259  0.47554  0.867117
#> var7  -0.218787 -0.04463  0.03039  0.10790  0.227401
#> var8   0.284365  0.97018  1.38005  1.79257  2.701030
#> var9   0.009273  0.69678  1.16980  1.53924  2.315710
#> var10 -0.600713 -0.14257  0.05590  0.24628  0.634237
#> var11 -0.023211 -0.01473 -0.01076 -0.00732 -0.001108
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-7-2.png" width="100%" />

#### Citation

To cite package`BayesMFSurv` in publications, please use:

``` r
citation(package = 'BayesMFSurv')
```
