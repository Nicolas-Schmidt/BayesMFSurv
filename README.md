
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

| Function                      | Description                                                                                                                                                      |
| ----------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `mfsurv`                      | fits a parametric Bayesian MF model via Markov Chain Monte Carlo (MCMC) to estimate the misclassification in the first stage and the hazard in the second stage. |
| `mfsurv.stats`                | A function to calculate the deviance information criterion (DIC) for fitted model objects of class mfsurv                                                        |
| `mfsurv.summary`              | Returns a summary of a mfsurv object via `coda::summary.mcmc`                                                                                                    |
| `mcmcsurv`                    | estimates a Bayesian Exponential or Weibull model via Markov Chain Monte Carlo (MCMC)                                                                            |
| `betas.post2`                 | log-posterior distribution of betas with pth element fixed as betas.p                                                                                            |
| `betas.slice.sampling2`       | slice sampling for betas                                                                                                                                         |
| `univ.betas.slice.sampling2`  | univariate slice sampling for betas.p                                                                                                                            |
| `gammas.post2`                | log-posterior distribution of gammas with pth element fixed as gammas.p                                                                                          |
| `gammas.slice.sampling2`      | slice sampling for gammas                                                                                                                                        |
| `univ.gammas.slice.sampling2` | univariate slice sampling for gammas.p                                                                                                                           |
| `lambda.post2`                | log-posterior distribution of lambda                                                                                                                             |
| `lambda.slice.sampling2`      | univariate slice sampling for lambda                                                                                                                             |

### Example

The data used to estimate the examples that follow comes from Reenock,
Bernhard and Sobek (2007) -DOI: 10.111/j.1468-2478.2007.00469.x-. The
RBS (2007) dataset uses continuous-time event history techniques to code
episodes of democratic breakdown in all democracies from 1961 to 1995.
In addition, it provides data on a number of economic and political
variables.

| Variable     | Description                                |
| ------------ | ------------------------------------------ |
| **calinv**   | inverse of caloric intake                  |
| **lnlevel**  | gross domestic product per capita (logged) |
| **calileve** | interaction calinv\*lnlevel                |
| **necon**    | economic growth                            |
| **presi**    | presidential regime                        |
| **tag**      | effective number of parties                |
| **rel**      | religious fractionalization                |
| **ethn**     | ethnic fractionalization                   |
| **prevdem**  | \# of previous democratic episodes         |
| **openc**    | trade openness                             |

``` r
library(BayesMFSurv)

# Baseline Bayesian misclassified failure (MF) model. 
# Misclassification stage only includes the intercept while the survival stage 
# includes all covariates described above.  

library(BayesMFSurv)

set.seed(95)
rbs <- na.omit(rbs)
Y   <- rbs$Y
X   <- as.matrix(cbind(1, rbs[,1:10]))
C   <- rbs$C
Z1  <- cbind(rep(1,nrow(rbs)))
Y0  <- rbs$Y0
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

#### Citation

To cite package`BayesMFSurv` in publications, please use:

``` r
citation(package = 'BayesMFSurv')
```
