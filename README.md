
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

``` r

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
stats(model1)
#> $DIC
#> [1] 466.1642
#> 
#> $Loglik
#> [1] 1243.269

summary(model1, parameter = c("betas"))
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
set.seed(95)
model2 <- mcmcsurv(Y = Y, Y0 = Y0, C =  C,  X = X, 
                   N = 15000, 
                   burn = 500, 
                   thin = 10, 
                   w = c(0.5, 0.5, 0.5),
                   m = 20, 
                   form = "Weibull")
```

``` r
stats(model2)
#> $DIC
#> [1] -937.9898
#> 
#> $Loglik
#> [1] 406.0937

summary(model2, parameter = c("betas"))
#> 
#> Iterations = 1:1450
#> Thinning interval = 1 
#> Number of chains = 1 
#> Sample size per chain = 1450 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>               Mean       SD  Naive SE Time-series SE
#> X1         2.51481 1.709420 0.0448916      0.3005889
#> Xcalinv   -0.06664 2.225281 0.0584387      0.0862545
#> Xlnlevel  -0.94028 0.199542 0.0052402      0.0335193
#> Xcalileve  0.02971 2.162636 0.0567936      0.0809021
#> Xnecon    -2.13686 1.751111 0.0459864      0.0696546
#> Xpresi     0.28542 0.291404 0.0076527      0.0102171
#> Xtag       0.04106 0.109635 0.0028791      0.0037634
#> Xrel       1.40353 0.628743 0.0165116      0.0210372
#> Xethn      1.14562 0.640030 0.0168080      0.0591352
#> Xprevdem   0.02349 0.293343 0.0077036      0.0089523
#> Xopenc    -0.01113 0.005499 0.0001444      0.0002753
#> 
#> 2. Quantiles for each variable:
#> 
#>               2.5%      25%      50%       75%     97.5%
#> X1        -0.07736  1.41928  2.30796  3.285902  7.801641
#> Xcalinv   -4.06635 -1.16363 -0.02407  1.114916  4.158685
#> Xlnlevel  -1.52077 -1.02394 -0.90893 -0.811584 -0.637564
#> Xcalileve -3.99743 -1.13507  0.08103  1.181787  4.188306
#> Xnecon    -6.13402 -3.18673 -1.97276 -0.962183  0.904850
#> Xpresi    -0.28510  0.10039  0.27996  0.480292  0.831018
#> Xtag      -0.18224 -0.03058  0.04745  0.118333  0.245317
#> Xrel       0.24135  0.96086  1.37781  1.808208  2.753970
#> Xethn     -0.12731  0.72600  1.15664  1.555914  2.436604
#> Xprevdem  -0.57760 -0.16734  0.03143  0.221004  0.572993
#> Xopenc    -0.02229 -0.01464 -0.01062 -0.007251 -0.001303
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-7-2.png" width="100%" />

#### Citation

To cite package`BayesMFSurv` in publications, please use:

``` r
citation(package = 'BayesMFSurv')
```
