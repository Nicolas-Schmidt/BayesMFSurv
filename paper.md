---
title: 'BayesMFSurv: An R Package to Estimate Bayesian Split-Population Survival Models With (and Without) Misclassified Failure Events'
authors:
- affiliation: 1
  name: Minnie M. Joo
- affiliation: 2
  name: Nicolás Schmidt
  orcid: 0000-0001-5083-5792
- affiliation: 3
  name: Sergio Béjar
  orcid: 0000-0002-9352-3892
- affiliation: 4
  name: Vineeta Yadav
- affiliation: 4
  name: Bumba Mukherjee
date: "02 February 2020"
output:
  html_document:
    df_print: paged
  pdf_document: default
bibliography: paper.bib

tags:
- R
- survival analysis
- bayesian statistics

affiliations:
- index: 1
  name: Dept. of Political Science, University of Massachusetts Lowell
- index: 2
  name: Dept. of Political Science. Universidad de la Republica, UY
- index: 3
  name: Dept. of Political Science, San Jose State University
- index: 4
  name: Dept. of Political Science, Pennsylvania State University
---

# Summary

Social Scientists and Biostatisticians often employ conventional parametric survival or mixture cure models (e.g. Weibull, Exponential) to analyze outcome variables in survival data that focus on the time until an event occurred or "failed" [@box1999modeling; @maller1996survival; @lee17]. An important assumption underlying these models is that researchers record the date -year, month or day- in which an event or observation of interest failed (i.e. "terminated") _accurately_. Yet events that are recorded as having failed at a given point of time can be inaccurately measured [@clark2003survival; @schober2018survival; @bagozzi2019bayesian]. Inaccurate measurement of this sort leads to a subset of _misclassified failure_ cases in survival data in which some subjects are recorded as having failed or experienced the event of interest even though they in actuality "live on" past their recorded-failure point. 

There are several scenarios where a subset of recorded failure events may persist beyond their recorded failure time, leading to misclassiﬁcation in event failures. For example, political scientists who analyze the duration of civil wars fought between rebel groups and governments often record end dates (“failures”) for specific conﬂicts based upon 24-month spells with fewer than 25 battle-deaths per year [@balch2000killing; @thyne2012information]. The aforementioned threshold is prone to measurement error, especially for lower-intensity civil wars in poor information environments that persist beyond their recorded end date. Other examples include the study of the duration of ancient civilizations [@cioffi1999evolution] and the time taken to detect cancer [@schober2018survival]. In both these latter examples, researchers typically do not have data on the precise time-point of a given failure due to the sands of time or because of inaccurate information. This leads them to —similar to the civil conflict example— underestimate the duration of particularly misclassified event failure cases.

Since these underestimates of duration are non-random, bias will arise in survival estimates of the phenomena mentioned above when using conventional survival models. Hence, the main motivation for developing the Bayesian Misclassiﬁed Failure (MF hereafter) split population survival model is to resolve methodological challenges resulting from misclassiﬁed event failures by accounting for the possibility that some failure events survive beyond their recorded failure time. The development of the Bayesian MF model is also driven by the fact that it permits researchers to identify when the end date of observations in survival data is misclassiﬁed therein providing substantive insights into this process. Further, there does not exist an R package that extracts posterior distribution of estimates from parametric cure (split-population) models, including the MF model, using Bayesian Markov Chain Monte Carlo (MCMC) methods. `BayesMFSurv` [@BayesMFSurv-v1] is an R package [@r-core] that contains functions and computationally intensive routines in `C++` to fit the parametric Weibull and Exponential (i) survival model and (ii) Misclassified Failure survival model via Bayesian MCMC methods using slice-sampling [@neal2003slice; @bagozzi2019bayesian].

# Motivation, Description, Applications

Numerous R packages offer functionalities to estimate conventional parametric and semi-parametric survival models via maximum likelihood estimation (MLE) or Bayesian MCMC methods [@OIsurv; @Zhouetal-spBayes; @survival-package; @wangetal-dynsurv]. Other R packages focus on estimation of parametric or semi-parametric cure survival models using MLE [@smcure-cai; @beger2017splitting; @rcure-r; @amdahl-flex]. To our knowledge, there does not exist an R package that fits parametric mixture cure models, including the MF survival model, via Bayesian MCMC (e.g. slice sampling) methods which offer a powerful yet flexible tool for estimating such models. Further, extant R packages that use Bayesian inference for survival analyses only focus on standard survival models that do not take into account latent misclassified failure events in survival data. Because misclassified failure events are right-censored events, there is a non-zero probability that these misclassified cases persisted beyond their recorded failure time. Failing to account for misclassified failure events in survival data that results from estimating standard survival or cure models will lead researchers to underestimate the duration of time of these events.

Since the underestimates of duration are non-random, bias will arise in survival estimates of these phenomena when researchers use standard survival or cure models. To address this misclassified failure challenge in survival data, our `BayesMFSurv` R package incorporates various functions listed below that fit @bagozzi2019bayesian's parametric MF survival model via Bayesian MCMC methods. This model estimates a system of two equations to account for the possibility that some unknown subset of failure events actually "lived on" beyond their recorded failure time. The first is a "splitting" equation that estimates the probability of a case being a misclassified failure, with or without covariates. The second equation is a standard parametric survival model, whose relevant failure and survival probabilities are estimated conditional on a case _not_ being a misclassified failure. These features of the model in `BayesMFSurv` account for a heterogeneous mixture of failure cases in survival data and address the non-random underestimates of duration for misclassified failure events. `BayesMFSurv` also incorporates time-varying covariates that are common in panel survival datasets. This model can be applied at least to the following survival datasets where misclassified failure cases are prevalent: civil war termination that determines civil war duration [@thyne2012information], time taken to detect onset of cancer [@schober2018survival], and collapse (and thus duration) of ancient civilizations or political regimes [@cioffi1999evolution; @reenock2007regressive].

# BayesMFSurv R Package

The R package `BayesMFSurv` contains four functions to fit the parametric (Weibull and Exponential) (i) standard survival model and (ii) MF survival model via Bayesian MCMC using a slice-sampling algorithm described in @bagozzi2019bayesian. Bayesian MCMC estimation is conducted by using the Multivariate Normal prior for these models' split and survival stage parameters, and the Gamma prior for the shape parameter. The four functions in `BayesMFSurv` are:   

* `mfsurv`: Fits a parametric MF model via Bayesian MCMC with slice-sampling to estimate the misclassification failure probability in the split (first) stage and hazard in the second (survival) stage. Slice-sampling, which is conducted by using the univariate slice sampler [@neal2003slice], is employed to draw the posterior sample of the model's split and survival stage parameters.
* `mcmcsurv`: Fits a standard parametric survival model via Bayesian MCMC with slice-sampling employed to draw the posterior sample of the model's survival stage parameters.
* `stats`: Calculates log-likelihood and deviance information criterion (DIC) for fitted model objects of class `mfsurv` _and_ `mcmcsurv`.
* `summary`: Summarizes Bayesian MCMC output -this includes the mean, standard deviation, standard error of the mean, and 95 credible interval- of each parameter's posterior distribution from the Bayesian standard and MF parametric survival models.

The ease and speed of estimating the Bayesian standard and MF parametric survival models in `BayesMFSurv` is improved by using `C++` to perform computationally intensive routines (e.g. slice-sampling) before pulling the output into R. Users can also illustrate trace-plots and kernel density plots for each parameter from `mcmcsurv` and `mfsurv` that fits the Bayesian standard and MF parametric models respectively. To illustrate the `BayesMFSurv` package's functionality, all the four functions listed above have been tested on a survival dataset of democratic regime failure [@reenock2007regressive] described and included in this package.

# Availability

`BayesMFSurv` is an open source software made available under the MIT license. It can be installed from its github repository using the `remotes` package: `remotes::install_github("Nicolas-Schmidt/BayesMFSurv")`.

# References
