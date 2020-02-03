---
title: 'BayesMFSurv: An R Package to Estimate Bayesian Split-Population Survival Models With (and Without) Misclassified Failure Events'
authors:
- affiliation: 1
  name: Minnie M. Joo
- affiliation: 2
  name: Nicolas Schmidt
  orcid: 0000-0001-5083-5792
- affiliation: 3
  name: Sergio Bejar
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

Social Scientists and Biostatisticians often employ conventional parametric survival or mixture cure models (e.g. Weibull, Exponential) to analyze outcome variables in survival data that focus on the time until an event occurred or "failed" [@box1999modeling; @maller1996survival; @lee17]. An important assumption underlying these models is that researchers record the date -year, month or day- in which an event or observation of interest failed (i.e. "terminated") _accurately_. Yet events that are recorded as having failed at a given point of time can be inaccurately measured [@clark2003survival; @schober2018survival; @bagozzi2019bayesian]. Inaccurate measurement of this sort leads to a subset of _misclassified failure_ cases in survival data in which some subjects are recorded as having failed or experienced the event of interest even though they in actuality "live on" past their recorded-failure point. Conventional survival and cure models yield biased parameter estimates in survival datasets that include latent misclassified failure cases as they do not statistically account for such cases [@schober2018survival; @bagozzi2019bayesian].

As described below, [@bagozzi2019bayesian] develop a Bayesian parametric Misclassified Failure (MF) survival model which -similar to standard cure models [@zhou2017spbayessurv; @box1999modeling; @maller1996survival]- consists of a system of two equations that addresses misclassified event failures in survival data. `BayesMFSurv` [@BayesMFSurv-v1] is an R package [@r-core] that contains functions and computationally intensive routines in `C++` to fit the parametric Weibull and Exponential (i) survival model and (ii) Misclassified Failure survival model via Bayesian Markov Chain Monte Carlo (MCMC) methods using slice-sampling [@neal2003slice; bagozzi2019bayesian].

# Motivation, Description, Applications

Numerous R packages offer functionalities to estimate conventional parametric and semi-parametric survival models via maximum likelihood estimation (MLE) or Bayesian MCMC methods [@OIsurv; @Zhouetal-spBayes; @survival-package; @wangetal-dynsurv]. Other R packages focus on estimation of parametric or semi-parametric cure survival models using MLE [@smcure-cai; @beger2017splitting; @rcure-r; @amdahl-flex]. To our knowledge, there does not exist an R package that fits parametric mixture cure models, including the MF survival model, via _BayesianMCMC_ (e.g. slice sampling) methods which offer a powerful yet flexible tool for estimating such models. Further, extant R packages that use Bayesian inference for survival analyses only focus on standard survival models that do not take into account latent misclassified failure events in survival data. Because misclassified failure events are right-censored events, there is a non-zero probability that these misclassified cases persisted beyond their recorded failure time. Failing to account for misclassified failure events in survival data that results from estimating standard survival or cure models will lead researchers to underestimate the duration of time of these events.

Since the underestimates of duration are non-random, bias will arise in survival estimates of these phenomena when researchers use standard survival or cure models. To address this misclassified failure challenge in survival data, our `BayesMFSurv` R package incorporates various functions listed below that fit [@bagozzi2019bayesian] parametric Misclassified Failure (MF) survival model via Bayesian MCMC methods. This model estimates a system of two equations to account for the possibility that some unknown subset of failure events actually "lived on" beyond their recorded failure time. The first is a "splitting" equation that estimates the probability of a case being a misclassified failure, with or without covariates. The second equation is a standard parametric survival model, whose relevant failure and survival probabilities are estimated conditional on a case _not_ being a misclassified failure. These features of the MF survival model in `BayesMFSurv` accounts for a heterogeneous mixture of failure cases in survival data and address the non-random underestimates of duration for misclassified failure events. The MF survival model in `BayesMFSurv` also incorporates time-varying covariates that are common in panel survival datasets and this model can be applied to the following survival datasets where misclassified failure cases are prevalent: civil war termination that determines civil war duration [@thyne2012information], time taken to detect onset of cancer [@schober2018survival], and collapse (and thus duration) of ancient civilizations or political regimes [@cioffi1999evolution; @reenock2007regressive].

# BayesMFSurv R Package

The R package `BayesMFSurv` contains four functions to fit the parametric (Weibull and Exponential) (i) standard survival model and (ii) Misclassified Failure (MF) survival model via Bayesian MCMC using a slice-sampling algorithm described in [@bagozzi2019bayesian]. Bayesian MCMC estimation is conducted by using the Multivariate Normal prior for these models' split and survival stage parameters, and the Gamma prior for the shape parameter. The four functions in `BayesMFSurv`are:   


*`mfsurv`: Fits a parametric MF model via Bayesian MCMC with slice-sampling to estimate the misclassification failure probability in the split (first) stage and hazard in the second (survival) stage. Slice-sampling, which is conducted by using the univariate slice sampler [@neal2003slice], is employed to draw the posterior sample of the model's split and survival stage parameters.
*`mcmcsurv`: Fits a standard parametric survival model via Bayesian MCMC with slice-sampling employed to draw the posterior sample of the model's survival stage parameters.
*`stats`: Calculates log-likelihood and deviance information criterion (DIC) for fitted model objects of class `mfsurv` _and_ `mcmcsurv`.
*`summary`: Summarizes Bayesian MCMC output -this includes the mean, standard deviation, standard error of the mean, and 95 credible interval- of each parameter's posterior distribution from the Bayesian standard and MF parametric survival models.


The ease and speed of estimating the Bayesian standard and MF parametric survival models in `BayesMFSurv` is improved by using `C++` to perform computationally intensive routines (e.g. slice-sampling) before pulling the output into R. Users can also illustrate trace-plots and kernel density plots for each parameter from `mcmcsurv` and `mfsurv` that fits the Bayesian standard and MF parametric models respectively. To illustrate the `BayesMFSurv` package's functionality, all the four functions listed above has been tested on a survival dataset of democratic regime failure [@reenock2007regressive] described and included in this package.

# Availability

`BayesMFSurv` is an open source software made available under the MIT license. It can be installed from its github repository using the `remotes` package: `remotes::install_git(remotes::install_github("Nicolas-Schmidt/BayesMFSurv")`.

# References
