

#' @title betas.slice.sampling2
#' @description slice sampling for betas
#' @param Sigma.b variance estimate of betas
#' @param Y the time (duration) dependent variable for the survival stage (t).
#' @param Y0 the elapsed time since inception until the beginning of time period (t-1).
#' @param X covariates for betas
#' @param betas current value of betas
#' @param alpha probability of true censoring
#' @param C censoring indicator
#' @param lambda current value of lambda
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param form type of parametric model (Exponential or Weibull)
#' @return One sample update using slice sampling
#' @export

betas.slice.sampling2 <- function(Sigma.b, Y, Y0, X, betas, alpha, C, lambda, w, m, form){

  p1 = length(betas)
  for (p in sample(1:p1, p1, replace = FALSE)){
    betas[p] = univ.betas.slice.sampling2(betas[p], p, Sigma.b, Y, Y0, X, betas, alpha, C, lambda, w, m, form = form)
  }
  return(betas)

}

#' @title univ.betas.slice.sampling2
#' @description univariate slice sampling for betas.p
#' @param betas.p current value of the pth element of betas
#' @param p pth element
#' @param Sigma.b variance estimate of betas
#' @param Y the time (duration) dependent variable for the survival stage (t).
#' @param Y0 the elapsed time since inception until the beginning of time period (t-1).
#' @param X covariates for betas
#' @param betas current value of betas
#' @param alpha probability of true censoring
#' @param C censoring indicator
#' @param lambda current value of lambda
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param lower lower bound on support of the distribution
#' @param upper upper bound on support of the distribution
#' @param form type of parametric model (Exponential or Weibull)
#' @return One sample update using slice sampling
#' @export
#'
univ.betas.slice.sampling2 <- function(betas.p, p, Sigma.b, Y, Y0, X, betas, alpha, C, lambda, w, m, lower = -Inf, upper = +Inf, form){

  b0 = betas.p
  b.post0 = betas.post2(b0, p, Sigma.b, Y, Y0, X, betas, alpha, C, lambda, form)
  if (exp(b.post0) > 0) { b.post0 = log(runif(1, 0, exp(b.post0)))}

  u = runif(1, 0, w)
  L = b0 - u
  R = b0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (betas.post2(L, p, Sigma.b, Y, Y0,X, betas, alpha, C, lambda, form) <= b.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (betas.post2(R, p, Sigma.b, Y, Y0, X, betas, alpha, C, lambda, form) <= b.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (betas.post2(L, p, Sigma.b, Y, Y0, X, betas, alpha, C, lambda, form) <= b.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (betas.post2(R, p, Sigma.b, Y, Y0, X, betas, alpha, C, lambda, form) <= b.post0) break
      R = R + w
      K = K - 1
    }
  }

  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }

  repeat
  {
    b1 = runif(1, L, R)
    b.post1 = betas.post2(b1, p, Sigma.b, Y,  Y0, X, betas, alpha, C, lambda, form)

    if (b.post1 >= b.post0) break
    if (b1 > b0) {
      R = b1
    } else {
      L = b1
    }
  }
  return(b1)

}


#' @title betas.post2
#' @description log-posterior distribution of betas with pth element fixed as betas.p
#' @param betas.p current value of the pth element of betas
#' @param p pth element
#' @param Sigma.b variance estimate of betas
#' @param Y the time (duration) dependent variable for the survival stage (t).
#' @param Y0 the elapsed time since inception until the beginning of time period (t-1).
#' @param X covariates for betas
#' @param betas current value of betas
#' @param alpha probability of true censoring
#' @param C censoring indicator
#' @param lambda current value of lambda
#' @param form type of parametric model (Exponential or Weibull)
#' @return log- posterior density of betas
#' @export

betas.post2 <- function(betas.p, p, Sigma.b, Y, Y0,  X, betas, alpha, C, lambda, form){

  betas[p] = betas.p
  if (form %in% "Weibull") {
    eXB = exp(X %*% betas)
  } else {
    eXB = exp(X %*% betas)
  }
  lprior = dmvnorm(betas, rep(0, length(betas)), Sigma.b, log = TRUE)
  lpost = llikWeibull2(Y, Y0, eXB, alpha, C, lambda) + lprior
  return(lpost)

}


#' @title gammas.post2
#' @description log-posterior distribution of gammas with pth element fixed as gammas.p
#' @param gammas.p current value of the pth element of gammas
#' @param p pth element
#' @param Sigma.g variance estimate of gammas
#' @param Y the time (duration) dependent variable for the survival stage (t).
#' @param Y0 the elapsed time since inception until the beginning of time period (t-1).
#' @param eXB exponentiated vector of covariates times betas
#' @param Z covariates for gammas
#' @param gammas current value of gammas
#' @param C censoring indicator
#' @param lambda current value of lambda
#' @param form type of parametric model (Exponential or Weibull)
#' @return log- posterior density of betas
#' @export

gammas.post2 <- function(gammas.p, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, form){

  gammas[p] = gammas.p
  if (form %in% "Weibull") {
    alpha = 1 / (1 + exp(-Z %*% gammas))
  } else {
    alpha = 1 / (1 + exp(-Z %*% gammas))
  }
  lprior = dmvnorm(gammas, rep(0, length(gammas)), Sigma.g, log = TRUE)
  lpost = llikWeibull2(Y, Y0, eXB, alpha, C, lambda) + lprior
  return(lpost)

}


#' @title gammas.slice.sampling2
#' @description slice sampling for gammas
#' @param Sigma.g variance estimate of gammas
#' @param Y the time (duration) dependent variable for the survival stage (t).
#' @param Y0 the elapsed time since inception until the beginning of time period (t-1).
#' @param eXB exponentiated vector of covariates times betas
#' @param Z covariates for gammas
#' @param gammas current value of gammas
#' @param C censoring indicator
#' @param lambda current value of lambda
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param form type of parametric model (Exponential or Weibull)
#' @return One sample update using slice sampling
#' @export

gammas.slice.sampling2 <- function(Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, w, m, form){

  p2 = length(gammas)
  for (p in sample(1:p2, p2, replace = FALSE)) {
    gammas[p] = univ.gammas.slice.sampling2(gammas[p], p, Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, w, m, form = form)
  }
  return(gammas)

}


#' @title univ.gammas.slice.sampling2
#' @description univariate slice sampling for gammas.p
#' @param gammas.p current value of the pth element of gammas
#' @param p pth element
#' @param Sigma.g variance estimate of gammas
#' @param Y the time (duration) dependent variable for the survival stage (t).
#' @param Y0 the elapsed time since inception until the beginning of time period (t-1).
#' @param eXB exponentiated vector of covariates times betas
#' @param Z covariates for gammas
#' @param gammas current value of gammas
#' @param C censoring indicator
#' @param lambda current value of lambda
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param lower lower bound on support of the distribution
#' @param upper upper bound on support of the distribution
#' @param form type of parametric model (Exponential or Weibull)
#' @return One sample update using slice sampling
#' @export

univ.gammas.slice.sampling2 <- function(gammas.p, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, w, m, lower = -Inf, upper = +Inf, form){

  g0 = gammas.p
  g.post0 = gammas.post2(g0, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, form)
  if (exp(g.post0) > 0) { g.post0 = log(runif(1, 0, exp(g.post0)))}

  u = runif(1, 0, w)
  L = g0 - u
  R = g0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (gammas.post2(L, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, form) <= g.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (gammas.post2(R, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, form) <= g.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (gammas.post2(L, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, form) <= g.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (gammas.post2(R, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, form) <= g.post0) break
      R = R + w
      K = K - 1
    }
  }

  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }

  repeat
  {
    g1 = runif(1, L, R)
    g.post1 = gammas.post2(g1, p, Sigma.g, Y, Y0, eXB, Z, gammas, C, lambda, form)

    if (g.post1 >= g.post0) break
    if (g1 > g0) {
      R = g1
    } else {
      L = g1
    }
  }
  return(g1)

}



#' @title lambda.post2
#' @description log-posterior distribution of lambda
#' @param Y the time (duration) dependent variable for the survival stage (t).
#' @param Y0 the elapsed time since inception until the beginning of time period (t-1).
#' @param eXB exponentiated vector of covariates times betas
#' @param alpha probability of true censoring
#' @param C censoring indicator
#' @param lambda current value of lambda
#' @param a shape parameter of gammas prior
#' @param b scale parameter of gammas prior
#' @return log- posterior density of betas
#' @export

lambda.post2 <- function(Y, Y0, eXB, alpha, C, lambda, a = 1, b = 1){

  lprior = dgamma(lambda, a, b, log = TRUE)
  lpost = llikWeibull2(Y, Y0, eXB, alpha, C, lambda) + lprior
  return(lpost)

}

#' @title lambda.slice.sampling2
#' @description univariate slice sampling for lambda
#' @param Y the time (duration) dependent variable for the survival stage (t).
#' @param Y0 the elapsed time since inception until the beginning of time period (t-1).
#' @param eXB exponentiated vector of covariates times betas
#' @param alpha probability of true censoring
#' @param C censoring indicator
#' @param lambda current value of lambda
#' @param w size of the slice in the slice sampling
#' @param m limit on steps in the slice sampling
#' @param lower lower bound on support of the distribution
#' @param upper upper bound on support of the distribution
#' @return One sample update using slice sampling
#' @export

lambda.slice.sampling2 <- function(Y, Y0, eXB, alpha, C, lambda, w, m, lower = 0.01, upper = +Inf){

  l0 = lambda
  l.post0 = lambda.post2(Y, Y0, eXB, alpha, C, l0)
  if (exp(l.post0) > 0) { l.post0 = log(runif(1, 0, exp(l.post0)))}

  u = runif(1, 0, w)
  L = l0 - u
  R = l0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (lambda.post2(Y, Y0, eXB, alpha, C, L) <= l.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (lambda.post2(Y, Y0, eXB, alpha, C, R) <= l.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (lambda.post2(Y, Y0, eXB, alpha, C, L) <= l.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (lambda.post2(Y, Y0, eXB, alpha, C, R) <= l.post0) break
      R = R + w
      K = K - 1
    }
  }

  if (L < lower) {
    L = lower
  }
  if (R > upper) {
    R = upper
  }

  repeat
  {
    l1 = runif(1, L, R)
    l.post1 = lambda.post2(Y, Y0, eXB, alpha, C, l1)
    if (l.post1 >= l.post0) break
    if (l1 > l0) {
      R = l1
    } else {
      L = l1
    }
  }
  return(l1)

}





