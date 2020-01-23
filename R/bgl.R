
#' @importFrom stats dgamma runif rgamma dnorm
# @import grDevices
# @import graphics
# @import RcppArmadillo
# @importFrom Rcpp sourceCpp
#' @importFrom MCMCpack riwish
#' @importFrom coda mcmc
#' @importFrom FastGP rcpp_rmvnorm rcpp_log_dmvnorm
#NULL


# @title betas.slice.sampling
# @description slice sampling for betas
#
# @param Sigma.b variance estimate of betas
# @param Y response variable
# @param Y0 add...
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Exponential or Weibull)
#
# @return One sample update using slice sampling
#
# @export
betas.slice.sampling = function(Sigma.b, Y, Y0, X, W, betas, delta, C, rho, w, m, form) {
  p1 = length(betas)
  for (p in sample(1:p1, p1, replace = FALSE)) {
    betas[p] = univ.betas.slice.sampling(betas[p], p, Sigma.b, Y, Y0,X, W, betas, delta, C, rho, w, m, form = form)
  }
  return(betas)
}

# @title univ.betas.slice.sampling
# @description univariate slice sampling for betas.p
#
# @param betas.p current value of the pth element of betas
# @param p pth element
# @param Sigma.b variance estimate of betas
# @param Y response variable
# @param Y0 add...
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
# @param form type of parametric model (Exponential or Weibull)
#
# @return One sample update using slice sampling
#
# @export
univ.betas.slice.sampling = function(betas.p, p, Sigma.b, Y,Y0, X, W, betas, delta, C, rho, w, m, lower = -Inf, upper = +Inf, form) {
  b0 = betas.p
  b.post0 = betas.post(b0, p, Sigma.b, Y, Y0, X, W, betas, delta, C, rho, form)
  if (exp(b.post0) > 0) { b.post0 = log(runif(1, 0, exp(b.post0)))}

  u = runif(1, 0, w)
  L = b0 - u
  R = b0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (betas.post(L, p, Sigma.b, Y, Y0, X, W, betas, delta, C, rho, form) <= b.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (betas.post(R, p, Sigma.b, Y, Y0,X, W, betas, delta, C, rho, form) <= b.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (betas.post(L, p, Sigma.b, Y,Y0, X, W, betas, delta, C, rho, form) <= b.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (betas.post(R, p, Sigma.b, Y,Y0, X, W, betas, delta, C, rho, form) <= b.post0) break
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
    b.post1 = betas.post(b1, p, Sigma.b, Y, Y0,X, W, betas, delta, C, rho, form)

    if (b.post1 >= b.post0) break
    if (b1 > b0) {
      R = b1
    } else {
      L = b1
    }
  }
  return(b1)
}

# @title gammas.slice.sampling
# @description slice sampling for gammas
#
# @param Sigma.g variance estimate of gammas
# @param Y response variable
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Exponential or Weibull)
#
# @return One sample update using slice sampling
#
# @export
gammas.slice.sampling = function(Sigma.g, Y, eXB, Z, gammas, C, rho, w, m, form) {
  p2 = length(gammas)
  for (p in sample(1:p2, p2, replace = FALSE)) {
    gammas[p] = univ.gammas.slice.sampling(gammas[p], p, Sigma.g, Y, eXB, Z, gammas, C, rho, w, m, form = form)
  }
  return(gammas)
}

# @title gammas.slice.sampling_2
# @description slice sampling for gammas
#
# @param Sigma.g variance estimate of gammas
# @param Y response variable
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param form type of parametric model (Exponential or Weibull)
#
# @return One sample update using slice sampling
#
# @export
gammas.slice.sampling_2 = function(Sigma.g, Y, eXB, Z, V, gammas, C, rho, w, m, form) {
  p2 = length(gammas)
  for (p in sample(1:p2, p2, replace = FALSE)) {
    gammas[p] = univ.gammas.slice.sampling_2(gammas[p], p, Sigma.g, Y, eXB, Z, V, gammas, C, rho, w, m, form = form)
  }
  return(gammas)
}
# @title univ.gammas.slice.sampling
# @description univariate slice sampling for gammas.p
#
# @param gammas.p current value of the pth element of gammas
# @param p pth element
# @param Sigma.g variance estimate of gammas
# @param Y response variable
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
# @param form type of parametric model (Exponential or Weibull)
#
# @return One sample update using slice sampling
#
# @export
univ.gammas.slice.sampling = function(gammas.p, p, Sigma.g, Y, eXB, Z, gammas, C, rho, w, m, lower = -Inf, upper = +Inf, form) {
  g0 = gammas.p
  g.post0 = gammas.post(g0, p, Sigma.g, Y, eXB, Z,gammas, C, rho, form)
  if (exp(g.post0) > 0) { g.post0 = log(runif(1, 0, exp(g.post0)))}


  u = runif(1, 0, w)
  L = g0 - u
  R = g0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (gammas.post(L, p, Sigma.g, Y, eXB, Z, gammas, C, rho, form) <= g.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (gammas.post(R, p, Sigma.g, Y, eXB, Z, gammas, C, rho, form) <= g.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (gammas.post(L, p, Sigma.g, Y, eXB, Z, gammas, C, rho, form) <= g.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (gammas.post(R, p, Sigma.g, Y, eXB, Z, gammas, C, rho, form) <= g.post0) break
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
    g.post1 = gammas.post(g1, p, Sigma.g, Y, eXB, Z, gammas, C, rho, form)

    if (g.post1 >= g.post0) break
    if (g1 > g0) {
      R = g1
    } else {
      L = g1
    }
  }
  return(g1)
}

# @title univ.gammas.slice.sampling_2
# @description univariate slice sampling for gammas.p
#
# @param gammas.p current value of the pth element of gammas
# @param p pth element
# @param Sigma.g variance estimate of gammas
# @param Y response variable
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
# @param form type of parametric model (Exponential or Weibull)
#
# @return One sample update using slice sampling
#
# @export
univ.gammas.slice.sampling_2 = function(gammas.p, p, Sigma.g, Y, eXB, Z, V, gammas, C, rho, w, m, lower = -Inf, upper = +Inf, form) {
  g0 = gammas.p
  g.post0 = gammas.post_2(g0, p, Sigma.g, Y, eXB, Z, V, gammas, C, rho, form)
  if (exp(g.post0) > 0) { g.post0 = log(runif(1, 0, exp(g.post0)))}


  u = runif(1, 0, w)
  L = g0 - u
  R = g0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (gammas.post_2(L, p, Sigma.g, Y, eXB, Z, V, gammas, C, rho, form) <= g.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (gammas.post_2(R, p, Sigma.g, Y, eXB, Z, V, gammas, C, rho, form) <= g.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (gammas.post_2(L, p, Sigma.g, Y, eXB, Z, V, gammas, C, rho, form) <= g.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (gammas.post_2(R, p, Sigma.g, Y, eXB, Z, V, gammas, C, rho, form) <= g.post0) break
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
    g.post1 = gammas.post_2(g1, p, Sigma.g, Y, eXB, Z, V, gammas, C, rho, form)

    if (g.post1 >= g.post0) break
    if (g1 > g0) {
      R = g1
    } else {
      L = g1
    }
  }
  return(g1)
}

# @title rho.slice.sampling
# @description univariate slice sampling for rho
#
# @param Y response variable
# @param Y0 add...
# @param eXB exponentiated vector of covariates times betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param w size of the slice in the slice sampling
# @param m limit on steps in the slice sampling
# @param lower lower bound on support of the distribution
# @param upper upper bound on support of the distribution
#
# @return One sample update using slice sampling
#
# @export
rho.slice.sampling = function(Y,Y0, eXB, delta, C, rho, w, m, lower = 0.01, upper = +Inf) {
  l0 = rho
  l.post0 = rho.post(Y, Y0,eXB, delta, C, l0)
  if (exp(l.post0) > 0) { l.post0 = log(runif(1, 0, exp(l.post0)))}


  u = runif(1, 0, w)
  L = l0 - u
  R = l0 + (w - u)
  if (is.infinite(m)) {
    repeat
    { if (L <= lower) break
      if (rho.post(Y, Y0,eXB, delta, C, L) <= l.post0) break
      L = L - w
    }
    repeat
    {
      if (R >= upper) break
      if (is.na(rho.post(Y, Y0,eXB, delta, C, R))) browser()
      if (rho.post(Y, Y0,eXB, delta, C, R) <= l.post0) break
      R = R + w
    }
  } else if (m > 1) {
    J = floor(runif(1, 0, m))
    K = (m - 1) - J

    while (J > 0) {
      if (L <= lower) break
      if (rho.post(Y, Y0,eXB, delta, C, L) <= l.post0) break
      L = L - w
      J = J - 1
    }

    while (K > 0) {
      if (R >= upper) break
      if (is.na(rho.post(Y, Y0,eXB, delta, C, R))) browser()
      if (rho.post(Y, Y0,eXB, delta, C, R) <= l.post0) break
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
    l.post1 = rho.post(Y, Y0,eXB, delta, C, l1)

    if (l.post1 >= l.post0) break
    if (l1 > l0) {
      R = l1
    } else {
      L = l1
    }
  }
  return(l1)
}

# @title betas.post
# @description log-posterior distribution of betas with pth element fixed as betas.p
#
# @param betas.p current value of the pth element of betas
# @param p pth element
# @param Sigma.b variance estimate of betas
# @param Y response variable
# @param Y0 add...
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param form type of parametric model (Exponential or Weibull)
#
# @return log- posterior density of betas
#
# @export
betas.post = function(betas.p, p, Sigma.b, Y, Y0,X, W, betas, delta, C, rho, form) {
  betas[p] = betas.p
  eXB = exp(X %*% betas + W)
  lprior = rcpp_log_dmvnorm(Sigma.b, rep(0, length(betas)), betas, FALSE)
  lpost = llikWeibull3(Y,Y0, eXB, delta, C, rho) + lprior
  return(lpost)
}

# @title gammas.post
# @description log-posterior distribution of gammas with pth element fixed as gammas.p
#
# @param gammas.p current value of the pth element of gammas
# @param p pth element
# @param Sigma.g variance estimate of gammas
# @param Y response variable
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param form type of parametric model (Exponential or Weibull)
#
# @return log- posterior density of betas
#
# @export
gammas.post = function(gammas.p, p, Sigma.g, Y, eXB, Z, gammas, C, rho, form) {
  gammas[p] = gammas.p
  delta = 1 / (1 + exp(-Z %*% gammas))
  lprior = rcpp_log_dmvnorm(Sigma.g, rep(0, length(gammas)), gammas, FALSE)
  lpost = llikWeibull3(Y, eXB, delta, C, rho) + lprior
  return(lpost)
}

# @title gammas.post_2
# @description log-posterior distribution of gammas with pth element fixed as gammas.p
#
# @param gammas.p current value of the pth element of gammas
# @param p pth element
# @param Sigma.g variance estimate of gammas
# @param Y response variable
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param form type of parametric model (Exponential or Weibull)
#
# @return log- posterior density of betas
#
# @export
gammas.post_2 = function(gammas.p, p, Sigma.g, Y, eXB, Z, V, gammas, C, rho, form) {
  gammas[p] = gammas.p
  delta = 1 / (1 + exp(-Z %*% gammas- V))
  lprior = rcpp_log_dmvnorm(Sigma.g, rep(0, length(gammas)), gammas, FALSE)
  lpost = llikWeibull3(Y, eXB, delta, C, rho) + lprior
  return(lpost)
}


# @title rho.post
# @description log-posterior distribution of rho
#
# @param Y response variable
# @param Y0 add...
# @param eXB exponentiated vector of covariates times betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param a shape parameter of gammas prior
# @param b scale parameter of gammas prior
#
# @return log- posterior density of betas
#
# @export
rho.post = function(Y, Y0,eXB, delta, C, rho, a = 1, b = 1) {
  lprior = dgamma(rho, a, b, log = TRUE)
  lpost = llikWeibull3(Y,Y0, eXB, delta, C, rho) + lprior
  return(lpost)
}


# @title lambda.gibbs.sampling
# @description log-posterior distribution of rho
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param W spatial random effects
# @param a shape parameter of gammas prior
# @param b scale parameter of gammas prior
#
# @return log- posterior density of betas
#
# @export
lambda.gibbs.sampling <- function(S, A, W, a = 1, b = 1) {
  S_uniq = unique(cbind(S, W))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  J = nrow(S_uniq)
  sums = 0
  for (j in 1:J) {
    adj = which(A[j,]==1)
    m_j = length(adj)
    W_j_bar = mean(S_uniq[which(S_uniq[,1] %in% adj),2])
    W_j = S_uniq[j, 2]
    sums = sums + m_j/2 * (W_j-W_j_bar)^2
  }
  lambda = rgamma(1, J + a, sums + b)
  return(lambda)
}


# @title W.post
# @description log-posterior distribution of W with sth element fixed as W.s
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y response variable
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
#
# @return log- posterior density of betas
#
# @export
W.post = function(S, A, lambda, Y,X, W, betas, delta, C, rho) {
  eXB = exp(X %*% betas + W)
  S_uniq = unique(cbind(S, W))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  lprior = 0
  for (s in 1:nrow(A)) {
  adj = which(A[s,] == 1)
  m_j = length(adj)
  W_j_bar = mean(S_uniq[which(S_uniq[,1] %in% adj),2])
  lprior = lprior + dnorm(S_uniq[s, 2], W_j_bar, sqrt(1/(lambda * m_j)), log = TRUE)
  }
  lpost = llikWeibull3(Y, eXB, delta, C, rho) + lprior
  return(lpost)
}

# @title W.MH.sampling
# @description slice sampling for W
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y response variable
# @param X covariates for betas
# @param W spatial random effects
# @param betas current value of betas
# @param delta probability of true censoring
# @param C censoring indicator
# @param rho current value of rho
# @param prop.var proposal variance for Metropolis-Hastings
#
# @return One sample update using slice sampling
#
# @export
W.MH.sampling = function(S, A, lambda, Y, X, W, betas, delta, C, rho, prop.var) {
  S_uniq = unique(cbind(S, W))
  W_old = S_uniq[order(S_uniq[,1]), 2]
  W_new = rcpp_rmvnorm(1, prop.var * diag(length(W_old)), W_old)
  W_new = W_new - mean(W_new)
  u = log(runif(1))
  temp = W.post(S, A, lambda, Y, X, W_new[S], betas, delta, C, rho) -
  		 W.post(S, A, lambda, Y, X, W_old[S], betas, delta, C, rho)
  alpha = min(0, temp)
  if (u <= alpha) {
  	W = W_new[S]
  } else {
  	W = W_old[S]
  }
  return(W)
}

# @title V.post
# @description log-posterior distribution of W with sth element fixed as W.s
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y response variable
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
#
# @return log- posterior density of betas
#
# @export
V.post = function(S, A, lambda, Y, eXB, Z, V, gammas, C, rho) {
  delta = 1 / (1 + exp(- Z %*% gammas - V))
  S_uniq = unique(cbind(S, V))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  lprior = 0
  for (s in 1:nrow(A)) {
  adj = which(A[s,] == 1)
  m_j = length(adj)
  V_j_bar = mean(S_uniq[which(S_uniq[,1] %in% adj),2])
  lprior = lprior + dnorm(S_uniq[s, 2], V_j_bar, sqrt(1/(lambda * m_j)), log = TRUE)
  }
  lpost = llikWeibull3(Y, eXB, delta, C, rho) + lprior
  return(lpost)
}

# @title V.MH.sampling
# @description slice sampling for W
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param lambda CAR parameter
# @param Y response variable
# @param eXB exponentiated vector of covariates times betas
# @param Z covariates for gammas
# @param V spatial random effects
# @param gammas current value of gammas
# @param C censoring indicator
# @param rho current value of rho
# @param prop.var proposal variance for Metropolis-Hastings
#
# @return One sample update using slice sampling
#
# @export
V.MH.sampling = function(S, A, lambda, Y, eXB, Z, V, gammas, C, rho, prop.var) {
  S_uniq = unique(cbind(S, V))
  V_old = S_uniq[order(S_uniq[,1]), 2]
  V_new = rcpp_rmvnorm(1, prop.var * diag(length(V_old)), V_old)
  V_new = V_new - mean(V_new)
  u = log(runif(1))
  temp = V.post(S, A, lambda, Y, eXB, Z, V_new[S], gammas, C, rho) -
  		 V.post(S, A, lambda, Y, eXB, Z, V_old[S], gammas, C, rho)
  alpha = min(0, temp)
  if (u <= alpha) {
  	V = V_new[S]
  } else {
  	V = V_old[S]
  }
  return(V)
}

# @title lambda.gibbs.sampling2
# @description log-posterior distribution of rho
#
# @param S spatial information (e.g. district)
# @param A adjacency information corresponding to spatial information
# @param W spatial random effects
# @param V spatial random effects
# @param a shape parameter of gammas prior
# @param b scale parameter of gammas prior
#
# @return log- posterior density of betas
#
# @export
lambda.gibbs.sampling2 <- function(S, A, W, V, a = 1, b = 1) {
  S_uniq = unique(cbind(S, W, V))
  S_uniq = S_uniq[order(S_uniq[,1]),]
  J = nrow(S_uniq)
  sums = 0
  for (j in 1:J) {
    adj = which(A[j,]==1)
    m_j = length(adj)
    W_j_bar = mean(S_uniq[which(S_uniq[,1] %in% adj),2])
    V_j_bar = mean(S_uniq[which(S_uniq[,1] %in% adj),3])
    W_j = S_uniq[j, 2]
    V_j = S_uniq[j, 3]
    sums = sums + m_j/2 * ((W_j-W_j_bar)^2 + (V_j-V_j_bar)^2)
  }
  lambda = rgamma(1, J + a, sums + b)
  return(lambda)
}



