#' @title mcmcsurv
#' @description Markov Chain Monte Carlo (MCMC) to run Bayesian survival (Weibull) model
#'
#' @param Y response variable
#' @param Y0 the elapsed time since inception until the beginning of time period (t-1).
#' @param C censoring indicator
#' @param X covariates for betas
#' @param N number of MCMC iterations
#' @param burn burn-in to be discarded
#' @param thin thinning to prevent from autocorrelation
#' @param w size of the slice in the slice sampling for (betas, gammas, rho)
#' @param m limit on steps in the slice sampling
#' @param form type of parametric model (Exponential or Weibull)
#'
#' @return chain of the variables of interest
#'
#' @export
mcmcsurv <- function(Y, Y0,C, X, N, burn, thin, w = c(1, 1, 1), m = 10, form) {
    p1 = dim(X)[2]
    # initial values
    betas = rep(0, p1)
    rho = 1
    W = rep(0, length(Y))
    delta = rep(0, length(Y))
    Sigma.b = 10 *p1  * diag(p1)
    betas.samp = matrix(NA, nrow = (N - burn) / thin, ncol = p1)
    rho.samp = rep(NA, (N - burn) / thin)
    for (iter in 1:N) {
        if (iter %% 5000 == 0) print(iter)
        if (iter > burn) {
            Sigma.b = riwish(1 + p1, betas %*% t(betas) + p1* diag(p1))
        }
        betas = betas.slice.sampling(Sigma.b, Y, Y0,X, W, betas, delta, C, rho, w[1], m, form = form)
        eXB = exp(X %*% betas+ W)
        if (form %in% "Weibull") {
            rho = rho.slice.sampling(Y,Y0, eXB, delta, C, rho, w[3], m)
        }
        if (iter > burn & (iter - burn) %% thin == 0) {
            betas.samp[(iter - burn) / thin, ] = betas
            rho.samp[(iter - burn) / thin] = rho
        }
    }

    colnames(betas.samp) <- paste0("X", colnames(X))#---------------------------
    out <- list(betas = betas.samp, rho = rho.samp, Y=Y, Y0=Y0, X=X, N=N, C=C, iterations = N, burn_in = burn,
                thinning = thin, betan = nrow(betas), distribution = form)
    class(out) <- "mcmcsurv"
    return(out)
}



#' @title summary.mcmcsurv
#' @description Returns a summary of a mfsurv object via \code{\link[coda]{summary.mcmc}}.
#' @param object an object of class \code{mfsurv}, the output of \code{\link{mfsurv}}.
#' @param parameter one of three parameters of the mfsurv output. Indicate either "betas" or "rho".
#' @param ... additional parameter
#' @return list. Empirical mean, standard deviation and quantiles for each variable.
#' @rdname mcmcsurv
#' @export

summary.mcmcsurv <- function(object, parameter = c("betas", "rho"), ...){

    if (parameter == "betas"){
        sum <- summary(mcmc(object$betas))
        return(sum)
    }
    if (parameter == "rho"){
        sum <- summary(mcmc(object$lambda))
        return(sum)
    }
}

