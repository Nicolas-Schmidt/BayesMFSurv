
#' @title mfsurv.stats
#' @description A function to calculate the deviance information criterion (DIC) for fitted model objects of class \code{mfsurv}
#' for which a log-likelihood can be obtained, according to the formula \emph{DIC = -2 * (L - P)},
#' where \emph{L} is the log likelihood of the data given the posterior means of the parameter and
#' \emph{P} is the  estimate of the effective number of parameters in the model.
#' @param object an object of class \code{mfsurv}, the output of \code{mfsurv()}.
#' @return list.
#' @export

mfsurv.stats <- function(object){
    .Deprecated("stats")
    return(stats(object = object))
}


#' @title summary.mfsurv
#' @description Returns a summary of a mfsurv object via \code{\link[coda]{summary.mcmc}}.
#' @param object an object of class \code{mfsurv}, the output of \code{\link{mfsurv}}.
#' @param parameter one of three parameters of the mfsurv output. Indicate either "betas", "gammas" or "lambda".
#' @param ... additional parameter
#' @return list. Empirical mean, standard deviation and quantiles for each variable.
#' @rdname mfsurv
#' @export

mfsurv.summary <- function(object, parameter){
    .Deprecated("summary")
    return(summary(object = object, parameter = parameter))
}

