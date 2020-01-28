#' @title mfsurv.stats
#' @description A function to calculate the deviance information criterion (DIC) for fitted model objects of class \code{mfsurv}
#' for which a log-likelihood can be obtained, according to the formula \emph{DIC = -2 * (L - P)},
#' where \emph{L} is the log likelihood of the data given the posterior means of the parameter and
#' \emph{P} is the  estimate of the effective number of parameters in the model.
#' @param object an object of class \code{mfsurv}, the output of \code{mfsurv()}.
#' @return list.
#' @examples
#' set.seed(95)
#' bgl <- Buhaugetal_2009_JCR
#' bgl <- subset(bgl, coupx == 0)
#' bgl <- na.omit(bgl)
#' Y   <- bgl$Y
#' X   <- as.matrix(cbind(1, bgl[,1:7]))
#' C   <- bgl$C
#' Z1  <- matrix(1, nrow = nrow(bgl))
#' Y0  <- bgl$Y0
#' model1 <- mfsurv(Y ~ X | C ~ Z1, Y0 = Y0,
#'                 N = 50,
#'                 burn = 20,
#'                 thin = 15,
#'                 w = c(0.1, .1, .1),
#'                 m = 5,
#'                 form = "Weibull",
#'                 na.action = 'na.omit')
#'
#' stats(model1)
#' @export

stats <- function(object){

    class.object <- c("mfsurv","mcmcsurv")
    if(!inherits(object, class.object)){stop("'object' must be of the class 'mfsurv' or 'mcmcsurv'.")}
    if(class(object) == class.object[1]){  #Calculate L
        X <- as.matrix(object$X)
        Z <- as.matrix(object$Z)
        Y <- as.matrix(object$Y)
        Y0 <- as.matrix(object$Y0)
        C <- as.matrix(object$C)
        data <- as.data.frame(cbind(object$Y, object$Y0, object$C, object$X, object$Z))
        theta_post = cbind(object$gammas, object$betas, object$lambda)
        theta_hat = apply(theta_post, 2, mean)
        L = llFun(theta_hat,Y, Y0,C,X,Z,data)$llik
        #Calculate P
        S = nrow(theta_post) #S = number of iterations
        #Add up the log likelihoods of each iteration
        llSum = 0
        sum1 = 0
        sum2 = 0
        sum3 = 0
        #l <- as.matrix(NA, nrow=S, ncol=1)
        for (s in 1:S) {
            theta_s = as.matrix(theta_post[s,])
            ll <- llFun.mfsurv(theta_s,Y,Y0,C,X,Z,data)
            llSum <- llSum + ll$llik
            sum1 <- sum1 + ll$one
            sum2 <- sum2 + ll$two
            sum3 <- sum3 + ll$three
        }
        P = 2 * (L - (1 / S * llSum))
        #Calculate DIC
        DIC <- -2 * (L - P)
        all <- sum1/S
        finite <- sum2/S
        small <- sum3/S
        list(DIC = DIC, Loglik = L)

    }else{
            X <- as.matrix(object$X)
            Y <- as.matrix(object$Y)
            Y0 <- as.matrix(object$Y0)
            C <- as.matrix(object$C)
            data <- as.data.frame(cbind(object$Y, object$Y0, object$C, object$X))
            theta_post = cbind(object$betas, object$lambda)
            theta_hat = apply(theta_post, 2, mean)
            L = llFun.mcmcsurv(theta_hat,Y, Y0,C,X,data)$llik
            #Calculate P
            S = nrow(theta_post) #S = number of iterations
            #Add up the log likelihoods of each iteration
            llSum = 0
            sum1 = 0
            sum2 = 0
            sum3 = 0
            #l <- as.matrix(NA, nrow=S, ncol=1)
            for (s in 1:S) {
                theta_s = as.matrix(theta_post[s,])
                ll <- llFun.mcmcsurv(theta_s,Y,Y0,C,X,data) ### llFun <--> llFun.mcmcsurv
                llSum <- llSum + ll$llik
                sum1 <- sum1 + ll$one
                sum2 <- sum2 + ll$two
                sum3 <- sum3 + ll$three
            }
            P = 2 * (L - (1 / S * llSum))
            #Calculate DIC
            DIC <- -2 * (L - P)
            all <- sum1/S
            finite <- sum2/S
            small <- sum3/S
            list(DIC = DIC, Loglik = L)
    }
}

