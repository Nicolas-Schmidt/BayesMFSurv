
#' @export
mfsurv.stats <- function(object){
    .Deprecated("stats")
    return(stats(object = object))
}

#' @export
mfsurv.summary <- function(object, parameter){
    .Deprecated("summary")
    return(summary(object = object, parameter = parameter))
}

