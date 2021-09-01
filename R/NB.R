#' Draw nonzero negative binomial random numbers
#'
#' @param  basemeans vector of means, one per draw
#' @param  size vector of size parameters (controlling the mean/variance
#'   relationship); one per draw
#' @param  seed optional seed to set before drawing
#' @return vector of negative binomial draws from specified distributions,
#'   where any zero draw is replaced with a 1. Length of return vector is
#'   equal to \code{length(basemeans)}.
#' @export
#' @examples
#'   randomNBs = NB(c(100, 4, 29), size=c(50, 2, 4), seed=21)
#'   randomNBs  # 115, 5, 15
#'
NB = function(basemeans, size, seed=NULL){
    if(!is.null(seed)) set.seed(seed)
    numreads = rnbinom(n = length(basemeans), mu = basemeans, size = size)
    numreads[numreads == 0] = 1
    return(numreads)
}
