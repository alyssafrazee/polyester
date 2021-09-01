#' Estimate zero-inflated negative binomial parameters from a real dataset
#' 
#' This function estimates the parameters of a zero inflated negative binomial
#' distribution based on a real count data set based on the method of moments.
#' The function also returns a spline fit of log mean to log size which can be
#' used when generating new simulated data. 
#'
#' @param counts A matrix of counts. If you want to simulate from a ballgown
#'   object, see \code{\link{fpkm_to_counts}}
#' @param threshold Only estimate parameters from transcripts with row means 
#'   greater than threshold
#' 
#' @return p0 A vector of probabilities that the count will be zero, one for 
#'   each gene/transcript.
#' @return mu The estimated negative binomial mean by method of moments for the
#'   non-zero counts
#' @return size The estimated negative binomial size by method of moments for
#'   the non-zero counts
#' @return fit A fit relating log mean to log size for use in simulating new
#'   data. 
#' 
#' @export
#' @author Jeff Leek
#' @examples
#'   library(ballgown)
#'   data(bg)
#'   countmat = fpkm_to_counts(bg, mean_rps=400000)
#'   params = get_params(countmat)
#'

get_params = function(counts, threshold=NULL){
  
    if(!is.null(threshold)){
        rowm = rowMeans(counts)
        index1 = which(rowm > threshold)
        counts = counts[index1,]
    }
  
    nsamples = dim(counts)[2]
    counts0 = counts==0
    nn0 = rowSums(!counts0)
    if(any(nn0 == 1)){
        # need more than 1 nonzero count to estimate variance
        counts = counts[nn0 > 1, ]
        nn0 = nn0[nn0 > 1]
        counts0 = counts==0
    }
    mu = rowSums((!counts0)*counts)/nn0
    s2 = rowSums((!counts0)*(counts - mu)^2)/(nn0-1)
    size = mu^2/(s2-mu + 0.0001)
    size = ifelse(size > 0, size, min(size[size > 0]))
    p0 = (nsamples-nn0)/nsamples

    lsize = log(size)
    lmu = log(mu + 0.0001)
    fit = smooth.spline(lsize ~ lmu)
    return(list(p0=p0, mu=mu, size=size, fit=fit))
}

