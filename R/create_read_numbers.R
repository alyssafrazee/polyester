#' Generate a simulated data set based on known model parameters
#'
#' @param mu Baseline mean expression for negative binomial model
#' @param fit Fitted relationship between log mean and log size
#' @param p0 A vector of the probabilities a count is zero
#' @param m Number of genes/transcripts to simulate (not necessary if mod,
#'   beta are specified)
#' @param n Number of samples to simulate (not necessary if mod, beta are
#'   specified)
#' @param mod  Model matrix you would like to simulate from without an intercept
#' @param beta set of coefficients for the model matrix (must have same number
#'   of columns as mod)
#' @param seed optional seed to set (for reproducibility)
#'
#' @return counts Data matrix with counts for genes in rows and samples in
#'   columns
#'
#' @export
#'
#' @author Jeff Leek
#' @examples
#'   library(ballgown)
#'   data(bg)
#'   countmat = fpkm_to_counts(bg, mean_rps=400000)
#'   params = get_params(countmat)
#'   Ntranscripts = 50
#'   Nsamples = 10
#'   custom_readmat = create_read_numbers(mu=params$mu, fit=params$fit,
#'     p0=params$p0, m=Ntranscripts, n=Nsamples, seed=103)
#'


create_read_numbers = function(mu, fit, p0, m=NULL, n=NULL, mod=NULL, beta=NULL,
    seed=NULL){

    if(!is.null(seed)){set.seed(seed)}
    if(is.null(mod) | is.null(beta)){
        cat("Generating data from baseline model.\n")
        if(is.null(m) | is.null(n)){
            stop(.makepretty("create_read_numbers error: if you don't specify
            mod and beta, you must specify m and n.\n"))
        }
        index = sample(1:length(mu),size=m)
        mus = mu[index]
        p0s = p0[index]
        mumat = log(mus + 0.001) %*% t(rep(1,n))
    } else {
        m = dim(beta)[1]
        n = dim(mod)[1]
        index = sample(1:length(mu),size=m)
        mus = mu[index]
        p0s = p0[index]

        ind = !apply(mod,2,function(x){all(x==1)})
        mod = cbind(mod[,ind])
        beta = cbind(beta[,ind])
        mumat = log(mus + 0.001) + beta %*% t(mod)
    }

    muvec = as.vector(mumat)
    sizevec = predict(fit,muvec)$y
    sizemat = matrix(sizevec,nrow=m)
    counts = sizemat*NA
    for(i in 1:m){
      counts[i,] = rbinom(n,prob=(1-p0s[i]),size=1)*
        rnbinom(n,mu=exp(mumat[i,]),size=exp(sizemat[i,]))
    }
    return(counts)
}
