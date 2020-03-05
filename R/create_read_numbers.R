#' Generate a simulated data set based on known model parameters
#'
#' @param mu Baseline mean expression for negative binomial model
#' @param fit Fitted relationship between log mean and log size
#' @param p0 A vector of the probabilities a count is zero 
#' @param m Number of genes/transcripts to simulate (not necessary if mod, 
#'   beta are specified, but can still be specified if preserve.beta.order=TRUE)
#' @param n Number of samples to simulate (not necessary if mod, beta are 
#'   specified)
#' @param mod  Model matrix you would like to simulate from without an intercept
#' @param beta set of coefficients for the model matrix (must have same number 
#'   of columns as mod). Natural log scale, i.e. ln(expression ratio), is assumed.
#' @param seed optional seed to set (for reproducibility)
#' @param preserve.names optional, use names(mu) as rownames for output
#' @param preserve.beta.order optional. Each beta value corresponds to a specific
#'   mu value. length(beta) must be equal to length(mu); m can be given to limit 
#'   number of genes to simulate.
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
#'     p0=params$p0, m=Ntranscripts, n=Nsamples, seed=103, preserve.names=TRUE)
#'   # plot simulated readcounts vs input mu
#'   plot(x=log10(params$mu[rownames(custom_readmat)]), y=log10(rowSums(custom_readmat)))
#'   
#'   mod = matrix(c(rep(0,10),rep(1,10))) # 10 test, 10 control samples
#'   beta = matrix(c(rep(1,20),rep(-1,20), rep(0,length(params$mu)-40))) 
#'   # 20 upregulated, 20 downregulated genes (2.72-fold)
#'   custom_DE = create_read_numbers(mu=params$mu, fit=params$fit, 
#'     p0=params$p0, mod=mod, beta=beta, seed=103, preserve.beta.order=TRUE, 
#'     preserve.names=TRUE, m=70)
#'   # get crude log2 fold-change with a 0.1 prior count
#'   l2fc <- cbind(log2((rowSums(custom_DE[,mod==1])+0.5)/(rowSums(custom_DE[,mod!=1])+0.5)))
#'   # plot fold-change (natural log scale)
#'   plot(l2fc, col = factor(beta[names(params$mu) %in% rownames(custom_DE),]), ylim=c(-3,3))  
#'   # 17 up; 14 down; 39 unregulated


create_read_numbers = function(mu, fit, p0, m=NULL, n=NULL, mod=NULL, beta=NULL, 
    seed=NULL, preserve.names = FALSE, preserve.beta.order = FALSE){
 
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
        if (preserve.beta.order) {
          m = ifelse(is.null(m), dim(beta)[1], m)
        } else {
          m = dim(beta)[1]
        }
        n = dim(mod)[1]
        index = sample(1:length(mu),size=m)
        mus = mu[index]
        p0s = p0[index]
  
        ind = !apply(mod,2,function(x){all(x==1)})
        mod = cbind(mod[,ind])
        beta = cbind(beta[,ind])
        if(preserve.beta.order){
          if(length(beta) == length(mu)){
            beta = beta[index,]
          } else {
            warning('Warning: can\'t preserve beta order with gene mus. length(beta) not equal to length(mu).')
          }
        }
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
    if (preserve.names){
      if (is.null(names(mu))){
        warning('Warning: preserve.names was set to TRUE, but elements of mu are not named. Naming the output rows by original index of mu.')
        nms = as.character(1:length(mu))
      } else {
        nms = names(mu)
      }
      rownames(counts) = nms[index]
    }
    counts <- counts[order(index),]
    return(counts)
}
