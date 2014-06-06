#' Generate a simulated data set based on known model parameters
#'
#' @param mu Baseline mean expression for negative binomial model
#' @param fit Fitted relationship between log mean and log size
#' @param mod  Model matrix you would like to simulate from without an intercept
#' @param beta set of coefficients for the model matrix (must have same number of columns as mod)
#' @param p0 A vector of the probabilities a count is zero 
#' @param seed optional seed to set (for reproducibility)
#' 
#' 
#' @return counts Data matrix with counts for genes in rows and samples in columns
#' 
#' @export
#' 
#' @author Jeff Leek

create_read_numbers = function(mu,fit,mod,beta,p0,seed = NULL){
 
  if(!is.null(seed)){set.seed(seed)}
  
  m = dim(beta)[1]
  n = dim(mod)[1]
  index = sample(1:length(mu),size=m)
  mus = mu[index]
  p0s = p0[index]
  
  ind = !apply(mod,2,function(x){all(x==1)})
  mod = cbind(mod[,ind])
  beta = cbind(beta[,ind])
  mumat = log(mus + 0.001) + beta %*% t(mod)
  
  muvec = as.vector(mumat)
  sizevec = predict(fit,muvec)$y
  p0vec = rep(p0s,m)
  
  counts = matrix(rbinom(m*n,prob=(1-p0vec),size=1)*rnbinom(m*n,size=exp(sizevec),mu=exp(muvec)),nrow=m)
  
  return(counts)
}
