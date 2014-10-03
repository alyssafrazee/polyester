#' generate a set of fragments from a set of transcripts
#'
#' Convert each sequence in a DNAStringSet to a "fragment" (subsequence)
#' @param tObj DNAStringSet of sequences from which fragments should be 
#'   extracted
#' @param fraglen Mean fragment length. 
#' @param fragsd Standard deviation of fragment length. Fragment lengths are 
#'   drawn from a normal 
#' distribution with mean \code{fraglen} and standard deviation \code{fragsd}. 
#' @export
#' @return DNAStringSet consisting of one randomly selected subsequence per 
#'   element of \code{tObj}.
#' @examples
#'   library(Biostrings)
#'   data(srPhiX174)
#' 
#'   ## get fragments with lengths drawn from normal distrubution
#'   set.seed(174)
#'   srPhiX174_fragments = generate_fragments(srPhiX174, fraglen=15, fragsd=3)
#'   srPhiX174_fragments
#'   srPhiX174
#' 
#'   ## get fragments with lengths drawn from an empirical distribution
#'   empirical_frags = generate_fragments(srPhiX174, empirical=TRUE)
#'   empirical_frags
#'   
generate_fragments = function(tObj, fraglen=250, fragsd=25, empirical=FALSE){
    L = width(tObj)
    if(empirical){
        data('empirical_density')
        fraglens = round(rlogspline(L, empirical_density))
    }else{
        fraglens = round(rnorm(L, mean=fraglen, sd=fragsd)) 
    }    
    s = which(fraglens < L)
    tObj[s] = subseq(tObj[s], 
        start = floor(runif(length(s), 
            min=rep(1,length(s)), max=L[s]-fraglens[s])), 
        width=fraglens[s])
    return(tObj)
}


