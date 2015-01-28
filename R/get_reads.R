#' get sequencing reads from fragments
#'
#' simulate the sequencing process by returning the sequence of one or both ends of provided 
#' fragments
#' 
#' @param tFrags DNAStringSet representing fragments
#' @param readlen Read length.
#' @param paired If \code{FALSE}, return only the first \code{readlen} bases of 
#'   each element of \code{tFrags} in the result; if \code{TRUE}, also return
#'   last \code{readlen} bases.
#' @export
#' @return DNAStringSet representing simulated RNA-seq reads
#' @seealso \code{\link{simulate_experiment}}, 
#'   \code{\link{simulate_experiment_countmat}}
#' @examples
#'   library(Biostrings)
#'   data(srPhiX174)
#'   set.seed(174)
#'   srPhiX174_reads = get_reads(srPhiX174, readlen=15, paired=FALSE)
#'   srPhiX174_reads  
#'   # set of single-end, 15bp reads, treating srPhiX174 as the fragments
get_reads = function(tFrags, readlen, paired=TRUE){
  
    # when fragments are shorter than reads:
    isShort = (width(tFrags) <= readlen)
    isLong = !isShort
      
    if(paired) {
      
        if(sum(isShort) > 0){
            x = tFrags[isShort] # left mate
            names(x) = paste0(seq(along=x), "a")
            rc = reverseComplement(x) # right mate
            names(rc) = paste0(seq(along=x), "b")
            out = c(x,rc)
            outInds = rep(1:length(x), each=2)
            outInds[seq(2, length(outInds), by=2)] = (1:length(x))+length(x)
            outShort = out[outInds] # puts pairs of reads next to each other
            names(outShort) = paste0(rep(names(tFrags)[isShort], each=2))
        }
    
        if(sum(isLong) > 0){
            x = tFrags[isLong]
            lr = subseq(x, start=1, end=readlen) # left mate
            names(lr) = paste0(seq(along=x), "a")
            rr = subseq(x, start=(width(x)-readlen+1), end=width(x))
            rr = reverseComplement(rr) # right mate
            names(rr) = paste0(seq(along=x), "b")
            out = c(lr, rr)
            outInds = rep(1:length(lr), each=2)
            outInds[seq(2, length(outInds), by=2)] = (1:length(lr))+length(lr)
            outLong = out[outInds] # puts pairs of reads next to each other
            names(outLong) = paste0(rep(names(tFrags)[isLong], each=2))    
        }
      
        if(sum(isLong) > 0 & sum(isShort) > 0) {
            theReads = c(outLong, outShort)
        } else if(sum(isLong) > 0) {
            theReads = outLong
        } else {
            theReads = outShort
        }
      
        return(theReads)
    
    } else { #single end
      theReads = tFrags
      theReads[isLong] = subseq(tFrags[isLong], start=1, end=readlen)
      return(theReads)
    }
    
}
