#' add sequencing error to simulated reads
#'
#' simulate sequencing error by randomly changing the sequenced nucleotide on
#'   some of the reads 
#' @param tFrags DNAStringSet representing sequencing reads
#' @param  error_rate error probability
#' @export
#' @return DNAStringSet equivalent to \code{tFrags} but with random sequencing
#'   errors inserted
#' @examples
#'   library(Biostrings)
#'   data(srPhiX174)
#'   set.seed(174)
#'   srPhiX174_withError = add_error(srPhiX174)
#'   #error was introduced in, e.g., position 10 of 2nd string in set.
add_error = function(tFrags, error_rate = 0.005){
    adj_error = error_rate*5/4 
    #^so you don't have to choose *another* nucleotide for an error: just 
    # choose *a* nucleotide.
    
    allSeq = unlist(tFrags)
    insertLocs = Rle(sample(c(TRUE,FALSE), size = length(allSeq), 
           replace=TRUE, prob = c(adj_error, 1-adj_error)))
  
    newletters = DNAString(
        paste(sample(c('A', 'C', 'G', 'T', 'N'), sum(insertLocs), 
            replace=TRUE), collapse="") )
    allSeq = replaceLetterAt(allSeq, insertLocs, newletters)
    
    eFrags = DNAStringSet(allSeq, 
        start=c(1, (cumsum(width(tFrags))+1)[-length(tFrags)]), 
        width=width(tFrags))

    names(eFrags) = names(tFrags)
    return(eFrags)
}

