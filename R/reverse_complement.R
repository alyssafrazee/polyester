#' reverse-complement some fragments
#'
#' randomly reverse-complement half of the sequences in a DNAStringSet
#' @param tObj DNAStringSet representing sequences.
#' @param seed optional seed to set before randomly selecting the sequences to
#' be reverse-complemented.
#' @export
#' @return DNAStringSet that is the same as \code{tObj}, but with about half
#' the sequences reverse-complemented.
#' @examples
#' library(Biostrings)
#' data(srPhiX174)
#' srPhiX174_halfrc = reverse_complement(srPhiX174, seed=174)
reverse_complement = function(tObj, library_type, strand_error_rate, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  if(library_type != "unstranded"){
    strand = sample(c(0,1), length(tObj), replace=TRUE, prob = c(1-strand_error_rate, strand_error_rate))
  } else {
    # default unstranded option
    strand = sample(c(0,1), length(tObj), replace=TRUE)
  }
  tObj[strand==0] = reverseComplement(tObj[strand==0])
  return(tObj)
}
