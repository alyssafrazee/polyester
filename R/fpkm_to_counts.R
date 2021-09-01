#' Turn FPKMs from a ballgown object into estimated counts for transcripts
#'
#' @param bg ballgown object created from real RNA-seq dataset
#' @param mat matrix of isoform-level FPKMs from which to derive counts. Rows
#'   should represent transcripts and columns should represent counts. Provide
#'   exactly one of \code{bg} or \code{mat}.
#' @param tlengths if using \code{mat} instead of \code{bg}, vector of
#'   transcript lengths. Entries correspond to the rows of \code{mat}. Lengths
#'   should only count the nucleotides within transcripts' exons.
#' @param mean_rps This should be the number of reads per sample in total for
#'   use in backing out the FPKM calculations.
#' @param threshold only estimate parameters from transcripts with mean FPKM
#'   measurements at least as large as \code{threshold}.
#'
#' @return A matrix of counts with the same number of rows and columns as the
#'   ballgown object
#'
#' @details If transcripts/exons are represented by \code{GRanges} or
#'   \code{GRangesList} objects, the \code{width} function is really useful
#'   in calculating transcript lengths.
#' @export
#' @author Jeff Leek
#' @examples
#'   library(ballgown)
#'   data(bg)
#'   countmat = fpkm_to_counts(bg, mean_rps=400000)
#'

fpkm_to_counts = function(bg=NULL, mat=NULL, tlengths=NULL, mean_rps=100e6,
    threshold=0){
    if(is.null(mat)){
        tmeas = as.matrix(ballgown::texpr(bg, 'FPKM'))
        tlengths = sapply(width(ballgown::structure(bg)$trans), sum)
    }else{
        tmeas = mat
        stopifnot(!is.null(tlengths))
    }
    index1 = which(rowMeans(tmeas) >= threshold)
    tlengths = tlengths[index1]
    counts = tlengths*tmeas[index1,]/1000
    counts = round(counts*mean_rps/1e6)
    return(counts)
}
