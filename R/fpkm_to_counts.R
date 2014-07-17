#' Turn FPKMs from a ballgown object into estimated counts for transcripts
#'
#' @param bg ballgown object created from real RNA-seq dataset
#' @param threshold only estimate parameters from transcripts with mean FPKM measurements larger 
#' than \code{threshold}
#' @param mean_rps This should be the number of reads per sample in total for use in backing out the
#' FPKM calculations
#' 
#' @return A matrix of counts with the same number of rows and columns as the ballgown object
#' 
#' @export
#' @author Jeff Leek
#' @examples \dontrun{
#'   require(ballgown)
#'   data(bg)
#'   countmat = fpkm_to_counts(bg, mean_rps=400000) 
#' }
#' 

fpkm_to_counts = function(bg, mean_rps=100e6, threshold=0){
    tmeas = as.matrix(ballgown::texpr(bg, "FPKM"))
    trowm = rowMeans(tmeas)
    index1 = which(trowm > threshold)
    tlengths = sapply(width(ballgown::structure(bg)$trans[index1]), sum)
    counts = tlengths*tmeas[index1,]/1000
    counts = round(counts*mean_rps/1e6)
    return(counts)
}
