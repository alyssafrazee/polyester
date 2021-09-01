#' add GC bias to a count matrix
#'
#' Given a matrix with rows corresponding to transcripts and sample-specific
#' GC bias models, bias the count matrix using the bias model.
#' 
#' @param readmat matrix of counts, with rows corresponding to features
#'   (transcripts) and columns corresponding to replicates
#' @param gcbias List of GC bias models to add to readmat. Must have length
#'   equal to the number of columns of \code{readmat}. List elements must
#'   either be integers 0 through 7, where 0 means no bias and 1-7 correspond
#'   to built-in GC bias models, or objects of class \code{loess} which can
#'   predict a deviation from overall mean count (on the log scale) given
#'   a GC percentage between 0 and 1. 
#' @param transcripts \code{DNAStringSet} object containing the sequences of
#'   the features (transcripts) corresponding to the rows of \code{readmat}.
#'   Length must be equal to the number of rows in \code{readmat}.
#' @details Designed for internal use in \code{simulate_experiment} functions.
#' @return matrix of the same size as \code{readmat}, but with counts for each
#'   replicate biased according to \code{gcbias}. 
#' @examples
#'   library(Biostrings)
#'   fastapath = system.file("extdata", "chr22.fa", package="polyester")
#'   numtx = count_transcripts(fastapath)
#'   transcripts = readDNAStringSet(fastapath)
#' 
#'   # create a count matrix:
#'   readmat = matrix(20, ncol=10, nrow=numtx)
#'   readmat[1:30, 1:5] = 40
#'
#'   # add biases randomly: use built-in bias models
#'   set.seed(137)
#'   biases = sample(0:7, 10, replace=TRUE)
#'   readmat_biased = add_gc_bias(readmat, as.list(biases), transcripts)
#' @export
#' 

add_gc_bias = function(readmat, gcbias, transcripts){

    stopifnot(length(transcripts) == nrow(readmat))
    stopifnot(length(gcbias) == ncol(readmat))
    
    GC = letterFrequency(transcripts, letters='GC', as.prob=TRUE)
    
    for(i in 1:ncol(readmat)){
        if(class(gcbias[[i]]) == 'loess'){
            fit = gcbias[[i]]
        }else if(gcbias[[i]] != 0){
            eval(parse(text=paste0('data(loessfit', gcbias[[i]], ')')))
            eval(parse(text=paste0('x <- loessfit', gcbias[[i]], '$x')))
            eval(parse(text=paste0('y <- loessfit', gcbias[[i]], '$y')))
            fit = loess(y~x, span=0.3)
        }
        shifts = predict(fit, GC)
        adjcounts = log2(readmat+1)[,i] + shifts
        readmat[,i] = round(2^adjcounts - 1)
    }
    return(readmat)
}