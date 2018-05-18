#' generate a set of fragments from a set of transcripts
#'
#' Convert each sequence in a DNAStringSet to a "fragment" (subsequence)
#' @param tObj DNAStringSet of sequences from which fragments should be 
#'   extracted
#' @param distr One of 'normal', 'empirical', or 'custom'. If 'normal', draw 
#'   fragment lengths from a normal distribution with mean \code{fraglen} and 
#'   standard deviation \code{fragsd}. If 'empirical', draw fragment lengths 
#'   from a fragment length distribution estimated from a real data set. If
#'   'custom', draw fragment lengths from a custom distribution, provided as
#'   the \code{custdens} argument, which should be a density fitted using 
#'   \code{\link{logspline}}. 
#' @param fraglen Mean fragment length, if drawing fragment lengths from a
#'   normal distribution. 
#' @param fragsd Standard deviation of fragment lengths, if drawing lengths
#'   from a normal distribution. Note: \code{fraglen} and \code{fragsd} are 
#'   ignored unless \code{distr} is 'normal'.
#' @param readlen Read length. Default 100. Used only to label read positions.
#' @param custdens If \code{distr} is 'custom', draw fragments from this
#'   density. Should be an object of class \code{logspline}.
#' @param bias One of 'none', 'rnaf', or 'cdnaf' (default 'none'). 'none' 
#'   represents uniform fragment selection (every possible fragment in a 
#'   transcript has equal probability of being in the experiment); 'rnaf'
#'   represents positional bias that arises in protocols using RNA
#'   fragmentation, and 'cdnaf' represents positional bias arising in protocols
#'   that use cDNA fragmentation (Li and Jiang 2012). Using the 'rnaf' model,
#'   coverage is higher in the middle of the transcript and lower at both ends,
#'   and in the 'cdnaf' model, coverage increases toward the 3' end of the
#'   transcript. The probability models used come from Supplementary Figure S3
#'   of Li and Jiang (2012).
#' @param frag_GC_bias See explanation in \code{\link{simulate_experiment}}.
#' @export
#' @return DNAStringSet consisting of one randomly selected subsequence per 
#'   element of \code{tObj}.
#' @details
#'   The empirical fragment length distribution was estimated using 7 randomly 
#'   selected RNA-seq samples from the GEUVADIS dataset ('t Hoen et al 2013),
#'   one sample from each laboratory that performed sequencing for that data 
#'   set. We used Picard's "CollectInsertSizeMetrics" 
#'   (http://broadinstitute.github.io/picard/), version 1.121, to estimate the
#'   insert size distribution based on the read alignments. 
#' @references
#'   't Hoen PA, et al (2013): Reproducibility of high-throughput mRNA and 
#'   small RNA sequencing across laboratories. Nature Biotechnology 31(11):
#'   1015-1022.
#'  
#'   Li W and Jiang T (2012): Transcriptome assembly and isoform expression
#'   level estimation from biased RNA-Seq reads. Bioinformatics 28(22):
#'   2914-2921.
#'
#' @seealso \code{\link{logspline}}
#'
#' @examples
#'   library(Biostrings)
#'   data(srPhiX174)
#' 
#'   ## get fragments with lengths drawn from normal distrubution
#'   set.seed(174)
#'   srPhiX174_fragments = generate_fragments(srPhiX174, fraglen=15, fragsd=3, 
#'       readlen=4)
#'   srPhiX174_fragments
#'   srPhiX174
#' 
#'   ## get fragments with lengths drawn from an empirical distribution
#'   empirical_frags = generate_fragments(srPhiX174, distr='empirical')
#'   empirical_frags
#'   
#'   ## get fragments with lengths from a normal distribution, but include
#'   ## positional bias from cDNA fragmentation:
#'   biased_frags = generate_fragments(srPhiX174, bias='cdnaf')
#'   biased_frags
#'  
generate_fragments = function(tObj, fraglen=250, fragsd=25, 
  readlen=100, distr='normal', custdens=NULL, bias='none',
  frag_GC_bias='none') {
  
    bias = match.arg(bias, c('none', 'rnaf', 'cdnaf'))
    distr = match.arg(distr, c('normal', 'empirical', 'custom'))
    L = width(tObj)
    if(distr == 'empirical'){
        data('empirical_density')
        fraglens = round(rlogspline(L, empirical_density))
    }else if(distr == 'normal'){
        fraglens = round(rnorm(L, mean=fraglen, sd=fragsd)) 
    }else{
        # distr == 'custom'
        if(is.null(custdens)){
            stop('must provide custom logspline density when distr is "custom"')
        }
        stopifnot(class(custdens) == 'logspline')
        fraglens = round(rlogspline(L, custdens))
    }
    s = which(fraglens < L)
    n = length(s)
    if(bias == 'none'){
        start_pos = floor(runif(n, min=rep(1,n), max=L[s]-fraglens[s]+2))
    }else if(bias == 'rnaf'){
        data(rnaf)
        starts_pct = sample(rnaf$pospct, size=n, prob=rnaf$prob, replace=TRUE)
        starts_pct[starts_pct==1] = 0.999
        start_pos = floor(starts_pct * (L[s]-fraglens[s]+2))
        start_pos[start_pos==0] = 1
    }else{
        # bias == 'cdnaf'
        data(cdnaf)
        starts_pct = sample(cdnaf$pospct, size=n, prob=cdnaf$prob, replace=TRUE)
        starts_pct[starts_pct==1] = 0.999
        start_pos = floor(starts_pct * (L[s]-fraglens[s]+2))
        start_pos[start_pos==0] = 1
    }
    tObj[s] = subseq(tObj[s], start=start_pos, width=fraglens[s])
    names(tObj)[s] = paste0(names(tObj[s]), ';mate1:', start_pos, '-', 
        start_pos+readlen-1, ';mate2:', start_pos+fraglens[s]-readlen, '-', 
        start_pos+fraglens[s]-1)
    nonseqinds = (1:length(tObj))[-s]
    names(tObj)[nonseqinds] = paste0(names(tObj[nonseqinds]), 
        ';mate1Start:1;mate2Start:1')

    # fragment GC bias coin flips (Bernoulli trials)
    gc <- as.numeric(letterFrequency(tObj, "GC", as.prob=TRUE))
    if (is.numeric(frag_GC_bias)) {
      gc.idx <- as.integer(cut(gc, breaks=c(-Inf,(0:99)/100+.005,Inf)))
      prob <- frag_GC_bias[gc.idx]
      stopifnot(all(prob >= 0 & prob <= 1))
      coinflip <- rbinom(length(tObj), 1, prob) # flip a coin
      tObj <- tObj[ coinflip == 1 ] # only return successes
    }

    return(tObj)
}


