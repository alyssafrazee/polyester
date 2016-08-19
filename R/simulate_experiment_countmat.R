#' Simulate RNA-seq experiment 
#'
#' create FASTA files containing RNA-seq reads simulated from provided 
#'   transcripts, with optional differential expression between two groups
#'   (designated via read count matrix)
#' @param fasta path to FASTA file containing transcripts from which to simulate
#'    reads. See details.
#' @param gtf path to GTF file or data frame containing transcript structures
#'   from which reads should be simulated. See details and 
#'   \code{\link{seq_gtf}}. 
#' @param seqpath path to folder containing one FASTA file (\code{.fa} 
#'   extension) or DNAStringSet containing one entry for each chromosome in 
#'   \code{gtf}. See details and \code{\link{seq_gtf}}.
#' @param readmat matrix with rows representing transcripts and columns 
#'   representing samples. Entry i,j specifies how many reads to simulate from
#'   transcript i for sample j.
#' @param outdir character, path to folder where simulated reads should be 
#'   written, without a slash at the end of the folder name. By default, reads
#'   written to the working directory.
#' @param paired If \code{TRUE}, paired-end reads are simulated; else single-end
#'   reads are simulated.
#' @param seed Optional seed to set before simulating reads, for 
#'   reproducibility.
#' @param ... Additional arguments to add nuance to the simulation, as described
#'   extensively in the details of \code{\link{simulate_experiment}}, or to pass
#'   to \code{seq_gtf}, if \code{gtf} is not \code{NULL}.
#' @return No return, but simulated reads are written to \code{outdir}.
#' @export
#' @details Reads can either be simulated from a FASTA file of transcripts
#'   (provided with the \code{fasta} argument) or from a GTF file plus DNA
#'   sequences (provided with the \code{gtf} and \code{seqpath} arguments). 
#'   Simulating from a GTF file and DNA sequences may be a bit slower: it took
#'   about 6 minutes to parse the GTF/sequence files for chromosomes 1-22, 
#'   X, and Y in hg19.
#' @references
#'   Li W and Jiang T (2012): Transcriptome assembly and isoform expression
#'   level estimation from biased RNA-Seq reads. Bioinformatics 28(22):
#'   2914-2921.
#' @examples \donttest{
#'   fastapath = system.file("extdata", "chr22.fa", package="polyester")
#'   numtx = count_transcripts(fastapath)
#'   readmat = matrix(20, ncol=10, nrow=numtx)
#'   readmat[1:30, 1:5] = 40
#' 
#'   simulate_experiment_countmat(fasta=fastapath, 
#'     readmat=readmat, outdir='simulated_reads_2', seed=5)
#'}

simulate_experiment_countmat = function(fasta=NULL, gtf=NULL, seqpath=NULL, 
    readmat, outdir='.', paired=TRUE, seed=NULL, ...){

    extras = list(...)

    if(!is.null(seed)) set.seed(seed)
    
    if(!is.null(fasta) & is.null(gtf) & is.null(seqpath)){
        transcripts = readDNAStringSet(fasta)
    }else if(is.null(fasta) & !is.null(gtf) & !is.null(seqpath)){
        transcripts = seq_gtf(gtf, seqpath, ...)
    }else{
        stop('must provide either fasta or both gtf and seqpath')
    }

    stopifnot(class(readmat) == 'matrix')
    stopifnot(nrow(readmat) == length(transcripts))

    # validate extra arguments
    extras = .check_extras(extras, paired, total.n=ncol(readmat))

    sysoutdir = gsub(' ', '\\\\ ', outdir)
    if(.Platform$OS.type == 'windows'){
        shell(paste('mkdir', sysoutdir))
    }else{
        system(paste('mkdir -p', sysoutdir))    
    }

    # do the sequencing
    sgseq(readmat, transcripts, paired, outdir, extras)

}
