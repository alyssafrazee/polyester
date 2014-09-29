#' Simulate RNA-seq experiment 
#'
#' create FASTA files containing RNA-seq reads simulated from provided 
#'   transcripts, with optional differential expression between two groups
#'   (designated via read count matrix)
#' @param fasta path to FASTA file containing transcripts from which to simulate
#'    reads. See details.
#' @param gtf path to GTF file containing transcript structures from which reads
#'   should be simulated. See details.
#' @param seqpath path to folder containing one FASTA file (\code{.fa} 
#'   extension) for each chromosome in \code{gtf}. See details. 
#' @param readmat matrix with rows representing transcripts and columns 
#'   representing samples. Entry i,j specifies how many reads to simulate from
#'   transcript i for sample j.
#' @param outdir character, path to folder where simulated reads should be 
#'   written, without a slash at the end of the folder name. By default, reads
#'   written to the working directory.
#' @param fraglen Mean RNA fragment length. Sequences will be read off the 
#'   end(s) of these fragments.
#' @param fragsd Standard deviation of fragment lengths. 
#' @param readlen Read length
#' @param error_rate Sequencing error rate. Must be between 0 and 1. A uniform 
#'   error model is assumed. 
#' @param error_model one of \code{'uniform'}, \code{'custom'}, 
#'   \code{'illumina4'}, \code{'illumina5'}, or \code{'roche454'} specifying
#'   which sequencing error model to use while generating reads. See 
#'   \code{?add_platform_error} for more information.
#' @param model_path If using a custom error model, the output folder you
#'   provided to \code{build_error_model.py}. Should contain either two files 
#'   suffixed _mate1 and _mate2, or a file suffixed _single.
#' @param model_prefix If using a custom error model, the prefix argument you
#'   provided to \code{build_error_model.py}. This is whatever comes before
#'   _mate1 and _mate2 or _single files in \code{model_path}.
#' @param paired If \code{TRUE}, paired-end reads are simulated; else single-end
#'   reads are simulated.
#' @param seed Optional seed to set before simulating reads, for 
#'   reproducibility.
#' @param ... Further arguments to pass to \code{seq_gtf}, if \code{gtf} is not
#'   \code{NULL}.
#' @return No return, but simulated reads are written to \code{outdir}.
#' @export
#' @details Reads can either be simulated from a FASTA file of transcripts
#'   (provided with the \code{fasta} argument) or from a GTF file plus DNA
#'   sequences (provided with the \code{gtf} and \code{seqpath} arguments). 
#'   Simulating from a GTF file and DNA sequences may be a bit slower: it took
#'   about 6 minutes to parse the GTF/sequence files for chromosomes 1-22, 
#'   X, and Y in hg19.
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
    readmat, outdir=".", fraglen=250, fragsd=25, readlen=100, error_rate=0.005,
    error_model='uniform', model_path=NULL, model_prefix=NULL, paired=TRUE,
    seed=NULL, ...){

    if(!is.null(seed)) set.seed(seed)
    
    if(!is.null(fasta) & is.null(gtf) & is.null(seqpath)){
        transcripts = readDNAStringSet(fasta)
    }else if(is.null(fasta) & !is.null(gtf) & !is.null(seqpath)){
        message('parsing gtf and sequences...')
        transcripts = seq_gtf(gtf, seqpath, ...)
        message('done parsing')
    }else{
        stop('must provide either fasta or both gtf and seqpath')
    }

    stopifnot(class(readmat) == 'matrix')
    stopifnot(nrow(readmat) == length(transcripts))

    # check error model
    error_model = match.arg(error_model, c('uniform', 'illumina4', 'illumina5',
        'roche454', 'custom'))
    if(error_model == 'uniform'){
        stopifnot(error_rate >= 0 & error_rate <= 1)
    }
    if(error_model == 'custom'){
        if(is.null(model_path) | is.null(model_prefix)){
            stop(.makepretty('with custom error models, you must provide both
                the path to the folder that holds your error model
                (model_path) and the prefix of your error model (model_prefix),
                where the prefix is whatever comes before _mate1 and _mate2
                (for paired reads) or _single (for single-end reads). (You
                provided prefix when running build_error_models.py)'))
        }
        if(paired){
            if(!file.exists(paste0(model_path, '/', model_prefix, '_mate1')) |
                !file.exists(paste0(model_path, '/', model_prefix, '_mate2'))){
                stop('could not find error model.')
            }
        }else if(!file.exists(paste0(model_path, '/', model_prefix, '_single'))){
                stop('could not find error model.')
        }
        path = paste0(model_path, '/', model_prefix)
    }
    if(error_model == 'roche454' & paired){
        stop(.makepretty('The Roche 454 error model is only available for
            single-end reads'))
    } 


    sysoutdir = gsub(' ', '\\\\ ', outdir)
    if(.Platform$OS.type == 'windows'){
        shell(paste('mkdir', sysoutdir))
    }else{
        system(paste('mkdir -p', sysoutdir))    
    }

    for(i in 1:ncol(readmat)){
        tObj = rep(transcripts, times=readmat[,i])
  
        #get fragments
        tFrags = generate_fragments(tObj, fraglen=fraglen, fragsd=fragsd)

        #reverse_complement some of those fragments
        rctFrags = reverse_complement(tFrags)

        #get reads from fragments
        reads = get_reads(rctFrags, readlen, paired)

        #add sequencing error
        if(error_model == 'uniform'){
            errReads = add_error(reads, error_rate)            
        }else if(error_model == 'custom'){
            errReads = add_platform_error(reads, 'custom', paired, path)
        }else{
            errReads = add_platform_error(reads, error_model, paired)
        }

        #write read pairs
        write_reads(errReads, readlen=readlen, 
            fname=paste0(outdir, '/sample_', sprintf('%02d', i)), 
            paired=paired)
    }
}