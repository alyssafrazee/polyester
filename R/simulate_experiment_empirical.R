#' Simulate RNA-seq experiment based on abundances from a data set
#'
#' Create fasta files representing reads from a two-group experiment, where
#' abundances and differential expression are estimated from a real data set
#'
#' @param bg Ballgown object containing estimated transcript abundances in FPKM.
#'   Reads will be simulated for the same number of replicates that are in 
#'   \code{bg}. Must provide exactly one of \code{bg} and \code{fpkmMat}.
#' @param fpkmMat transcript-by-sample matrix containing abundances (in FPKM)
#'   estimated from a real data set. MUST have row names identifying
#'   transcripts. The number of columns is the number of samples that will be
#'   simulated.
#' @param mean_rps Number of reads per sample to use in converting FPKM 
#'   measurements to counts. Should be somewhat close to the number of reads
#'   per sample in the experiment that generated the estimated FPKMs. Defaults
#'   to 5 million (5e6). 
#' @param grouplabels vector indicating the group labels for each replicate
#'   in the experiment. Must be convertible to a factor with exactly 2 levels.
#' @param decut A transcript will be recorded as truly differentially expressed
#'   if its fold change between the two groups is more extreme than
#'   \code{decut}, in either direction.
#' @param outdir character, path to folder where simulated reads should be 
#'   written, without a slash at the end of the folder name. By default, reads
#'   written to the working directory.
#' @param ncores the number of cores to be utilized for parallel generation
#'   of simulated reads. Note that if more than one core is specified,
#'   the code will parallelize by replicate, so if num_reps == 1, this
#'   will not be utilized.
#' @param ... Additional arguments to pass to simulate_experiment_countmat
#' @importFrom limma lmFit
#' @export
#' @return No return, but reads are written to \code{outdir}.
#' @examples \dontrun{
#' 
#'   library(ballgown) 
#'   data(bg)
#'   bg = subset(bg, "chr=='22'")
#'   
#'   # load gtf file:
#'   gtfpath = system.file('extdata', 'bg.gtf.gz', package='polyester')
#'   gtf = subset(gffRead(gtfpath), seqname=='22')
#'  
#'   # load/download chromosome sequence (just for this example)
#'   system('wget https://www.dropbox.com/s/04i6msi9vu2snif/chr22seq.rda')
#'   load('chr22seq.rda')
#'   names(chr22seq) = '22'
#'
#'   # simulate reads based on this experiment's FPKMs
#'   simulate_experiment_empirical(bg, grouplabels=pData(bg)$group, gtf=gtf, 
#'      seqpath=chr22seq, mean_rps=5000, outdir='simulated_reads_3', seed=1247)
#'}
simulate_experiment_empirical = function(bg=NULL, fpkmMat=NULL,
    mean_rps=5e6, grouplabels=NULL, decut=1.5, outdir='.', ncores = 1L, ...){

    if(!xor(is.null(bg), is.null(fpkmMat))){
        stop('must provide exactly one of bg or fpkmMat')
    }
    if(is.null(fpkmMat)){
        if(!('all' %in% bg@meas | !'FPKM' %in% bg@meas)){
            stop('bg does not contain transcript FPKMs')
        }
        fpkmMat = ballgown::texpr(bg, 'FPKM')
        rownames(fpkmMat) = ballgown::texpr(bg, 'all')$t_name
    }
    if(is.null(bg)){
        if(is.null(rownames(fpkmMat))){
            stop('fpkmMat must have row names identifying transcripts')
        }
    }

    # get correct ordering + transcript lengths
    extras = list(...)
    if('fasta' %in% names(extras)){
        tt = readDNAStringSet(extras$fasta)
    }else{
        if(!('gtf' %in% names(extras) & 'seqpath' %in% names(extras))){
            stop('must provide either fasta or gtf + seqpath')
        }
        tt = seq_gtf(extras$gtf, extras$seqpath)
    }
    
    # If for some reason the names of the transcripts are missing,
    # stop.
    if (is.null(names(tt))) {
      stop(.makepretty('the transcripts do not appear to have names;
           something went wrong with the import process or the provided
           file has empty names for each transcript; please provide unique
           names that match the row names in the ballgown object/matrix.'))
    }

    fastanames = names(tt)
    matrixnames = rownames(fpkmMat)
    if(length(fastanames) != length(matrixnames)){
        stop(.makepretty('fasta and provided ballgown object/matrix must
            annotate the same number of transcripts'))
    }
    if(any(sort(fastanames) != sort(matrixnames))){
        stop(.makepretty(paste0('ballgown object/FPKM matrix does not
            contain the same transcripts as fasta. fasta example: ', 
            sort(fastanames)[1], '; bg/fpkmMat example: ', 
            sort(matrixnames)[1])))
    }
    o = match(fastanames, matrixnames)
    fpkmMat = fpkmMat[o,]
    stopifnot(all(rownames(fpkmMat) == names(tt))) #should never be reached

    tlengths = width(tt)

    countmat = fpkm_to_counts(mat=fpkmMat, tlengths=tlengths, mean_rps=mean_rps)
    logcountmat = log2(countmat+1)
    stopifnot(length(unique(grouplabels)) == 2)
    x = model.matrix(~grouplabels)
    fit = lmFit(logcountmat, x)
    truebetas = 2^(fit$coefficients[,2])
    refcat = levels(as.factor(unique(grouplabels)))[1]
    coefcat = levels(as.factor(unique(grouplabels)))[2]
    print(paste0('fold change direction: ', coefcat, '/', refcat))

    # prep output directory
    sysoutdir = gsub(' ', '\\\\ ', outdir)
    if(.Platform$OS.type == 'windows'){
        shell(paste('mkdir', sysoutdir))
    }else{
        system(paste('mkdir -p', sysoutdir))    
    }

    # write table of transcript statuses
    isDE = truebetas > decut | truebetas < (1/decut)
    sim_info = data.frame(transcriptid=rownames(countmat), 
        foldchange=as.numeric(truebetas), DEstatus=isDE)
    write.table(sim_info, quote=FALSE, row.names=FALSE, sep='\t',
        file=paste0(outdir, '/sim_tx_info.txt'))

    # write table of sample information
    rep_info = data.frame(
        rep_id=paste0('sample_', sprintf('%02d', 1:ncol(fpkmMat))),
        group=grouplabels)
    write.table(rep_info, row.names=FALSE, quote=FALSE, sep='\t', 
        file=paste0(outdir, '/sim_rep_info.txt'))

    simulate_experiment_countmat(readmat=countmat, outdir=outdir, ncores = ncores, ...)
}

