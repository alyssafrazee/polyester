.makepretty = function(x){
    msg = gsub('\n', ' ', x)
    msg = gsub('    ', '', msg)
    msg
}

#' simulate RNA-seq experiment using negative binomial model
#'
#' create FASTA files containing RNA-seq reads simulated from provided 
#'   transcripts, with optional differential expression between two groups
#' @param fasta path to FASTA file containing transcripts from which to simulate
#'   reads. See details.
#' @param gtf path to GTF file containing transcript structures from which reads
#'   should be simulated. See details.
#' @param seqpath path to folder containing one FASTA file (\code{.fa} 
#'   extension) for each chromosome in \code{gtf}. See details. 
#' @param num_reps How many biological replicates should be in each group? If
#'   \code{num_reps} is a single integer, \code{num_reps} replicates will be 
#'   simulated in each group. Otherwise, \code{num_reps} can be a length-2 
#'   vector, where \code{num_reps[1]} and \code{num_reps[2]} replicates will be
#'   simulated in each of the two groups.
#' @param fraglen Mean RNA fragment length. Sequences will be read off the 
#'   end(s) of these fragments.
#' @param fragsd Standard deviation of fragment lengths. 
#' @param readlen Read length.
#' @param lib_sizes Library size factors for the biological replicates. If 
#'   \code{num_reps} is an integer, \code{lib_sizes} should have length 
#'   \code{2*num_reps}. If \code{num_reps} is a length-2 vector, 
#'   \code{lib_sizes} should have length \code{sum(num_reps)}. For each
#'   replicate, once the number of reads to simulate from each transcript for 
#'   that replicate is known, all read numbers across all transcripts from that
#'   replicate are multiplied by the corresponding entry in \code{lib_sizes}. 
#' @param error_rate Sequencing error rate. Must be between 0 and 1. Only used
#'   if error_model is \code{'uniform'}. 
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
#' @param paired If \code{TRUE}, paired-end reads are simulated; else 
#'   single-end reads are simulated.
#' @param reads_per_transcript baseline mean number of reads to simulate 
#'   from each transcript. Can be an integer, in which case this many reads
#'   are simulated from each transcript, or an integer vector whose length
#'   matches the number of transcripts in \code{fasta}. 
#' @param fold_changes Vector of multiplicative fold changes between groups,
#'   one entry per transcript in \code{fasta}. A fold change > 1 means the 
#'   transcript is overexpressed in the first \code{num_reps} (or 
#'   \code{num_reps[1]}) samples. Fold change < 1 means transcript is 
#'   overexpressed in the last \code{num_reps} (or \code{num_reps[2]}) samples.
#'   The change is in the mean number of reads generated from the transcript, 
#'   between groups.
#' @param size the negative binomial \code{size} parameter (see 
#'   \code{\link{NegBinomial}}) for the number of reads drawn per transcript. 
#'   If left blank, defaults to \code{reads_per_transcript / 3}. Negative 
#'   binomial variance is mean + mean^2 / size. Can either be left at default,
#'   a vector of the same length as number of transcripts in \code{fasta},
#'   if the two groups should have the same size parameters, or a list with 2
#'   elements, each of which is a vector with length equal to the number of 
#'   transcripts in \code{fasta}, which represent the size parameters for each
#'   transcript in groups 1 and 2, respectively.
#' @param outdir character, path to folder where simulated reads should be 
#'   written, with *no* slash at the end. By default, reads are 
#'   written to current working directory.
#' @param write_info If \code{TRUE}, write a file matching transcript IDs to 
#'   differential expression status into the file \code{outdir/sim_tx_info.txt}
#'   and a file matching biological replicate IDs to group membership and 
#'   library size into the file \code{outdir/sim_rep_info.txt}. 
#' @param transcriptid optional vector of transcript IDs to be written into 
#'   \code{sim_info.txt} and used as transcript identifiers in the fasta files.
#'   Defaults to \code{names(readDNAStringSet(fasta))}. This option is useful
#'   if default names are very long or contain special characters.
#' @param seed Optional seed to set before simulating reads, for 
#'   reproducibility.
#' @param ... additional arguments to pass to \code{seq_gtf} if using 
#'   \code{gtf} and \code{seqpath}
#' @return No return, but simulated reads and a simulation info file are written
#'   to \code{outdir}.
#' @export
#' @details Reads can either be simulated from a FASTA file of transcripts 
#'   (provided with the \code{fasta} argument) or from a GTF file plus DNA 
#'   sequences (provided with the \code{gtf} and \code{seqpath} arguments). 
#'   Simulating from a GTF file and DNA sequences may be a bit slower: it took 
#'   about 6 minutes to parse the GTF/sequence files for chromosomes 1-22, X, 
#'   and Y in hg19.
#' @examples \donttest{
#'   ## simulate a few reads from chromosome 22
#' 
#'   fastapath = system.file("extdata", "chr22.fa", package="polyester")
#'   numtx = count_transcripts(fastapath)
#'   set.seed(4)
#'   fold_changes = sample(c(0.5, 1, 2), size=numtx, 
#'      prob=c(0.05, 0.9, 0.05), replace=TRUE)
#'   library(Biostrings)
#'   # remove quotes from transcript IDs:
#'   tNames = gsub("'", "", names(readDNAStringSet(fastapath))) 
#' 
#'   simulate_experiment(fastapath, reads_per_transcript=10, 
#'      fold_changes=fold_changes, outdir='simulated_reads', 
#'      transcriptid=tNames, seed=12)
#'}
simulate_experiment = function(fasta=NULL, gtf=NULL, seqpath=NULL, num_reps=10, 
    fraglen=250, fragsd=25, readlen=100, lib_sizes=NULL, error_rate=0.005,
    error_model='uniform', model_path=NULL, model_prefix=NULL, paired=TRUE, 
    reads_per_transcript=300, fold_changes, size=NULL, outdir=".", 
    write_info=TRUE, transcriptid=NULL, seed=NULL, ...){

    if(!is.null(fasta) & is.null(gtf) & is.null(seqpath)){
        transcripts = readDNAStringSet(fasta)
    }else if(is.null(fasta) & !is.null(gtf) & !is.null(seqpath)){
        message('parsing gtf and sequences...')
        transcripts = seq_gtf(gtf, seqpath, ...)
        message('done parsing')
    }else{
        stop('must provide either fasta or both gtf and seqpath')
    }

    # check argument lengths:
    ntx = length(transcripts)
    if(length(fold_changes)!=ntx){
        stop(.makepretty('fold_changes must be a vector with one entry per
            transcript in fasta or gtf; use count_transcripts to find out how
            many transcripts are in fasta or gtf.'))
    }
    if(length(reads_per_transcript)!=1 & length(reads_per_transcript)!=ntx){
        stop(.makepretty('reads_per_transcript must be a single number or a
            vector with one entry per transcript in fasta or gtf; use
            count_transcripts to find out how many transcripts are in fasta or
            gtf.'))
    }
    if(length(size)>2 & length(size)!=ntx){
        stop(.makepretty('size must be a single number or a vector with one
            entry per transcript in fasta or gtf; use count_transcripts to find
            out how many transcripts are in fasta or gtf.'))
    }
    if(length(transcriptid)!=0 & length(transcriptid)!=ntx){
        stop(.makepretty('transcriptid must be a character vector with one entry
            per transcript in fasta or gtf; use count_transcripts to find out
            how many transcripts are in fasta or gtf.'))
    }
    nreps = ifelse(length(num_reps)==1, 2*num_reps, sum(num_reps))
    if(is.null(lib_sizes)){
        lib_sizes = rep(1, nreps)
    }else{
        if(length(lib_sizes) != nreps){
            stop(.makepretty('lib_sizes must have length equal to total number
                of reps in the experiment, which is either 2*num_reps if
                num_reps is a single number, or sum(num_reps) if num_reps is
                a length-2 vector.'))
        }
        stopifnot(all(lib_sizes > 0))
    }

    # check error model:
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

    L = width(transcripts)
    if(!is.null(transcriptid)){
        names(transcripts) = transcriptid
    }
        
    if(!is.null(seed)) set.seed(seed)

    ## get number of reps per group
    stopifnot(length(num_reps)==1 | length(num_reps)==2)
    if(length(num_reps)==2){
        n1 = num_reps[1]
        n2 = num_reps[2]
    }else{
        n1 = n2 = num_reps
    }

    ## get baseline means for each group from fold changes:
    if(exp(mean(log(fold_changes))) != 1){
        warning(.makepretty('provided fold changes mean that one group will
            have more overall expression than the other, so some of the DE
            signal might be lost due to library size adjustment.'))
    }
    basemeans1 = ifelse(fold_changes > 1, 
        reads_per_transcript*fold_changes, reads_per_transcript)
    basemeans1 = round(basemeans1)
    basemeans2 = ifelse(fold_changes < 1,
        reads_per_transcript/fold_changes, reads_per_transcript)
    basemeans2 = round(basemeans2)

    if (is.null(size)) {
        nbdp1 = basemeans1/3
        nbdp2 = basemeans2/3
    }
    else if(class(size)=='list'){
        stopifnot(length(size)==2)
        stopifnot(length(size[[1]]) == length(basemeans1))
        stopifnot(length(size[[2]]) == length(basemeans2))
        nbdp1 = size[[1]]
        nbdp2 = size[[2]]
    }else{
        if(length(size) == 1 | length(size) == length(basemeans1)){
            nbdp1 = nbdp2 = size
        }else{
            stop(.makepretty('size must either be NULL, length 1,
                or a vector with the same number of elements as
                transcripts you have.'))
        }
    }

    numreadsList = vector("list", n1+n2)
    for(i in 1:n1){
        numreadsList[[i]] = NB(basemeans1, nbdp1)
    }
    for(i in (1:n2)+n1){
        numreadsList[[i]] = NB(basemeans2, nbdp2)
    }

    ## simulate reads for each sample:
    #######################################
    sysoutdir = gsub(' ', '\\\\ ', outdir)
    if(.Platform$OS.type == 'windows'){
        shell(paste('mkdir', sysoutdir))
    }else{
        system(paste('mkdir -p', sysoutdir))    
    }
    for(i in 1:(n1+n2)){
        
        tObj = rep(transcripts, times=ceiling(numreadsList[[i]]*lib_sizes[i]))
        
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
            fname=paste0(outdir, '/sample_', sprintf('%02d', i)), paired=paired)
    }

    ## write out simulation information, if asked for:
    if(write_info){
        if(is.null(transcriptid)){
            transcriptid = names(transcripts)
        }
        sim_info = data.frame(transcriptid=transcriptid, 
            foldchange=fold_changes, DEstatus=fold_changes!=1)
        write.table(sim_info, row.names=FALSE, quote=FALSE, sep="\t", 
            file=paste0(outdir, '/sim_tx_info.txt'))

        rep_info = data.frame(
            rep_id=paste0('sample_', sprintf('%02d', 1:(n1+n2))),
            group=rep(c(1,2), times=c(n1, n2)))
        if(is.null(lib_sizes)){
            rep_info$lib_size = 1
        }else{
            rep_info$lib_size = lib_sizes
        }
        write.table(rep_info, row.names=FALSE, quote=FALSE, sep='\t', 
            file=paste0(outdir, '/sim_rep_info.txt'))
    }
}
