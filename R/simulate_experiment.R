.makepretty = function(x){
    msg = gsub('\n', ' ', x)
    msg = gsub('    ', '', msg)
    msg
}

.check_error_model = function(extras, paired){

    # make sure it's an available model
    error_model = match.arg(extras$error_model,
        c('uniform', 'illumina4', 'illumina5', 'custom'))

    # check uniform model --> error rate
    if(error_model == 'uniform'){
        if('error_rate' %in% names(extras)){
            error_rate = extras$error_rate
            stopifnot(is.numeric(error_rate))
            stopifnot(error_rate >= 0 & error_rate <= 1)
        }
    }

    # check paths and such for custom model
    if(error_model == 'custom'){
        if(!('model_path' %in% names(extras)) |
            !('model_prefix' %in% names(extras))){
            stop(.makepretty('with custom error models, you must provide both
                the path to the folder that holds your error model
                (model_path) and the prefix of your error model (model_prefix),
                where the prefix is whatever comes before _mate1 and _mate2
                (for paired reads) or _single (for single-end reads). You
                provided prefix when running build_error_models.py.'))
        }
        model_path = extras$model_path
        model_prefix = extras$model_prefix
        if(paired){
            if(!file.exists(paste0(model_path, '/', model_prefix, '_mate1')) |
               !file.exists(paste0(model_path, '/', model_prefix, '_mate2'))){
               stop('could not find error model')
            }
        }else{
            if(!file.exists(paste0(model_path, '/', model_prefix, '_single'))){
                stop('could not find error model')
            }
        }
    }
}

.check_fold_changes = function(fold_changes, num_reps, transcripts){

    # make sure fold change matrix is compatible with experiment size
    if(length(num_reps) == 1 | length(num_reps) == 2){
        stopifnot(is.numeric(fold_changes))
    }else{
        stopifnot(is.matrix(fold_changes))
        if(ncol(fold_changes) != length(num_reps)){
            stop(.makepretty('wrong number of columns in fold change matrix:
                need same number of columns as number of groups.'))
        }
        if(nrow(fold_changes) != length(transcripts)){
            stop(.makepretty('wrong number of rows in fold change matrix: need
                same number of rows as number of simulated transcripts. see
                count_transcripts to find out that number.'))
        }
    }
}

.write_info = function(extras, transcripts, num_reps, fold_changes, outdir,
    group_ids, counts_matrix){

    if(!('transcriptid' %in% names(extras))){
        extras$transcriptid = names(transcripts)
    }

    if(is.numeric(fold_changes)){
        sim_info = data.frame(transcriptid=extras$transcriptid,
            foldchange=fold_changes, DEstatus=fold_changes!=1)
    }else{
        fcv = apply(fold_changes, 1, function(x){
            paste(x, collapse=';')
        })
        DEstatus = rowSums(fold_changes) != ncol(fold_changes)
        sim_info = data.frame(transcriptid=extras$transcriptid,
            foldchange=fcv, DEstatus=DEstatus)
    }

    write.table(sim_info, row.names=FALSE, quote=FALSE, sep="\t",
            file=paste0(outdir, '/sim_tx_info.txt'))

    rep_info = data.frame(
        rep_id=paste0('sample_', sprintf('%02d', 1:sum(num_reps))),
        group=group_ids, lib_sizes=extras$lib_sizes)

    write.table(rep_info, row.names=FALSE, quote=FALSE, sep='\t',
                file=paste0(outdir, '/sim_rep_info.txt'))

    rownames(counts_matrix) <- names(transcripts)
    colnames(counts_matrix) <- rep_info$rep_id
    save(counts_matrix, file=paste0(outdir, '/sim_counts_matrix.rda'))
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
#' @param outdir character, path to folder where simulated reads should be
#'   written, with *no* slash at the end. By default, reads are
#'   written to current working directory.
#' @param num_reps How many biological replicates should be in each group? The
#'   length \code{num_reps} determines how many groups are in the experiment.
#'   For example, \code{num_reps = c(5,6,5)} specifies a 3-group experiment with
#'   5 samples in group 1, 6 samples in group 2, and 5 samples in group 3.
#'   Defaults to a 2-group experiment with 10 reps per group (i.e.,
#'   \code{c(10,10)}).
#' @param reads_per_transcript baseline mean number of reads to simulate
#'   from each transcript. Can be an integer, in which case this many reads
#'   are simulated from each transcript, or an integer vector whose length
#'   matches the number of transcripts in \code{fasta}. Default 300. You can
#'   also leave \code{reads_per_transcript} empty and set \code{meanmodel=TRUE}
#'   to draw baseline mean numbers from a model based on transcript length.
#' @param size the negative binomial \code{size} parameter (see
#'   \code{\link{NegBinomial}}) for the number of reads drawn per transcript.
#'   It can be a matrix (where the user can specify the size parameter per
#'   transcript, per group), a vector (where the user can specify the size per
#'   transcript, perhaps relating to reads_per_transcript), or a single number,
#'   specifying the size for all transcripts and groups.
#'   If left NULL, defaults to \code{reads_per_transcript * fold_changes / 3}.
#'   Negative binomial variance is mean + mean^2 / size.
#' @param fold_changes Matrix specifying multiplicative fold changes
#'   between groups. There is no default, so you must provide this argument.
#'   In real data sets, lowly-expressed transcripts often show high fold
#'   changes between groups, so this can be kept in mind when setting
#'   \code{fold_changes} and \code{reads_per_transcript}. This argument must
#'   have the same number of columns as there are groups as
#'   specified by \code{num_reps}, and must have the same number of rows as
#'   there are transcripts in \code{fasta}. A fold change of X in matrix entry
#'   i,j means that for replicate j, the baseline mean number of reads
#'   (reads_per_transcript[i]) will be multiplied by X. Note that the
#'   multiplication happens before the negative binomial value
#'   (for the number of reads that *actually will* be
#'   drawn from transcript i, for replicate j) is drawn. This argument is
#'   ignored if \code{length(num_reps)} is 1 (meaning you only have 1 group in
#'   your simulation).
#' @param paired If \code{TRUE}, paired-end reads are simulated; else
#'   single-end reads are simulated. Default \code{TRUE}
#' @param reportCoverage whether to write out coverage information to
#'   \code{sample_coverages.rda} file in the \code{outdir}.
#'   defaults to \code{FALSE}
#' @param ... any of several other arguments that can be used to add nuance
#'   to the simulation. See details.
#'
#' @details Reads can either be simulated from a FASTA file of transcripts
#'   (provided with the \code{fasta} argument) or from a GTF file plus DNA
#'   sequences (provided with the \code{gtf} and \code{seqpath} arguments).
#'   Simulating from a GTF file and DNA sequences may be a bit slower: it took
#'   about 6 minutes to parse the GTF/sequence files for chromosomes 1-22, X,
#'   and Y in hg19.
#'
#'   Several optional parameters can be passed to this function to adjust the
#'   simulation. The options are:
#'
#'   \itemize{
#'   \item \code{readlen}: read length. Default 100.
#'   \item \code{lib_sizes}: Library size factors for the biological replicates.
#'   \code{lib_sizes} should have length equal to the total number of
#'   replicates in the experiment, i.e., \code{sum(num_reps)}. For each
#'   replicate, once the number of reads to simulate from each transcript for
#'   that replicate is known, all read numbers across all transcripts from that
#'   replicate are multiplied by the corresponding entry in \code{lib_sizes}.
#'   \item \code{distr} One of 'normal', 'empirical', or 'custom', which
#'   specifies the distribution from which to draw RNA fragment lengths. If
#'   'normal', draw fragment lengths from a normal distribution. You can provide
#'   the mean of that normal distribution with \code{fraglen} (defaults to 250)
#'   and the standard deviation of that normal distribution with \code{fragsd}
#'   (defaults to 25). You can provide a single number for each, or a
#'   vector with length equal to the total number of samples.
#'   If 'empirical', draw fragment lengths
#'   from a fragment length distribution estimated from a real data set. If
#'   'custom', draw fragment lengths from a custom distribution, which you can
#'   provide as the \code{custdens} argument. \code{custdens} should be a
#'   density fitted using \code{\link{logspline}}.
#'   \item \code{error_model}: The error model can be one of:
#'     \itemize{
#'     \item \code{'uniform'}: errors are distributed uniformly across reads.
#'     You can also provide an \code{'error_rate'} parameter, giving the overall
#'     probability of making a sequencing error at any given nucleotide. This
#'     error rate defaults to 0.005.
#'     \item \code{'illumina4'} or \code{'illumina5'}: Empirical error models.
#'     See \code{?add_platform_error} for more information.
#'     \item \code{'custom'}: A custom error model you've estimated from an
#'     RNA-seq data set using \code{GemErr}. See \code{?add_platform_error}
#'     for more info. You will need to provide both \code{model_path} and
#'     \code{model_prefix} if using a custom error model. \code{model_path} is
#'     the output folder you provided to \code{build_error_model.py}. This path
#'     should contain either two files suffixed _mate1 and _mate2, or a file
#'     suffixed _single. \code{model_prefix} is the 'prefix' argument you
#'     provided to \code{build_error_model.py} and is whatever comes before the
#'     _mate1/_mate2 or _single files in \code{model_path}.
#'     }
#'   \item \code{bias} One of 'none', 'rnaf', or 'cdnaf'. 'none'
#'   represents uniform fragment selection (every possible fragment in a
#'   transcript has equal probability of being in the experiment); 'rnaf'
#'   represents positional bias that arises in protocols using RNA
#'   fragmentation, and 'cdnaf' represents positional bias arising in protocols
#'   that use cDNA fragmentation (Li and Jiang 2012). Using the 'rnaf' model,
#'   coverage is higher in the middle of the transcript and lower at both ends,
#'   and in the 'cdnaf' model, coverage increases toward the 3' end of the
#'   transcript. The probability models used come from Supplementary Figure S3
#'   of Li and Jiang (2012). Defaults to 'none' if you don't provide this.
#'   \item \code{gcbias} list indicating which samples to add GC bias to, and
#'   from which models. Should be the same length as \code{sum(num_reps)};
#'   entries can be either numeric or of class \code{loess}. A numeric entry of
#'   0 indicates no GC bias. Numeric entries 1 through 7 correspond to the
#'   7 empirical GC models that ship with Polyester, estimated from GEUVADIS
#'   HapMap samples NA06985, NA12144, NA12776, NA18858, NA20542, NA20772,
#'   and NA20815, respectively. The code used to derive the empirical GC models
#'   is available at
#'   \url{https://github.com/alyssafrazee/polyester/blob/master/make_gc_bias.R}.
#'   A loess entry should be a loess prediction model
#'   that takes a GC content percent value (between 0 and 1) a transcript's
#'   deviation from overall mean read count based on that GC value. Counts for
#'   each replicate will be adjusted based on the GC bias model specified for
#'   it. Numeric and loess entries can be mixed. By default, no bias is
#'   included.
#'   \item \code{frag_GC_bias} Either a matrix of dimensions 101 x \code{sum(num_reps)}
#'   or 'none'. The default is 'none'. If specified, the matrix contains the probabilities
#'   (a number in the range [0,1]) that a fragment will appear in the output given its GC content.
#'   The first row corresponds to a fragment with GC content of 0 percent, the second row
#'   1 percent, the third row 2 percent, etc., and the last row 100 percent.
#'   The columns correspond to different probabilites for each sample. Internally,
#'   a coin flip (a Bernoulli trial) determines if each fragment is kept, depending on its GC content.
#'   Note that the final library size will depend on the elements of the matrix, and
#'   it might make sense to scale up the \code{lib_size} of the samples with low
#'   probabilites in the matrix in the range of the transcriptome GC content distribution.
#'   Note that the \code{count_matrix} written to \code{outdir} contains the counts before
#'   applying fragment GC bias.
#'   \item \code{strand_specific} defaults to \code{FALSE}, which means fragments are
#'   generated with equal probability from both strands of the transcript sequence.
#'   set to \code{TRUE} for strand-specific simulation (1st read forward strand,
#'   2nd read reverse strand with respect to transcript sequence).
#'   \item \code{meanmodel}: set to TRUE if you'd like to set
#'   \code{reads_per_transcripts} as a function of transcript length. We
#'   fit a linear model regressing transcript abundance on transcript length,
#'   and setting \code{meanmodel=TRUE} means we will use transcript lengths
#'   to draw transcript abundance based on that linear model. You can see our
#'   modeling code at \url{http://htmlpreview.github.io/?https://github.com/alyssafrazee/polyester_code/blob/master/length_simulation.html}
#'   \item \code{write_info}: set to FALSE if you do not want files of
#'   simulation information written to disk. By default, transcript fold
#'   changes and expression status, replicate library sizes and group
#'   identifiers, and an R data object of the counts matrix (before
#'   application of fragment GC bias) are written to \code{outdir}.
#'   \item \code{seed}: specify a seed (e.g. \code{seed=142} or some other
#'   integer) to set before randomly drawing read numbers, for reproducibility.
#'   \item \code{transcriptid}: optional vector of transcript IDs to be written
#'   into \code{sim_info.txt} and used as transcript identifiers in the output
#'   fasta files. Defaults to \code{names(readDNAStringSet(fasta))}. This
#'   option is useful if default names are very long or contain special
#'   characters.
#'   \item \code{gzip}: pass \code{gzip=TRUE} to write gzipped fasta files as
#'   output (by default, fasta output files are not compressed when written to
#'   disk).
#'   \item \code{exononly}: (passed to \code{\link{seq_gtf}}) if \code{TRUE}
#'   (as it is by default), only create transcript sequences from the features
#'   labeled \code{exon} in \code{gtf}.
#'   \item \code{idfield}: (passed to \code{\link{seq_gtf}})in the
#'   \code{attributes} column of \code{gtf}, what is the name of the field
#'   identifying transcripts? Should be character. Default
#'   \code{"transcript_id"}.
#'   \item \code{attrsep}: (passed to \code{\link{seq_gtf}}) in the
#'   \code{attributes} column of \code{gtf}, how are attributes separated?
#'   Default \code{"; "}.
#'   }
#'
#' @references
#'   't Hoen PA, et al (2013): Reproducibility of high-throughput mRNA and
#'   small RNA sequencing across laboratories. Nature Biotechnology 31(11):
#'   1015-1022.
#'
#'   Li W and Jiang T (2012): Transcriptome assembly and isoform expression
#'   level estimation from biased RNA-Seq reads. Bioinformatics 28(22):
#'   2914-2921.
#'
#'   McElroy KE, Luciani F and Thomas T (2012): GemSIM: general,
#'   error-model based simulator of next-generation sequencing data. BMC
#'   Genomics 13(1), 74.
#'
#' @return No return, but simulated reads and a simulation info file are written
#'   to \code{outdir}. Note that reads are written out transcript by transcript
#'   and so need to be shuffled if used as input to quantification algorithms such
#'   as eXpress or Salmon.
#'
#' @export
#'
#' @examples \donttest{
#'   ## simulate a few reads from chromosome 22
#'
#'   fastapath = system.file("extdata", "chr22.fa", package="polyester")
#'   numtx = count_transcripts(fastapath)
#'   set.seed(4)
#'   fold_change_values = sample(c(0.5, 1, 2), size=2*numtx,
#'      prob=c(0.05, 0.9, 0.05), replace=TRUE)
#'   fold_changes = matrix(fold_change_values, nrow=numtx)
#'   library(Biostrings)
#'   # remove quotes from transcript IDs:
#'   tNames = gsub("'", "", names(readDNAStringSet(fastapath)))
#'
#'   simulate_experiment(fastapath, reads_per_transcript=10,
#'      fold_changes=fold_changes, outdir='simulated_reads',
#'      transcriptid=tNames, seed=12)
#'}
#'
simulate_experiment = function(fasta=NULL, gtf=NULL, seqpath=NULL,
    outdir='.', num_reps=c(10,10), reads_per_transcript=300, size=NULL,
    fold_changes, paired=TRUE, reportCoverage=FALSE, ...){

    extras = list(...)

    # validate extra arguments/set sane defaults
    extras = .check_extras(extras, paired, total.n=sum(num_reps))

    # read in the annotated transcripts to sequence from
    if(!is.null(fasta) & is.null(gtf) & is.null(seqpath)){
        transcripts = readDNAStringSet(fasta)
    }else if(is.null(fasta) & !is.null(gtf) & !is.null(seqpath)){
        message('parsing gtf and sequences...')
        # parse out any extra seq_gtf arguments from the ... args
        if('exononly' %in% names(extras)){
            exononly = extras$exononly
        }else{
            exononly = TRUE
        }
        if('idfield' %in% names(extras)){
            idfield = extras$idfield
        }else{
            idfield = 'transcript_id'        
        }
        if('attrsep' %in% names(extras)){
            attrsep = extras$attrsep
        }else{
            attrsep = '; '
        }
        transcripts = seq_gtf(gtf, seqpath, feature='transcript',
            exononly=exononly, idfield=idfield, attrsep=attrsep)
        message('done parsing')
    }else{
        stop('must provide either fasta or both gtf and seqpath')
    }

    # check fold change matrix dimensions:
    .check_fold_changes(fold_changes, num_reps, transcripts)

    # get baseline means for each group, incl. fold changes:
    if('meanmodel' %in% names(extras)){
        b0 = -3.0158
        b1 = 0.8688
        sigma = 4.152
        logmus = b0 + b1*log2(width(transcripts)) + rnorm(length(transcripts),0,sigma)
        reads_per_transcript = 2^logmus-1
        reads_per_transcript = pmax( reads_per_transcript, 1e-6 )
    }

    if(length(num_reps) == 1){
        fold_changes = matrix(rep(1, length(transcripts)))
    } else if(length(num_reps) == 2) {
        # This means fold_changes is a numeric vector, per the check function
        if(length(reads_per_transcript) == 1) {
            basemeans = matrix(reads_per_transcript, ncol=2, nrow=length(transcripts))
            basemeans[,2] = fold_changes * basemeans[,1]
        } else if(length(reads_per_transcript) == length(transcripts)){
            basemeans = matrix(c(reads_per_transcript, reads_per_transcript), nrow=length(reads_per_transcript))
            basemeans = basemeans*fold_changes
        } else {
            stop('reads_per_transcript is the wrong length.')
        }
    } else {
        basemeans = reads_per_transcript * fold_changes
    }

    if(is.null(size)){
        size = basemeans / 3
    }else if(class(size) == 'numeric'){
        if(is.matrix(basemeans)){
            num_rows = nrow(basemeans)
            num_cols = ncol(basemeans)
        } else {
            num_rows = length(basemeans)
            num_cols = 1
        }
        size = matrix(size, nrow=num_rows, ncol=num_cols)
    }else if(class(size) == 'matrix'){
        if(!is.matrix(basemeans)){
            stop('If you provide a matrix for size, you also need a matrix for reads_per_transcript.')
        }
        stopifnot(nrow(size) == nrow(basemeans))
        stopifnot(ncol(size) == ncol(basemeans))
    }else{
        stop('size must be a number, numeric vector, or matrix.')
    }

    # create matrix of transcripts & number of reads to simulate
    if('seed' %in% names(extras)){
        set.seed(extras$seed)
    }
    group_ids = rep(1:length(num_reps), times=num_reps)
    numreadsList = vector("list", sum(num_reps))
    numreadsList = lapply(1:sum(num_reps), function(i){
        group_id = group_ids[i]
        NB(as.matrix(basemeans)[,group_id], as.matrix(size)[,group_id])
    })
    readmat = matrix(unlist(numreadsList), ncol=sum(num_reps))
    readmat = t(extras$lib_sizes * t(readmat))
    if('gcbias' %in% names(extras)){
        stopifnot(length(extras$gcbias) == sum(num_reps))
        gcclasses = unique(sapply(extras$gcbias, class))
        if(sum(gcclasses %in% c('numeric', 'loess')) < length(extras$gcbias)){
            stop(.makepretty('gc bias parameters must be integers 0 through 7
                or loess objects'))
        }
        if(any(extras$gcbias!=0)){
            readmat = add_gc_bias(readmat, extras$gcbias, transcripts)
        }
    }

    # prep output directory
    sysoutdir = gsub(' ', '\\\\ ', outdir)
    if(.Platform$OS.type == 'windows'){
        shell(paste('mkdir', sysoutdir))
    }else{
        system(paste('mkdir -p', sysoutdir))
    }

    # do the actual sequencing
    sgseq(readmat, transcripts, paired, outdir, extras, reportCoverage)

    # write out simulation information, if asked for:
    if(!('write_info' %in% names(extras))){
        write_info=TRUE
    }

    if(write_info){
      # save the *unbiased* counts matrix: the counts for each
      # transcript and each sample *before* fragment GC bias is applied
      counts_matrix <- readmat
      .write_info(extras, transcripts, num_reps, fold_changes, outdir,
                  group_ids, counts_matrix)
    }
}

# simulate_experiment(fastapath, reads_per_transcript=10, outdir='~/Desktop/tmp', num_reps=1)
# simulate_experiment(fastapath, reads_per_transcript=10, outdir='~/Desktop/tmp', num_reps=c(1,1), fold_changes=fc2[,2:3])
