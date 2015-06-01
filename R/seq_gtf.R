#' @title Get transcript sequences from GTF file and sequence info
#'
#' @description Given a GTF file (for transcript structure) and DNA sequences, 
#'   return a DNAStringSet of transcript sequences
#' @param gtf one of path to GTF file, or data frame representing a canonical
#'   GTF file.
#' @param seqs one of path to folder containing one FASTA file (\code{.fa} 
#'   extension) for each chromosome in \code{gtf}, or named DNAStringSet
#'   containing one DNAString per chromosome in \code{gtf}, representing its 
#'   sequence. In the latter case, \code{names(seqs)} should contain the
#'   same entries as the \code{seqnames} (first) column of \code{gtf}.
#' @param feature one of \code{'transcript'} or \code{'exon'} (default 
#'   transcript), depending on desired return.
#' @param exononly if \code{TRUE} (as it is by default), only create transcript
#'   sequences from the features labeled \code{exon} in \code{gtf}.
#' @param idfield in the \code{attributes} column of \code{gtf}, what is the
#'   name of the field identifying transcripts? Should be character. Default
#'   \code{"transcript_id"}.
#' @param attrsep in the \code{attributes} column of \code{gtf}, how are
#'   attributes separated? Default \code{"; "}.
#' @export
#' @references \url{http://www.ensembl.org/info/website/upload/gff.html}
#' @return 
#'   If feature is \code{'transcript'}, DNAStringSet containing 
#'   transcript sequences, with names corresponding to \code{idfield} in
#'   \code{gtf}. If feature is \code{'exon'}, DNAStringSet containing exon
#'   sequences from \code{gtf}, named by exon location (chr, start, end, 
#'   strand).
#' @examples  \dontrun{
#'   library(Biostrings)
#'   system('wget https://www.dropbox.com/s/04i6msi9vu2snif/chr22seq.rda')
#'   load('chr22seq.rda')
#'   data(gtf_dataframe)
#'   chr22_processed = seq_gtf(gtf_dataframe, chr22seq)
#'}
seq_gtf = function(gtf, seqs, feature='transcript', exononly=TRUE, 
    idfield='transcript_id', attrsep="; "){

    feature = match.arg(feature, c('transcript', 'exon'))

    gtfClasses = c("character", "character", "character", "integer", 
        "integer", "character", "character", "character", "character")
    if(is.character(gtf)){
        # read transcript structure from file:s
        gtf_dat = read.table(gtf, sep="\t", as.is=TRUE, quote="", header=FALSE, 
            comment.char="#", nrows= -1, colClasses=gtfClasses)
    } else if(is.data.frame(gtf)){
        # do what we can to check whether gtf really does represent a 
        # canonical GTF
        stopifnot(ncol(gtf) == 9)
        if(!all(unlist(lapply(gtf, class)) == gtfClasses)){
            stop("one or more columns of gtf have the wrong class")
        }
        gtf_dat = gtf
        rm(gtf)
    } else {
        stop("gtf must be a file path or a data frame")
    }

    colnames(gtf_dat) = c("seqname", "source", "feature", "start", "end",
        "score", "strand", "frame", "attributes")
    stopifnot(!any(is.na(gtf_dat$start)), !any(is.na(gtf_dat$end)))

    if(exononly){
        gtf_dat = gtf_dat[gtf_dat[,3]=="exon",]
    }

    # makes sure all chromosomes are present:
    chrs = unique(gtf_dat$seqname)
    if(is.character(seqs)){
        fafiles = list.files(seqs)
        lookingFor = paste0(chrs, '.fa')
    } else {
        fafiles = names(seqs)
        lookingFor = chrs
    }
    if(!(all(lookingFor %in% fafiles))){
        stop("all chromosomes in gtf must have corresponding sequences in seqs")
    }
    
    seqlist = lapply(chrs, function(chr){
        dftmp = gtf_dat[gtf_dat[,1] == chr,]
        if(is.character(seqs)){
            fullseq = readDNAStringSet(paste0(seqs, '/', chr, '.fa'))
        } else {
            fullseq = seqs[which(names(seqs) == chr)]
        }
        if(feature == 'exon'){
            dftmp = dftmp[!duplicated(dftmp[,c(1,4,5,7)]),] #unique exons
        }
        these_seqs = subseq(rep(fullseq, times=nrow(dftmp)), 
            start=dftmp$start, end=dftmp$end)
        if(feature == 'transcript'){
            names(these_seqs) = getAttributeField(dftmp$attributes, idfield, 
                attrsep=attrsep)
            if(substr(names(these_seqs)[1],1,1) == '"'){
                x = names(these_seqs)
                names(these_seqs) = substr(x, 2, nchar(x)-1)
            }
        }else{
            names(these_seqs) = paste0(dftmp[,1], ':', dftmp[,4], '-', 
                dftmp[,5], '(', dftmp[,7], ')')
        }
        revstrand = which(dftmp$strand == '-')
        these_seqs[revstrand] = reverseComplement(these_seqs[revstrand])
        these_seqs
    })

    full_list = do.call(c, seqlist)

    if(feature == 'exon'){
        return(full_list)
    }else{
        split_list = split(full_list, names(full_list))
        return(DNAStringSet(lapply(split_list, unlist)))
    }
}

