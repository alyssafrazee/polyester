#' @title Get transcript sequences from GTF file and sequence info
#'
#' @description Given a GTF file (for transcript structure) and DNA sequences, return a DNAStringSet
#' of transcript sequences
#' @param gtf path to GTF file
#' @param seqpath path to folder containing one FASTA file (\code{.fa} extension) for each 
#' chromosome in \code{gtf}
#' @param exononly if \code{TRUE} (as it is by default), only create transcript sequences from the 
#' features labeled \code{exon} in \code{gtf}
#' @param idfield in the \code{attributes} column of \code{gtf}, what is the name of the field 
#' identifying transcripts? Should be character. Default \code{"transcript_id"}.
#' @param attrsep in the \code{attributes} column of \code{gtf}, how are attributes separated? 
#' Default \code{"; "}.
#' @export
#' @references \url{http://www.ensembl.org/info/website/upload/gff.html}
#' @return DNAStringSet containing transcript sequences, with names corresponding to \code{idfield}
#' in \code{gtf}
seq_gtf = function(gtf, seqpath, exononly=TRUE, idfield="transcript_id", attrsep="; "){
    
    gtf_dat = read.table(gtf, sep = "\t", as.is = TRUE, quote = "", header = FALSE, 
        comment.char = "#", nrows = -1, colClasses = c("character", "character", "character", 
            "integer", "integer", "character", "character", "character", "character"))

    colnames(gtf_dat) = c("seqname", "source", "feature", "start", "end", "score", "strand", 
        "frame", "attributes")

    stopifnot(!any(is.na(gtf_dat$start)), !any(is.na(gtf_dat$end)))

    if(exononly){
        gtf_dat = gtf_dat[gtf_dat[,3]=="exon",]
    }

    chrs = unique(gtf_dat$seqname)
    fafiles = list.files(seqpath)
    if(!(all(paste0(chrs, '.fa') %in% fafiles))){
        stop("all chromosomes in the GTF file must have .fa files in seqpath")
    }

    seqlist = lapply(chrs, function(chr){
        dftmp = gtf_dat[gtf_dat[,1]==chr,]
        fullseq = readDNAStringSet(paste0(seqpath, '/', chr, '.fa'))
        these_seqs = subseq(rep(fullseq, times=nrow(dftmp)), start=dftmp$start, end=dftmp$end)
        names(these_seqs) = getAttributeField(dftmp$attributes, idfield, attrsep=attrsep)
        these_seqs
    })

    full_list = do.call(c, seqlist)
    split_list = split(full_list, names(full_list))
    DNAStringSet(lapply(split_list, unlist)) #took 340 sec on whole human transcriptome hg19
}

