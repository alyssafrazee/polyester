# count the transcripts in a FASTA or GTF file
#'
#' determine how many transcripts are annotated in a FASTA or GTF file
#' @param f character, path to a file in FASTA or GTF format
#' @param fasta TRUE if \code{f} is a fasta file; FALSE if \code{f} is a GTF
#'   file
#' @param identifier if \code{f} is a GTF file, how are transcripts identified
#'   in the attributes field (9th column) of the file? Default
#'   \code{transcript_id}.
#' @param attrsep if \code{f} is a GTF file, how are attributes separated in
#'   the attributes field (9th column) of the file? Default "; ".
#' @export
#' @return Number of transcripts annotated in \code{f}
#' @examples
#' fastapath = system.file("extdata", "chr22.fa", package="polyester")
#' count_transcripts(fastapath) #918
#'
count_transcripts = function(f, fasta=TRUE, identifier='transcript_id',
    attrsep="; "){
    flength = nchar(f)
    if(fasta){
        if(substr(f, flength-1, flength) != 'fa' &
            substr(f, flength-4, flength) != 'fasta'){
                stop(.makepretty('if f is a FASTA file, f must have extension
                    ".fa" or ".fasta"'))
        }
        return(length(readDNAStringSet(f)))
    }

    # else, it's a GTF file
    ext = substr(f, flength-2, flength)
    if(ext != 'gtf' & ext != 'gff'){
        stop('if f is a GTF file, f must have extension ".gtf" or ".gff"')
    }
    gtf_info = read.table(f, sep="\t", as.is=TRUE, quote="",
        header=FALSE, comment.char="#", nrows= -1,
        colClasses=c(rep('NULL', 8), 'character'))
    transcripts = getAttributeField(gtf_info[,1], field=identifier, attrsep)
    return(length(unique(transcripts)))
}
