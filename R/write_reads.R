#' write sequencing reads to disk
#'
#' given a DNAStringSet representing simulated sequencing reads, write FASTA 
#'   files to disk representing the simulated reads.
#' @param reads DNAStringSet representing sequencing reads
#' @param fname file path/prefix specifying where sequencing reads should be
#'   written. Should not contain ".fasta" (this is appended automatically). 
#' @param readlen maximum length of the reads in \code{reads}.
#' @param paired If \code{TRUE}, reads are assumed to be in pairs: i.e., read 1 
#'   and read 2 in \code{reads} are the left and right mate (respectively) of a
#'   read pair; same with read 3 and read 4, etc. The odd-numbered reads are 
#'   written to \code{fname_1.fasta} and the even-numbered
#'   reads are written to \code{fname_2.fasta}. If \code{FALSE}, reads are 
#'   assumed to be single-end and just one file, \code{fname.fasta}, is written.
#' @param gzip If \code{TRUE}, gzip the output fasta files.
#' @export
#' @details The \code{\link{get_reads}} function returns a DNAStringSet object
#'   representing sequencing reads that can be directly passed to 
#'   \code{write_reads}. If output other than that from \code{get_reads} is used
#'   and \code{paired} is \code{TRUE}, make sure \code{reads} is ordered 
#'   properly (i.e., that mate pairs appear together and that the left mate 
#'   appears first).
#' @return No return, but FASTA file(s) containing the sequences in \code{reads} 
#'   are written to \code{fname.fasta} (if \code{paired} is FALSE) or 
#'   \code{fname_1.fasta} and \code{fname_2.fasta} if \code{paired} is TRUE.
#' @seealso \code{\link{get_reads}}
#' @examples
#'   library(Biostrings)
#'   data(srPhiX174) # pretend srPhiX174 represents a DNAStringSet of *reads*
#'   readlen = unique(width(srPhiX174)) #35
#'   write_reads(srPhiX174, fname='./srPhiX174', readlen=readlen, paired=FALSE,
#'       gzip=FALSE)
#'
write_reads = function(reads, fname, readlen, paired=TRUE, gzip){
    compress = ifelse(is.null(gzip), FALSE, gzip)
    if(paired){
        lefts = reads[seq(1, length(reads), by=2)]
        rights = reads[seq(2, length(reads), by=2)]
        names(lefts) = paste0('read', 1:length(lefts), '/', names(lefts))
        names(rights) = paste0('read', 1:length(rights), '/', names(rights))
        left_filepath = paste0(fname, '_1.fasta')
        right_filepath = paste0(fname, '_2.fasta')
        if(compress){
            left_filepath = paste0(left_filepath, '.gz')
            right_filepath = paste0(right_filepath, '.gz')
        }
        writeXStringSet(lefts, filepath=left_filepath, 
            format="fasta", width=readlen, compress=compress)
        writeXStringSet(rights, filepath=right_filepath, 
            format="fasta", width=readlen, compress=compress)
    }else{
        outf = paste0(fname, '.fasta')
        if(compress){
            outf = paste0(outf, '.gz')
        }
        names(reads) = paste0('read', 1:length(reads), '/', names(reads))
        writeXStringSet(reads, filepath=outf, format="fasta", width=readlen,
            compress=compress)
    }
}
