#' Internal sequencing function
#' @importFrom parallel mclapply
sgseq = function(readmat, transcripts, paired, outdir, extras, reportCoverage=FALSE, ncores=1L){
  #report theoretically perfect coverage if reportCoverage=TRUE, will write a file
  if(reportCoverage==TRUE){
    templates = unique(transcripts)
    coverage_matrices = list()
    for(i in 1:length(templates)){coverage_matrices = c(coverage_matrices, list(matrix(0, ncol=dim(readmat)[2], width(templates)[i])))}
    names(coverage_matrices) = names(templates)
  }
  
  null <- parallel::mclapply(seq_len(ncol(readmat)), function(i) {
    ##$ begin small chunk regarding fragment GC bias or not
    if (is.matrix(extras$frag_GC_bias)) {
      frag_GC_bias <- extras$frag_GC_bias[,i]
    } else {
      frag_GC_bias <- 'none'
    }
    ### end 
    tObj = rep(transcripts, times=readmat[,i])
    iterations = ceiling(length(tObj) / 1e6L)
    offset = 1L
    for(iteration in seq_len(iterations)) {
      tSubset = tObj[offset:min(offset+999999L, length(tObj))] ## corrected value of integer added to offset to avoid duplicating reads
      tFrags = generate_fragments(tSubset, extras$fraglen[i], extras$fragsd[i],
                                  extras$readlen, extras$distr, extras$custdens,
                                  extras$bias, frag_GC_bias)

      if (!extras$strand_specific) {
        #reverse_complement some of those fragments
        tFrags = reverse_complement(tFrags)
      }
      
      #get reads from fragments
      reads = get_reads(tFrags, extras$readlen, paired)

      if(reportCoverage==TRUE){
            read_info = unique(names(reads))
            read_info_split = strsplit(read_info, ";mate1:|;mate2:")
            read_info_matrix = matrix(unlist(read_info_split), ncol=3, byrow=T)
            for(j in 1:dim(read_info_matrix)[1]){
                  read = read_info_matrix[j,]
                  target = which(names(coverage_matrices)==read[1])
                  # ML: changing these to strsplit (str_split requires stringr depends or imports)
                  coords1 = unlist(strsplit(read[2], "-"))
                  coords2 = unlist(strsplit(read[3], "-"))
                  coverage_matrices[[target]][coords1[1]:coords1[2],i]=coverage_matrices[[target]][coords1[1]:coords1[2],i]+1
                  coverage_matrices[[target]][coords2[1]:coords2[2],i]=coverage_matrices[[target]][coords2[1]:coords2[2],i]+1
                  save(coverage_matrices, file=file.path(outdir, 'sample_coverages.rda') )
            }
      }

      #add sequencing error
      if(extras$error_model == 'uniform'){
          errReads = add_error(reads, extras$error_rate)
      }else if(extras$error_model == 'custom'){
          errReads = add_platform_error(reads, 'custom', paired, extras$path)
      }else{
          errReads = add_platform_error(reads, extras$error_model, paired)
      }

      #write read pairs
      write_reads(errReads, readlen=extras$readlen,
          fname=sprintf('%s/sample_%04d', outdir,i), paired=paired,
          gzip=extras$gzip, offset=offset, shuffle = extras$shuffle)
      offset = offset + 1e6L
    }
  }, mc.cores = ifelse(ncol(readmat) >= ncores, ncores, ncol(readmat)))
  invisible(NULL)
}
