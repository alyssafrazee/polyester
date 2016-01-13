## internal sequencing function

sgseq = function(readmat, transcripts, paired, outdir, extras, reportCoverage=F){
      #report theoretically perfect coverage if reportCoverage=T, will write a file
      if(reportCoverage==T){
            templates = unique(transcripts)
            coverage_matrices = list()
            for(i in 1:length(templates)){coverage_matrices = c(coverage_matrices, list(matrix(0, ncol=dim(readmat)[2], width(templates)[i])))}
            names(coverage_matrices) = names(templates)
      }
      
  for(i in seq_len(ncol(readmat))) {
    tObj = rep(transcripts, times=readmat[,i])
    iterations = ceiling(length(tObj) / 1e6L)
    offset = 1L
    for(iteration in seq_len(iterations)) {
      tSubset = tObj[offset:min(offset+1e6L, length(tObj))]
      tFrags = generate_fragments(tSubset, extras$fraglen, extras$fragsd,
                                  extras$readlen, extras$distr, extras$custdens,
                                  extras$bias)      
      #reverse_complement some of those fragments
      rctFrags = reverse_complement(tFrags)

      #get reads from fragments
      reads = get_reads(rctFrags, extras$readlen, paired)

      if(reportCoverage==T){
            read_info = unique(names(reads))
            read_info_split = strsplit(read_info, ";mate1:|;mate2:")
            read_info_matrix = matrix(unlist(read_info_split), ncol=3, byrow=T)
            for(j in 1:dim(read_info_matrix)[1]){
                  read = read_info_matrix[j,]
                  target = which(names(coverage_matrices)==read[1])
                  coords1 = unlist(str_split(read[2], "-"))
                  coords2 = unlist(str_split(read[3], "-"))
                  coverage_matrices[[target]][coords1[1]:coords1[2],i]=coverage_matrices[[target]][coords1[1]:coords1[2],i]+1
                  coverage_matrices[[target]][coords2[1]:coords2[2],i]=coverage_matrices[[target]][coords2[1]:coords2[2],i]+1
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
          fname=paste0(outdir, '/sample_', sprintf('%02d', i)), paired=paired,
          gzip=extras$gzip, offset=offset)
      offset = offset + 1e6L
    }
  }
  save(coverage_matrices, file=paste0(outdir, '/sample_coverages.rda') )
  
}
