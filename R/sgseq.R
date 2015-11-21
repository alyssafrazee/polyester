## internal sequencing function

sgseq = function(readmat, transcripts, paired, outdir, extras){

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
}
