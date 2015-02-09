## internal sequencing function

sgseq = function(readmat, transcripts, paired, outdir, extras){

    for(i in 1:ncol(readmat)){

        tObj = rep(transcripts, times=readmat[,i])
        
        #get fragments
        tFrags = generate_fragments(tObj, extras$fraglen, extras$fragsd, 
            extras$distr, extras$custdens, extras$bias)

        #reverse_complement some of those fragments
        rctFrags = reverse_complement(tFrags, library_type, strand_error_rate)

        #get reads from fragments
        reads = get_reads(rctFrags, extras$readlen, paired)

        #add sequencing error
        if(extras$error_model == 'uniform'){
            errReads = add_error(reads, extras$error_rate)
        }else if(error_model == 'custom'){
            errReads = add_platform_error(reads, 'custom', paired, extras$path)
        }else{
            errReads = add_platform_error(reads, extras$error_model, paired)
        }

        #write read pairs
        write_reads(errReads, readlen=extras$readlen, 
            fname=paste0(outdir, '/sample_', sprintf('%02d', i)), paired=paired)
    }
}
