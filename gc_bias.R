## code for GC bias model implemented in Polyester

library(ballgown)
library(polyester)
library(Biostrings)
load('../polyester_paper/experiments/fpkm.rda') #transcript expression data. URL: http://files.figshare.com/1625419/fpkm.rda
load('../polyester_paper/experiments/sequences.rda') #transcript sequences

# here is the code used to make "sequences": (sequences.rda)
# this is the path to this index: ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz
# (this download will give you up to "Homo_sapiens")
seqpath = '/amber2/scratch/jleek/iGenomes-index/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/'
sequences = seq_gtf(gtfdf, seqpath)
names(sequences) = substr(names(sequences), 2, nchar(names(sequences))-1)
save(sequences, file='sequences.rda') 

# calculate GC content:
GC = letterFrequency(sequences, letters='GC', as.prob=TRUE)

# filter to reasonable expression
expressed = exprfilter(fpkm, cutoff=1)

# estimate transcript-level counts
txcounts = fpkm_to_counts(expressed, mean_rps=1e7)

# put counts in same order as sequence data
txcounts = txcounts[match(names(sequences), texpr(expressed,'all')$t_name),]
rownames(txcounts) = names(sequences)

# use the same 7 reps we used to estimate fragment distribution
# those reps are:
samples = c("NA06985", "NA18858", "NA20772", "NA12144", "NA20815", "NA12776", "NA20542")
cols_samples = ballgown:::ss(colnames(txcounts), pattern='\\.', slot=2)
counts = txcounts[,cols_samples %in% samples]
counts = counts[GC >= 0.25,] #remove some outliers driving loess
GC = GC[GC >= 0.25]
lcounts = log2(counts+1)

# this code creates Supplementary Figure 1 in the Polyester paper
pdf('GC_effects.pdf', height=12, width=6)
par(mfrow=c(4,2))
for(i in 1:7){
    sample_id = ballgown:::ss(colnames(counts)[i], pattern='\\.', slot=2)
    plot(GC, lcounts[,i] - mean(lcounts[,i]), col='#00000050', pch=19, 
        ylab='Difference from average count (log scale)', xlab='GC %', 
            main=paste0('Model ', i, ': Sample ', sample_id))
    o = order(GC)
    loessfit = loess(lcounts[,i] - mean(lcounts[,i]) ~ GC, span=0.3)
    lines(GC[o], predict(loessfit)[o], col='deeppink', lwd=3)
    legend('topright', lwd=2, col='deeppink', 'loess fit')
}
dev.off()


# This code saves the 7 loess models (deviation-wise) as data sets
lcounts = log2(counts+1)
for(s in samples){
    i = which(samples == s)
    fit = loess(lcounts[,i]-mean(lcounts[,i]) ~ GC, span=0.3)
    eval(parse(text=paste0('loessfit', i, ' <- fit')))
    save(list=paste0('loessfit',i), file=paste0('data/', s, '_GC.rda'), compress='xz')
}
