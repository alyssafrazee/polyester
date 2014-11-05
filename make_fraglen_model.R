## get the empirical fragment length distributions from GEUVADIS samples
## (maybe 1 from each lab?)
## save empirical distribution as a dataset in Polyester package
## I worked in a folder with the same root directory as the "polyester" R package

library(ballgown)
load('fpkm.rda') #link: http://files.figshare.com/1625419/fpkm.rda

pd = pData(fpkm)

library(RCurl)
QCurl = getURL('https://raw.githubusercontent.com/alyssafrazee/ballgown_code/master/GEUVADIS_preprocessing/GD667.QCstats.masterfile.txt')
qcstats = read.table(textConnection(QCurl), sep='\t', header=TRUE)
pd$lab = qcstats$SeqLabNumber[match(pd$SampleID, rownames(qcstats))]

# samples to analyze: 1 from each lab
set.seed(14)
lab1id = sample(pd$dirname[pd$lab == 'F1'], size=1)
lab2id = sample(pd$dirname[pd$lab == 'F2'], size=1)
lab3id = sample(pd$dirname[pd$lab == 'F3'], size=1)
lab4id = sample(pd$dirname[pd$lab == 'F4'], size=1)
lab5id = sample(pd$dirname[pd$lab == 'F5'], size=1)
lab6id = sample(pd$dirname[pd$lab == 'F6'], size=1)
lab7id = sample(pd$dirname[pd$lab == 'F7'], size=1)
selected = c(lab1id, lab2id, lab3id, lab4id, lab5id, lab6id, lab7id)
selected
# "NA06985" "NA18858" "NA20772" "NA12144" "NA20815" "NA12776" "NA20542"
# 3 CEU, 1 YRI, 3 TSI

# download the bam files from each of those samples:
for(s in selected){
    system(paste0('wget http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-6/',
        s, '_accepted_hits.bam'), wait=TRUE)
} #UGH this takes so long, about 20-40 minutes per file, depending on connection speed.
# could do in parallel
# (but more difficult to script) 

# run "CollectInsertSizeMetrics", from Picard:
system('sh insert_sizes.sh', wait=TRUE)
## SCRIPT CONTENTS:
# #!/bin/sh
# for sample in NA06985 NA18858 NA20772 NA12144 NA20815 NA12776 NA20542
# do
#     java -Xmx2g -jar picard-tools-1.121/CollectInsertSizeMetrics.jar HISTOGRAM_FILE=${sample}_histogram INPUT=${sample}_accepted_hits.bam OUTPUT=${sample}_metrics
# done
# takes about 25 minutes

# fit a distribution to the insert sizes:
library(IRanges)
library(logspline)
frag_sizes = Rle()
set.seed(212)
for(s in selected){
    ind = which(selected == s)
    metrics = file(paste0(s, '_metrics'))
    info = strsplit(readLines(metrics), '\t')
    close(metrics)
    dat = t(as.data.frame(info[-c(1:11, length(info))])) #extract sizes/counts
    rownames(dat) = NULL
    dat = matrix(as.numeric(dat), nrow=nrow(dat), ncol=ncol(dat))
    dat = as.data.frame(dat)
    names(dat) = c('size', 'count')
    allnums = Rle(values=dat$size, lengths=dat$count)
    samplednums = sort(sample(allnums, 1e5))
    frag_sizes = sort(c(frag_sizes, samplednums))
}
empirical_density = logspline(as.integer(frag_sizes), 
    lbound=min(frag_sizes), ubound=max(frag_sizes))
save(empirical_density, file='../polyester/data/empirical_density.rda', compress='xz')

# we can then draw from the empirical density in the fragmentation function!

# plot for the paper:
d = density(as.integer(frag_sizes))
pdf('empirical_density.pdf')
    x = seq(75, max(frag_sizes), by=1)
    plot(x, dnorm(x, 250, 25), col='blue', 
        type='l', lwd=2, xlab='Fragment Length', ylab='Density')
    lines(d, col='red', lwd=2)
    #with(d, polygon(x=c(min(frag_sizes), d$x, max(frag_sizes)), 
    #                y=c(0, d$y, 0), col="gray"))
    legend('topright', col=c('blue', 'red'), c('N(250, 25)', 'Empirical'), 
        lwd=c(2,2))
    title('Normal and empirical fragment length distributions')
dev.off()
