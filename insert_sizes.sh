#!/bin/sh
#$ -cwd -l mf=10G,h_vmem=10G

## estimate insert sizes from each of these BAM files using picard

for sample in NA06985 NA18858 NA20772 NA12144 NA20815 NA12776 NA20542
do
    java -Xmx2g -jar picard-tools-1.121/CollectInsertSizeMetrics.jar HISTOGRAM_FILE=${sample}_histogram INPUT=${sample}_accepted_hits.bam OUTPUT=${sample}_metrics
done