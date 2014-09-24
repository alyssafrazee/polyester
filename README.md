# Introduction

Polyester is an R package designed to simulate an RNA sequencing experiment. Given a set of annotated transcripts, polyester will simulate the steps of an RNA-seq experiment (fragmentation, reverse-complementing, and sequencing) and produce files containing simulated RNA-seq reads. Simulated reads can be analyzed using any of several downstream analysis tools. 

In particular, Polyester was designed to simulate a case/control experiment with biological replicates. Users are able to set differential transcript expression between cases and controls. This allows users to create datasets with known differential expression, which means they can the accuracy of statistical methods for differential expression detection.

Polyester was developed with several specific features in mind:  
* Simulation of differential expression at the transcript level
* Ability to set differential expression signal strength
* Simulation of small datasets, since large RNA-seq datasets can require lots of time and computing resources to analyze
* Generation of raw RNA-seq reads (as opposed to read alignments or transcript-level abundance estimates)
* Transparency/open-source code

# Installation

Start R and run:


```r
source("http://bioconductor.org/biocLite.R")
biocLite("polyester")
```

Polyester depends on the `Biostrings` and `IRanges` libraries from Bioconductor.

# Required Input
You will need either:

* a reference FASTA file containing names and sequences of transcripts from which reads should be simulated. Known transcripts from human chromosome 22 (hg19 build) are available in the `data` subdirectory of this package. 
* or a file in [GTF format](http://www.ensembl.org/info/website/upload/gff.html) denoting transcript structures, along with one FASTA file of the DNA sequence for each chromosome in the GTF file. All the FASTA files should be in the same directory. DNA sequences for some organisms can be downloaded [here](http://ccb.jhu.edu/software/tophat/igenomes.shtml) (sequences are in the `<organism>/<source>/<build>/Sequence/Chromosomes` subdirectory, e.g., `Homo_sapiens/UCSC/hg19/Sequence/Chromosomes`).

# Simulating reads

Simulating an RNA-seq experiment with Polyester requires just one function call. You can choose either `simulate_experiment()` or `simulate_experiment_countmat()`. 

### examples

A FASTA file called `chr22.fa` is provided with `polyester`. This file contains sequences for 918 transcripts on chromosome 22, as annotated in hg19. For this example, we will only simulate from the first 20 of these transcripts, and we will set the first 2 transcripts to be overexpressed in group A and the next 2 transcripts to be overexpressed in group B, each at a fold change of 4. A small experiment like this will only take a few seconds to run, even with many reps. A larger experiment (say, with all 918 transcripts) will run also run in a reasonable amount of time (minutes, not hours), with the exact timing depending on number of reads generated and number of reps.

To simulate a two-group experiment with 10 biological replicates in each group where the first 3 transcripts are differentially expressed with a fold change of 4, you can use code like this:


```r
library(polyester)
library(Biostrings)

fasta_file = system.file('extdata', 'chr22.fa', package='polyester')
fasta = readDNAStringSet(fasta_file)
small_fasta = fasta[1:20]
writeXStringSet(small_fasta, 'chr22_small.fa')
fold_changes = c(4, 4, 1/4, 1/4, rep(1, 16))
outdir = 'simulated_reads'

# ~20x coverage ----> reads per transcript = length/readlength * 20
# "width" is operating on a DNAStringSet (from Biostrings)
readspertx = round(20 * width(small_fasta) / 100)
simulate_experiment('chr22_small.fa', reads_per_transcript=readspertx, 
    num_reps=10, fold_changes=fold_changes, outdir=outdir) 
```

The `simulate_experiment` function draws the number of reads to simulate from each transcript from a negative binomial model. See below for details. Depending on your use case, it may be important to account for transcript length when deciding on the baseline mean number of reads to simulate from that transcript (as we did above with `readspertx`).

For more flexibility, you can use the `simulate_experiment_countmat` function. For example, we may want to simulate timecourse data. To do this, we can explicitly specify the number of reads for each transcript (rows), at each timepoint (columns). We will again only simulate from 20 transcripts.


```r
# set up matrix:
num_timepoints = 12
countmat = matrix(readspertx, nrow=length(small_fasta), ncol=num_timepoints)

# add spikes in expression at certain timepoints to certain transcripts:
up_early = c(1,2) 
up_late = c(3,4)
countmat[up_early, 2] = 3*countmat[up_early, 2]
countmat[up_early, 3] = round(1.5*countmat[up_early, 3])
countmat[up_late, 10] = 6*countmat[up_late, 10]
countmat[up_late, 11] = round(1.2*countmat[up_late, 11])

# simulate reads:
simulate_experiment_countmat('chr22_small.fa', readmat=countmat, 
    outdir='timecourse_reads') 
```

In this scenario, we simulated 12 total timepoints. We also added differential expression: transcripts 1 and 2 are overexpressed (compared to baseline) at timepoints 2 and 3, with a fold change of 3 at timepoint 2 and a fold change of 1.5 at timepoint 3. Similarly, transcripts 3 and 4 are overexpressed at timepoints 10 and 11, with a fold change of 6 at timepoint 10 and a fold change of 1.2 at timepoint 11. 

### More on the negative binomial model
The `simulate_experiment` function draws the number of reads to simulate from each transcript from a negative binomial distribution. For this function, you need to specify:
* `num_reps`: Number of biological replicates per experimental group (default: 10; can specify different numbers of replicates in the groups)
* `fold_changes`: A fold change for each transcript. This fold change represents the multiplicative change in the mean number of reads generated from each transcript, between the two experimental groups.
* `reads_per_transcript`: The baseline mean number of reads for each transcript. 
    - Fold changes compare the mean number of reads in group 1 to group 2. So a fold change of 0.5 means group 2's baseline mean number of reads for this transcript is twice that of group 1.
    - Long transcripts usually produce more reads in RNA-seq experiments than short ones, so you may want to specify `reads_per_transcript` as a function of transcript length
    - Default is 300 (regardless of transcript length).
* `size`: controls the per-transcript mean/variance relationship. In the negative binomial distribution, the mean/variance relationship is: ```mean = mean + (mean^2) / size```. You can specify the size for each transcript. By default, size is defined as 1/3 of the transcript's mean, which (in our experience) creates a somewhat idealized, low-variance situation.  Decrease the value of `size` to introduce more variance into your simulations.

### More on the count-matrix model 
The `simulate_experiment_readmat` function takes a count matrix as an argunent. Each row of this matrix represents a transcript, and each column represents a sample in the experiment. Entry `i,j` of the matrix specifies how many reads should be sampled from transcript `i` for sample `j`, allowing you to precisely and flexibly define the (differential) transcript expression structure for the experiment.

### other simulation parameters that can be set:
For both `simulate_experiment` and `simulate_experiment_countmat`, you can change these parameters:
* `fraglen`: Mean fragment length (default 250)
* `fragsd`: Standard devation of fragment lengths (default 25)
* `readlen`: Read length (default 100)
* `error_rate`: Sequencing error rate: probability that the sequencer records the wrong nucleotide at any given base (default 0.005, uniform error model assumed)
* `paired`: Whether the reads should be paired-end (default TRUE)

[This review paper](http://genomebiology.com/2010/11/12/220) (Oshlack, Robinson, and Young, _Genome Biology_ 2010, open access) provides a good overview of the RNA sequencing process, and might be particularly useful for understanding where some of these simulation parameters come into play.

If you'd like to explore specific steps in the sequencing process (fragmentation, reverse-complementing, error-adding), the functions called within `simulate_experiment` are also available and individually documented in Polyester.

See `?simulate_experiment` and `?simulate_experiment_countmat` for details on how to change these parameters.

### Using real data to guide simulation

To create a count matrix that resembles a real dataset, use the `create_read_numbers` function. To run this example, you will need to install the Ballgown package from Bioconductor if you do not already have it:


```r
source("http://bioconductor.org/biocLite.R")
biocLite("ballgown")
```


```r
library(ballgown)
```


```r
library(ballgown)
data(bg)
countmat = fpkm_to_counts(bg, threshold=0.01, mean_rps=400000)
params = get_params(countmat)
Ntranscripts = 50
Nsamples = 10
custom_readmat = create_read_numbers(mu=params$mu, fit=params$fit, p0=params$p0, m=Ntranscripts, n=Nsamples, seed=103)
```

```
## Generating data from baseline model.
```

The Ballgown package here is optional: the mean/variance relationship for each transcript can be estimated from any matrix of counts using `get_params`. You can add differential expression to the output from `create_read_numbers` (here, `custom_readmat`) and pass the resulting matrix to `simulate_experiment_countmat`.

## Output
A call to `simulate_experiment` or `simulate_experiment_countmat` will write FASTA files to the directory specified by the `outdir` argument. Reads in the FASTA file will be labeled with the transcript from which they were simulated.

If `paired` is true, you'll get two FASTA files per biological replicate (left mates are designated by the suffix `_1.fasta`; right mates by `_2.fasta`). If single-end reads are generated (`paired=FALSE`) you'll get one FASTA file per replicate. 

Files will be named `sample_01` through `sample_N` where `N` is the total number of replicates. The first `num_reps` (or `num_reps[1]`) samples belong to the same group in the two-group experiment scenario. 

In `simulate_experiment`, by default, a table called `sim_info.txt` is written to `outdir`, which will contain transcript IDs, fold changes, and whether or not that transcript was set to be differentially expressed. This file could be useful for downstream analysis. If the transcript names in the FASTA file cause problems down the line (e.g., a dangling single quote from a `5'-end` label), you can specify your own transcript names with the `transcriptid` argument. You will need to keep track of this information separately if you use `simulate_experiment_countmat.`

# Bug reports
Report bugs as issues on our [GitHub repository](https://github.com/alyssafrazee/polyester). 

# Session Information


```r
sessionInfo()
```

```
## R version 3.1.0 (2014-04-10)
## Platform: x86_64-apple-darwin10.8.0 (64-bit)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] ballgown_0.99.5      Biostrings_2.33.14   XVector_0.5.8       
##  [4] IRanges_1.99.28      S4Vectors_0.2.4      BiocGenerics_0.11.5 
##  [7] polyester_0.99.1     knitr_1.6            devtools_1.5        
## [10] BiocInstaller_1.15.5
## 
## loaded via a namespace (and not attached):
##  [1] annotate_1.43.5          AnnotationDbi_1.27.10   
##  [3] BatchJobs_1.3            BBmisc_1.7              
##  [5] Biobase_2.25.0           BiocParallel_0.99.19    
##  [7] bitops_1.0-6             brew_1.0-6              
##  [9] checkmate_1.4            codetools_0.2-9         
## [11] DBI_0.3.0                digest_0.6.4            
## [13] evaluate_0.5.5           fail_1.2                
## [15] foreach_1.4.2            formatR_1.0             
## [17] genefilter_1.47.6        GenomeInfoDb_1.1.19     
## [19] GenomicAlignments_1.1.29 GenomicRanges_1.17.40   
## [21] grid_3.1.0               httr_0.5                
## [23] iterators_1.0.7          lattice_0.20-29         
## [25] limma_3.21.15            Matrix_1.1-4            
## [27] memoise_0.2.1            mgcv_1.8-3              
## [29] nlme_3.1-117             RColorBrewer_1.0-5      
## [31] RCurl_1.95-4.3           Rsamtools_1.17.33       
## [33] RSQLite_0.11.4           rtracklayer_1.25.16     
## [35] sendmailR_1.1-2          splines_3.1.0           
## [37] stringr_0.6.2            survival_2.37-7         
## [39] sva_3.11.4               tools_3.1.0             
## [41] whisker_0.3-2            XML_3.98-1.1            
## [43] xtable_1.7-3             zlibbioc_1.11.1
```


