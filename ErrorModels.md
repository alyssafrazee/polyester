## Protocol/code for building platform-specific error models

The empirical error models that ship with Polyester were created using the [GemSIM framework](http://www.biomedcentral.com/1471-2164/13/74)<sup>1</sup>. GemSIM ships with error models estimated from real datasets sequenced on three different platforms: Illumina Genome Analyzer IIx with Illumina Sequencing Kit v4 chemistry (referred to as 'Illumina v4'), Illumina Genome Analyzer IIx with TrueSeq SBS Kit v5-GA (referred to as 'Illumina v5'), and Roche/454 FLX Titanium (referred to as 'Roche 454'). 

From these error models, we estimated error probabilities for each nucleotide in a read based on the read position and on the reference nucleotide at that base. Separate error probabilities were estimated for each of 4 specific sequencing errors (i.e., wrongly sequencing the base as each of the three incorrect nucleotides, or 'N'). 

To create the error models that ship with Polyester, we used the code in `build_error_models.py`. This code depends on [numpy](http://www.numpy.org/). I ran this in Python 2.7.5 on OSX 10.9.

```
python build_error_models.py /path/to/GemSIM_v1.6/models /path/to/output
```

We then saved these models as R objects in `data/` (so they can be used directly with Polyester) using the following R code:

```r
modfolder = '/path/to/output'
platforms = c('ill100v4_mate1', 'ill100v4_mate2', 
    'ill100v4_single', 'ill100v5_mate1', 'ill100v5_mate2', 
    'ill100v5_single', 'r454ti_single')
for(platform in platforms){
    i = which(platforms == platform)
    model = read.table(paste(modfolder, platform, sep='/'), header=TRUE)
    names(model)[1] = 'refbase'
    model$refbase = substr(model$refbase, 4, 4)
    eval(parse(text=paste0('model', i, ' <- model')))
    save(list=paste0('model',i), file=paste0('data/', platform, '.rda'), 
        compress='xz')
}
```

### References
McElroy KE, Luciani F, Thomas T (2012). GemSIM: general, error-model based simulator of next-generation sequencing data. _BMC Genomics_ *13*(1), 74.