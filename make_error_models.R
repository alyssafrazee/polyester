## create error models to be released wtih polyester
## from the GemSIM/GemErr estimation model

modfolder = '../polyester_paper/error_models'
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

## figures for paper:
## put one example in true manuscript,
## include the others in supplementary data.
library(RSkittleBrewer)
library(usefulstuff)
colrs = RSkittleBrewer('tropical')
getColor = function(nt){
    switch(nt, A='black', C=colrs[1], G=colrs[2], T='deeppink', N=colrs[4])
}
nts = c('A','C','G','T','N')

plot_nt = function(model, nt){
    d = subset(model, refbase==nt)
    wrongnts = nts[-which(nts==nt)]
    errCols = paste0('read', wrongnts)
    mnum = as.matrix(model[,2:6])
    ymax = max(mnum[mnum<0.8])
    plot(1:100, as.matrix(d[errCols[1]])[-1], col=makeTransparent(getColor(wrongnts[1])),
        pch=19, cex=0.5, xlab='Read Position', ylab='Error Probability', 
        type='o', ylim=c(0,ymax))
    for(i in 2:4){
        points(1:100, as.matrix(d[errCols[i]])[-1], pch=19, cex=0.5, type='o',
            col=makeTransparent(getColor(wrongnts[i])))
    }
    title(nt)
    legend('topleft', wrongnts, pch=19, col=sapply(wrongnts, getColor), cex=0.7)
}

plot_model = function(model, file){
    pdf(file)
    par(mfrow=c(2,2))
    for(nt in c('A','C','G','T')){
        plot_nt(model, nt)
    }
    dev.off()
}

plot_model(model1, 'illumina4_mate1.pdf')
plot_model(model2, 'illumina4_mate2.pdf')
plot_model(model3, 'illumina4_single.pdf')
plot_model(model4, 'illumina5_mate1.pdf')
plot_model(model5, 'illumina5_mate2.pdf')
plot_model(model6, 'illumina5_single.pdf')





