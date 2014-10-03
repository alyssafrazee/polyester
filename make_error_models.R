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
