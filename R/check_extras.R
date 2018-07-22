# function to make sure everything is compatible and assign sane defaults
# internal

.check_extras = function(extras, paired, total.n){

    if(!('distr' %in% names(extras))){
        extras$distr = 'normal'
    }else{
        extras$distr = match.arg(extras$distr, 
            c('normal', 'empirical', 'custom'))
        if(extras$distr == 'custom' & !('custdens' %in% names(extras))){
            stop(.makepretty('to use custom fragment distribution, provide
                "custdens", a logspline object representing the distribution.'))
        }
    }

    # I don't love this--fraglen and fragsd aren't needed unless distr is normal.
    # but we store them anyway. should code better?
    if (!('fraglen' %in% names(extras))) {
        extras$fraglen = rep(250, total.n)
    } else {
      if (length(extras$fraglen) == 1) {
        extras$fraglen = rep(extras$fraglen, total.n)
      } else {
        stopifnot(length(extras$fraglen) == total.n)
      }
    }
    if (!('fragsd' %in% names(extras))) {
        extras$fragsd = rep(25, total.n)
    } else {
      if (length(extras$fragsd) == 1) {
        extras$fragsd = rep(extras$fragsd, total.n)
      } else {
        stopifnot(length(extras$fragsd) == total.n)
      }
    }

    if(!('readlen' %in% names(extras))){
        extras$readlen = 100
    }

    if(!('bias' %in% names(extras))){
        extras$bias = 'none'
    }else{
        extras$bias = match.arg(extras$bias, c('none', 'rnaf', 'cdnaf'))
    }

    if(!('error_model' %in% names(extras))){
        extras$error_model = 'uniform'
    }
    .check_error_model(extras, paired)

    if(!('error_rate' %in% names(extras))){
        extras$error_rate = 0.005
    }
    if(extras$error_model == 'custom'){
        extras$path = paste0(extras$model_path, '/', extras$model_prefix)
    }#this should work beause we already checked stuff.

    if(!('bias' %in% names(extras))){
        extras$bias = 'none'
    }else{
        extras$bias = match.arg(extras$bias, c('none', 'rnaf', 'cdnaf'))
    }

    if(!('lib_sizes' %in% names(extras))){
        extras$lib_sizes = rep(1, total.n)
    }else{
        stopifnot(is.numeric(extras$lib_sizes))
        stopifnot(length(extras$lib_sizes) == total.n)
    }

    if (!('frag_GC_bias' %in% names(extras))) {
      extras$frag_GC_bias <- 'none'
    } else {
      stopifnot(is.matrix(extras$frag_GC_bias))
      stopifnot(nrow(extras$frag_GC_bias) == 101)
      stopifnot(ncol(extras$frag_GC_bias) == total.n)
      stopifnot(all(extras$frag_GC_bias >= 0 & extras$frag_GC_bias <= 1))
    }

    if (!('strand_specific' %in% names(extras))) {
      extras$strand_specific <- FALSE
    }

    if (!('shuffle' %in% names(extras))) {
       extras$shuffle <- FALSE
    }
    return(extras)

}
