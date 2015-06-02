# function to make sure everything is compatible and assign sane defaults
# internal

.check_extras = function(extras, paired){
  
  if(!('library_type' %in% names(extras))){
    extras$library_type = 'unstranded'
  }else{
    extras$library_type = match.arg(extras$library_type,
                                    c('unstranded', 'firststrand', 'secondstrand'))
  }
  
  if(!('strand_error_rate' %in% names(extras))){
    extras$strand_error_rate = 0
  }
  
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
  if(!('fraglen' %in% names(extras))){
    extras$fraglen = 250
  }
  if(!('fragsd' %in% names(extras))){
    extras$fragsd = 25
  }## I don't love this--these arguments aren't needed unless distr is normal.
  # but we store them anyway. should code better?
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
  
  return(extras)
  
}