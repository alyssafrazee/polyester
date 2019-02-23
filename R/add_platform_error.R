process_custom = function(f){
    model = read.table(f, header=TRUE)
    names(model)[1] = 'refbase'
    model$refbase = substr(model$refbase, 4, 4)
    return(model)
}

#' @title Simulate sequencing error using empirical error model
#'
#' @description Given a sequencing platform and a set of sequencing reads,
#'   add sequencing errors to the reads given a known error profile from 
#'   the platform. 
#' @param tFrags DNAStringSetList containing error-free sequencing reads. If
#'   simulating a paired-end experiment, mate-pairs should appear next to each
#'   other in tFrags.
#' @param platform Which sequencing platform should the error model be estimated
#'   from? Currently supports \code{'illumina4'}, \code{'illumina5'},
#'   \code{'roche454'}, and \code{'custom'}.
#' @param paired Does \code{tFrags} contain paired end reads, with mate pairs 
#'   next to each other? (TRUE if yes.)
#' @param path if platform is \code{'custom'}, provide the path to the error model. 
#'   After processing the error model with \code{build_error_models.py}, you
#'   will have either two files (ending in _mate1 and _mate2, if your model was
#'   for paired-end reads) or one file (ending in _single, if your model was
#'   for single-end reads). The \code{path} argument should be the path to 
#'   the error model \emph{up to but not including} _mate1/_mate2/_single.
#' @export
#' @seealso \code{\link{add_error} for uniform error}
#' @return \code{DNAStringSet} object that is the same as \code{tFrags} except
#'   but with sequencing error added.
#' @details This function adds sequencing error to a set of reads based on the 
#'   position in the read and the true nucleotide at that location. 
#'   Position-specific probabilities of making each possible sequencing error
#'   (reading a T when it should have been A, reading a G when it should have
#'   been T, etc.) were calculated for each of three platforms using the 
#'   empirical error models available with the GemSIM software (see references).
#'   Users can also estimate an error model from their own data using GemSIM
#'   and can use that error model with Polyester as described in the vignette.
#'   (You will need to run a Python script available at the Polyester 
#'   GitHub repository to process the error model).
#' @references 
#'   McElroy KE, Luciani F and Thomas T (2012): GemSIM: general, 
#'   error-model based simulator of next-generation sequencing data. BMC 
#'   Genomics 13(1), 74.
#' @examples
#'   library(Biostrings)
#'   # pretend the srPhiX174 DNAStringSet represents 35bp single-end 
#'   # sequencing reads:
#'   data(srPhiX174) 
#'   set.seed(718)
#'   data_with_errors = add_platform_error(srPhiX174, 'illumina4', paired=FALSE)
#'   
#'   # the 17th read in this set has an error at position 20:
#'   data_with_errors[17][[1]][20] # N
#'   srPhiX174[17][[1]][20] # T
#'   
#'   # 101 reads total have at least one sequencing error:
#'   sum(data_with_errors != srPhiX174)
#'  
add_platform_error = function(tFrags, platform, paired, path=NULL){
    platform = match.arg(platform, c('illumina4', 'illumina5', 'roche454', 
        'custom'))
    
    if(platform == 'roche454' & paired){
        stop(.makepretty('The Roche 454 error model is only available for
            single-end reads'))
    } ### also checked in simulate_experiment (to fail early/loudly)

    if(max(width(tFrags)) > 101){
        stop(.makepretty('built-in platform-specific error models only
            available for reads < 102bp right now'))
    }
    if(platform == 'custom' & is.null(path)){
        stop('must provide path to model files if using custom error model')
    }

    #model = NULL #model will be overwritten later
    nt = c('A', 'T', 'G', 'C', 'N')

    if(paired){
        if(length(tFrags) %% 2 != 0){
            stop(.makepretty('Odd number of rows in tFrags. Paired-end reads
                should be in order, in complete pairs
                (e.g. read1a, read1b, read2a, read2b,..., readNa,
                readNb)'))
        }
        if(platform == 'illumina4'){
            data("ill100v4_mate1")
            m1 = model1
            data("ill100v4_mate2")
            m2 = model2
        }else if(platform == 'illumina5'){
            data("ill100v5_mate1")
            m1 = model4
            data("ill100v5_mate2")
            m2 = model5           
        }else{
            # platform is 'custom'
            m1 = process_custom(paste0(path, '_mate1'))
            m2 = process_custom(paste0(path, '_mate2'))
        }
        m1reads = tFrags[seq(1, length(tFrags), by=2)]
        m2reads = tFrags[seq(2, length(tFrags), by=2)]
        for(mate in c('left', 'right')){
            if(mate=='left'){
                reads = m1reads
                m = m1
            }else{
                reads = m2reads
                m = m2
            }
            L = width(reads)
            ## reformat error model into 3D array
            ## only include the actual error columns
            ## remove refBase and pos columns
            ## accessing matrix elements is faster than with a data.frame
            m <- as.matrix(m[1:(5*max(L)), 2:6])
            ## this gives us max(L) 5x5 error matrix for each position
            dim(m) <- c(5, max(L), 5)
            ## rename rows in each position's matrix with the base
            dimnames(m) <- list(nt, NULL, NULL)
            locations = matrix(FALSE, nrow=length(reads), ncol=max(L))
            reads = padAndClip(reads, IRanges(1, max(L)), Lpadding.letter="N",
                Rpadding.letter="N") #pad to rectangular
            replacements = rep("", length(reads))
            error_names = rep("", length(reads))
            for(pos in 1:max(L)){
                p = as.character(subseq(reads, pos, pos))
                ## retrieve error matrix for position
                errMat = m[, pos, ]
                ## loop through bases and sample new ones
                errRes <- lapply(nt, function(base){
                  sample(nt, sum(p == base),
                         replace = TRUE,
                         prob = errMat[base, ])
                })
                names(errRes) <- nt
                ## remove entries for bases that aren't present
                errRes <- errRes[which(sapply(errRes, length) > 0)]
                ## replace the old bases with the new bases
                errBases <- p
                for (base in names(errRes)) {
                  errBases[p == base] <- errRes[[base]]
                }
                ei = which(errBases != p)
                locations[ei, pos] = TRUE
                replacements[ei] = paste0(replacements[ei], errBases[ei])
                mutations = paste0(p[ei], pos, errBases[ei])
                error_names[ei] = sprintf("%s,%s", error_names[ei], mutations)
            }
            newReads = replaceLetterAt(reads, locations, replacements)
            if(mate == 'left'){
                m1reads = subseq(newReads, 1, L)
                m1err_names = error_names
            }else{
                m2reads = subseq(newReads, 1, L)
                m2err_names = error_names
            }
        }
        out = c(m1reads, m2reads)
        error_names = c(m1err_names, m2err_names)
        npairs = length(m1reads)
        outInds = rep(1:npairs, each=2)
        outInds[seq(2, length(outInds), by=2)] = (1:npairs)+npairs
        ret = out[outInds] # puts pairs of reads next to each other again
        error_names = error_names[outInds]
        names(ret) = names(tFrags)
        errors = which(error_names != "")
        error_names = gsub("^,", "", error_names)
        names(ret)[errors] = sprintf("%s;errors:%s", names(ret)[errors],
                                     error_names[errors])
        return(ret)
    } else {
        if(platform == 'illumina4'){
            data("ill100v4_single")
            model = model3
        }else if(platform == 'illumina5'){
            data("ill100v5_single")
            model = model6
        }else if(platform == 'roche454'){
            data("r454ti_single")
            model = model7
        }else{
            model = process_custom(paste0(path, '_single'))
        }
        L = width(tFrags)
        ## reformat error model into 3D array
        ## accessing matrix elements is faster than with a data.frame
        model <- as.matrix(model[1:(5*max(L)), 2:6])
        ## this gives us max(L) 5x5 error matrix for each position
        dim(model) <- c(5, max(L), 5)
        ## rename the rows in each position's matrix with the base
        dimnames(model) <- list(nt, NULL, NULL)
        locations = matrix(FALSE, nrow=length(tFrags), ncol=max(L))
        reads = padAndClip(tFrags, IRanges(1, max(L)), Lpadding.letter="N",
            Rpadding.letter="N")
        #^ pad to rectangular
        replacements = rep("", length(reads))
        error_names = rep("", length(reads))
        for(pos in 1:max(L)){
            p = as.character(subseq(reads, pos, pos))
            ## retrieve error matrix for position
            errMat = model[, pos, ]
            ## loop through bases and sample new ones
            errRes <- lapply(nt, function(base){
              sample(nt, sum(p == base),
                     replace = TRUE,
                     prob = errMat[base, ])
            })
            names(errRes) <- nt
            ## remove entries for bases that aren't present
            errRes <- errRes[which(sapply(errRes, length) > 0)]
            ## replace the old bases with the new bases
            errBases <- p
            for (base in names(errRes)) {
              errBases[p == base] <- errRes[[base]]
            }
            ei = which(errBases != p)
            locations[ei, pos] = TRUE
            replacements[ei] = paste0(replacements[ei], errBases[ei])
            mutations = paste0(p[ei], pos, errBases[ei])
            error_names[ei] = sprintf("%s,%s", error_names[ei], mutations)
        }
        newReads = replaceLetterAt(reads, locations, replacements)
        ret = subseq(newReads, 1, L)
        names(ret) = names(tFrags)
        errors = which(error_names != "")
        error_names = gsub("^,", "", error_names)
        names(ret)[errors] = sprintf("%s;errors:%s", names(ret)[errors],
                                     error_names[errors])
        return(ret)
    }
}
