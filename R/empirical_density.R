#' @name empirical_density
#' @title Estimated distribution of fragment lengths
#' @description Empirical fragment length distribution was estimated using 7
#'   randomly selected RNA-seq samples from the GEUVADIS dataset ('t Hoen et
#'   al 2013). One sample was selected from each of the 7 laboratories that 
#'   performed the sequencing. We used Picard's "CollectInsertSizeMetrics" tool
#'   (http://broadinstitude.github.io/picard/), version 1.121, to estimate
#'   the fragment size distribution based on read alignments. Code we used
#'   to estimate this distribution is available at
#'   \url{https://github.com/alyssafrazee/polyester/blob/master/make_fraglen_model.R}.
#' @docType data
#' @format \code{logspline} object (created with \code{\link{logspline}}) 
#'   specifying the empirical density of fragment lengths in the 7 GEUVADIS
#'   samples.
#' @references 
#'   't Hoen PA, et al (2013): Reproducibility of high-throughput mRNA and 
#'   small RNA sequencing across laboratories. Nature Biotechnology 31(11):
#'   1015-1022.
NULL 