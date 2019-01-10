#' @include RNAmodR.R
NULL

.norm_min_quality <- function(input,x){
  minQuality <- x@minQuality
  if(!is.null(input[["minQuality"]])){
    minQuality <- input[["minQuality"]]
    if(!is.integer(minQuality) | minQuality < 1L){
      stop("'minQuality' must be integer with a value higher than 0L.",
           call. = FALSE)
    }
  }
  minQuality
}