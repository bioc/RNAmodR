#' @include RNAmodR.R
#' @include class-RNAmodR.R
NULL

#' @name RNAmodR-mod-class
#' 
#' @title RNAmodR mod class
#' 
#' @description 
#' Virtual class for modification detection.
#' 
#' @slot analysisType 
#' @slot modType
#' 
#' @import methods
#' @export
setClass("mod",
         contains = "VIRTUAL",
         slots = c(analysisType = "character",
                   modType = "character"),
         prototype = list(analysisType = "default")
)
#' @rdname RNAmodR-mod-class
#'
#' @param object a RNAmodR mod class object 
#' 
#' @export
setMethod(
  f = "show", 
  signature = signature(object = "mod"),
  definition = function(object) {
    NULL
  }
)

#' @rdname mod-accessors
#' @aliases getAnalysisType getModType
#'
#' @param object a RNAmodR mod class object 
#' 
#' @return character defining the plot type for this modification class
#' @export
#'
#' @examples
#' \donttest{
#' getAnalysisType(mod)
#' getModType(mod)
#' }
setMethod(
  f = "getAnalysisType", 
  signature = signature(object = "mod"),
  definition = function(object) {
    return(object@analysisType)
  }
)
#' @rdname mod-accessors
#'
#' @return character defining the modification type for this modification
#' class
#' @export
setMethod(
  f = "getModType", 
  signature = signature(object = "mod"),
  definition = function(object) {
    return(object@modType)
  }
)
