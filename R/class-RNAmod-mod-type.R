#' @include RNAmod.R
NULL



#' @name mod
#' 
#' @title mod
#'
#' @return
#' @export
#'
#' @examples
setClass("mod",
         contains = "VIRTUAL",
         slots = c(analysisType = "character",
                   modType = "character"),
         prototype = list(analysisType = "default")
)

setMethod(
  f = "show", 
  signature = signature(object = "mod"),
  definition = function(object) {
    
  }
)


#' @rdname mod-accessors
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
