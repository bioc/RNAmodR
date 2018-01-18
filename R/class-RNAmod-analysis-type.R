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
setClass("analysis",
         contains = "VIRTUAL",
         slots = c(plotType = "character"),
         prototype = list(plotType = "default")
)

setMethod(
  f = "show", 
  signature = signature(object = "analysis"),
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
#' getPlotType(mod)
#' }
setMethod(
  f = "getPlotType", 
  signature = signature(object = "analysis"),
  definition = function(object) {
    return(object@plotType)
  }
)