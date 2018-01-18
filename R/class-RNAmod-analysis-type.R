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
         slots = c(plotType = "character",
                   data = "list",
                   modifications = "list"),
         prototype = list(plotType = "default",
                          data = list(),
                          modifications = list())
)
setMethod(
  f = "show", 
  signature = signature(object = "analysis"),
  definition = function(object) {
    
  }
)

#' @rdname mod-accessors
#'
#' @return character defining the plot type for this analysis class
#' @export
#'
#' @examples
#' \donttest{
#' getPlotType(analysis)
#' }
setMethod(
  f = "getPlotType", 
  signature = signature(object = "analysis"),
  definition = function(object) {
    return(object@plotType)
  }
)

#' @rdname mod-accessors
#'
#' @return a list of data as a lists of list(replicate) of DataFrame(transcript)
#' @export
#'
#' @examples
#' \donttest{
#' getPositions(analysis)
#' }
setMethod(
  f = "getPositions", 
  signature = signature(object = "analysis"),
  definition = function(object) {
    return(object@data)
  }
)

#' @rdname mod-accessors
#'
#' @return a list of data as a lists of list(replicate) of DataFrame(transcript)
#' @export
#'
#' @examples
#' \donttest{
#' getModifications(analysis)
#' }
setMethod(
  f = "getModifications", 
  signature = signature(object = "analysis"),
  definition = function(object) {
    return(object@modifications)
  }
)

#' @rdname mod-accessors
#'
#' @return returns the total number of modifications detected independent of the
#' type of modification
#' @export
#'
#' @examples
#' \donttest{
#' getNumberOfModifications(analysis)
#' }
setMethod(
  f = "getNumberOfModifications", 
  signature = signature(object = "analysis"),
  definition = function(object) {
    sum(unlist(lapply(object@modifications,nrow)))
  }
)




# analysis class handling ------------------------------------------------------

# load classes for modification analysis
.load_analysis_classes <- function(analysis){
  analysisClasses <- vector(mode = "list", length = length(analysis))
  names(analysisClasses) <- analysis
  for(i in seq_along(analysis)){
    className <- paste0("mod_",analysis[[i]])
    # try to create modification detection classes
    tryCatch(
      class <- new(className),
      error = function(e) stop("Class for '",
                               analysis[[1]],
                               "' analysis does not exist (",className,").",
                               call. = FALSE)
    )
    # if( !existsMethod("convertReadsToPositions",signature(class(class),
    #                                                       "numeric",
    #                                                       "GRanges",
    #                                                       "DataFrame") ) )
    #   stop("Function convertReadsToPositions() not defined for ",class(class))
    # if( !existsMethod("parseMod",signature(class(class),
    #                                        "GRanges",
    #                                        "FaFile",
    #                                        "list") ) )
    #   stop("Function parseMod() not defined for ",class(class))
    # if( !existsMethod("mergePositionsOfReplicates",signature(class(class),
    #                                                          "GRanges",
    #                                                          "FaFile",
    #                                                          "list") ) )
    #   stop("Function mergePositionsOfReplicates() not defined for ",
    #        class(class))
    analysisClasses[[i]] <- class
  }
  return(analysisClasses)
}


# modification class handling --------------------------------------------------

# load classes for modification analysis
.load_mod_classes <- function(modifications){
  modClasses <- vector(mode = "list", length = length(modifications))
  for(i in seq_along(modifications)){
    className <- paste0("mod_",modifications[[i]])
    # try to create modification detection classes
    tryCatch(
      class <- new(className),
      error = function(e) stop("Class for detecting ",
                               modifications[[1]],
                               " does not exist (",className,").",
                               call. = FALSE)
    )
    # if( !existsMethod("convertReadsToPositions",signature(class(class),
    #                                                       "numeric",
    #                                                       "GRanges",
    #                                                       "DataFrame") ) )
    #   stop("Function convertReadsToPositions() not defined for ",class(class))
    # if( !existsMethod("parseMod",signature(class(class),
    #                                        "GRanges",
    #                                        "FaFile",
    #                                        "list") ) )
    #   stop("Function parseMod() not defined for ",class(class))
    # if( !existsMethod("mergePositionsOfReplicates",signature(class(class),
    #                                                          "GRanges",
    #                                                          "FaFile",
    #                                                          "list") ) )
    #   stop("Function mergePositionsOfReplicates() not defined for ",
    #        class(class))
    modClasses[[i]] <- class
  }
  names(modClasses) <- modifications
  return(modClasses)
}
