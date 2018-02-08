#' @include RNAmodR.R
#' @include class-RNAmodR.R
NULL

#' @name RNAmodR-analysis-class
#' @title RNAmodR analysis class
#' @description
#' The class is virtual and has to be extended from, eg. \code{analysis_default}
#' The analysis class is the main function for detecting modifications. It works
#' usually in a three step manner.
#' \itemize{
#'   \item{"1."}{\code{convertReadsToPositions}: This function aggregates the
#'   information of reads per bam file, eg. per replicate and stores them for
#'   subsequent analysis}
#'   \item{"2."}{\code{parseMod}: This function searches for a defined pattern. 
#'   For this \code{mod} class objects will be created on the fly based on the 
#'   names of the modification one wants to detect and which can understand the
#'   data aggregated by \code{convertReadsToPositions}.}
#'   \item{"3."}{\code{mergePositionsOfReplicates}: This function aggregates the
#'   base data further for storage. This involves usually generating a per
#'   position mean and sd or equivalent.}
#' }
#' 
#' @slot plotType 
#' @slot data 
#' @slot modifications 
#' @return a RNAmodR analysis class object
#' @import methods
#' @export
setClass("analysis",
         contains = "VIRTUAL",
         slots = c(dataLabel = "character",
                   dataFormat = "function",
                   plotType = "character",
                   data = "list",
                   modifications = "list"),
         prototype = list(data = list(),
                          modifications = list())
)
#' @rdname RNAmodR-analysis-class
#' @param object a RNAmodR analysis object. 
#' @export
setMethod(
  f = "show", 
  signature = signature(object = "analysis"),
  definition = function(object) {
    NULL
  }
)

#' @rdname analysis-accessors
#' @aliases getPlotType getPositions getModifications
#' @title Accessor for \code{analysis} class objects
#' 
#' @description
#' The accessor function to \code{analysis} class objects can be used to access
#' the data saved in slots of the object. See examples for available functions.
#' 
#' @param object an analysis object 
#' @return character defining the plot type for this analysis class
#' @export
#' @examples
#' \donttest{
#' getPlotType(analysis)
#' getPositions(analysis)
#' getModifications(analysis)
#' }
setMethod(
  f = "getPlotType", 
  signature = signature(object = "analysis"),
  definition = function(object) {
    return(object@plotType)
  }
)

.get_plot_types_for_modifications <- function(modifications){
  assertive::assert_all_are_non_missing_nor_empty_character(modifications)
  modifications <- unique(modifications)
  l <- lapply(.load_analysis_classes(.get_analysis_type(modifications)),
              getPlotType)
  if(length(l) != length(modifications))
    stop("Incompatibble modification and analysis classes.")
  names(l) <- modifications
  l
}

#' @rdname analysis-accessors
#' @return a list of data as a lists of list(replicate) of DataFrame(transcript)
#' @export
setMethod(
  f = "getPositions", 
  signature = signature(object = "analysis"),
  definition = function(object) {
    return(object@data)
  }
)

#' @rdname analysis-accessors
#' @return a list of data as a lists of list(replicate) of DataFrame(transcript)
#' @export
setMethod(
  f = "getModifications", 
  signature = signature(object = "analysis"),
  definition = function(object) {
    return(object@modifications)
  }
)

#' @rdname analysis-accessors
#' @return the description of the type of position data, which can be used eg. 
#' as a label for printing
#' @export
setMethod(
  f = "getDataLabel", 
  signature = signature(object = "analysis"),
  definition = function(object) {
    return(object@dataLabel)
  }
)

.get_data_label_for_modifications <- function(modifications){
  assertive::assert_all_are_non_missing_nor_empty_character(modifications)
  modifications <- unique(modifications)
  l <- lapply(.load_analysis_classes(.get_analysis_type(modifications)),
              getDataLabel)
  if(length(l) != length(modifications))
    stop("Incompatibble modification and analysis classes.")
  names(l) <- modifications
  l
}

#' @rdname analysis-accessors
#' @return a function for formating the label of the position data
#' @export
setMethod(
  f = "getDataFormat", 
  signature = signature(object = "analysis"),
  definition = function(object) {
    return(object@dataFormat)
  }
)

.get_data_format_for_modifications <- function(modifications){
  assertive::assert_all_are_non_missing_nor_empty_character(modifications)
  modifications <- unique(modifications)
  l <- lapply(.load_analysis_classes(.get_analysis_type(modifications)),
              getDataFormat)
  if(length(l) != length(modifications))
    stop("Incompatibble modification and analysis classes.")
  names(l) <- modifications
  l
}


# analysis class handling ------------------------------------------------------

# load classes for modification analysis
.load_analysis_classes <- function(analysis){
  analysisClasses <- vector(mode = "list", length = length(analysis))
  names(analysisClasses) <- analysis
  for(i in seq_along(analysis)){
    className <- paste0("analysis_",analysis[[i]])
    # try to create modification detection classes
    tryCatch(
      class <- new(className),
      error = function(e) stop("Class for '",
                               analysis[[i]],
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
                               modifications[[i]],
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
