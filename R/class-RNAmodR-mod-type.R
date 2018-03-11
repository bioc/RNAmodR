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
         slots = c(dataType = "character",
                   modType = "character",
                   positionOffset = "numeric")
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
#' @aliases getDataType getModType
#'
#' @param object a RNAmodR mod class object 
#' 
#' @return character defining the modification type for this modification
#' class
#' @export
#'
#' @examples
#' \donttest{
#' getDataType(mod)
#' getModType(mod)
#' }
setMethod(
  f = "getModType", 
  signature = signature(object = "mod"),
  definition = function(object) {
    return(object@modType)
  }
)

.check_modification_type <- function(modifications){
  assertive::assert_all_are_non_missing_nor_empty_character(modifications)
  modifications <- unique(modifications)
  l <- lapply(.load_mod_classes(modifications),
              getModType)
  if(unlist(l) != modifications){
    stop("Modification types do not match.")
  }
  names(l) <- modifications
  l
}

#' @rdname mod-accessors
#'
#' @return character defining the plot type for this modification class
#' 
#' @export
setMethod(
  f = "getDataType", 
  signature = signature(object = "mod"),
  definition = function(object) {
    return(object@dataType)
  }
)

.get_analysis_type <- function(modifications){
  assertive::assert_all_are_non_missing_nor_empty_character(modifications)
  modifications <- unique(modifications)
  l <- lapply(.load_mod_classes(modifications),
              getDataType)
  if(length(l) != length(modifications)){
    stop("Modification types do not match.")
  }
  names(l) <- modifications
  l
}

#' @rdname mod-accessors
#'
#' @return a position offset for scoring a modified positions, e.g. +1 for
#' a 5'-read stop at tlast unmodified position
#' 
#' @export
setMethod(
  f = "getPositionOffset", 
  signature = signature(object = "mod"),
  definition = function(object) {
    return(object@positionOffset)
  }
)

.get_position_offset <- function(modifications){
  assertive::assert_all_are_non_missing_nor_empty_character(modifications)
  modifications <- unique(modifications)
  l <- lapply(.load_mod_classes(modifications),
              getPositionOffset)
  if(length(l) != length(modifications)){
    stop("Modification types do not match.")
  }
  names(l) <- modifications
  l
}