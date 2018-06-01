#' @include RNAmodR.R
#' @include class-RNAmodR-args.R
NULL

#' @rdname RNAmodR-args-class
#'
#' @param x a RNAmodRargs class object
#' @param identifier an optional identifier for subsetting the arguments
#'
#' @return a data.frame with the selected arguments
#' @export
setMethod(
  f = "getArgs",
  signature = signature(x = "RNAmodRargs"),
  definition = function(x,
                        identifier,
                        param){
    # input check
    if(!missing(identifier)) {
      assertive::assert_is_a_string(identifier)
      .checkValueValidity(identifier,
                          unique(x@args$Identifier))
    }
    if(!missing(param)) {
      assertive::assert_is_a_string(param)
      .checkValueValidity(param,
                          unique(x@args$Param))
    }
    # if ModPipeline not given return all args
    if(missing(identifier) && missing(param)) {
      return(x@args)
    }
    if(missing(param)) {
      return(x@args[x@args$Identifier == identifier,])
    }
    if(missing(identifier)) {
      return(x@args[x@args$Param == param,"Value"])
    }
    return(x@args[x@args$Identifier == identifier &
                    x@args$Param == param,"Value"])
  }
)

#' @rdname RNAmodR-args-class
#'
#' @param x a RNAmodRargs class object
#'
#' @return a list of RNAmodR quantifier classes
#' @export
setMethod(
  f = "loadQuantifier",
  signature = signature(x = "RNAmodRargs"),
  definition = function(x){
    quantifier <- unique(getArgs(x, param = "data_type"))
    
  }
)

#' @rdname RNAmodR-args-class
#'
#' @param x a RNAmodRargs class object
#'
#' @return a list of RNAmodR identifier class
#' @export
setMethod(
  f = "loadIdentifier",
  signature = signature(x = "RNAmodRargs"),
  definition = function(x){
    # get modifictions
    modifications <- unique(x@args$Identifier)
    #
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
)
