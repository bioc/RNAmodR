#' @include RNAmodR.R
#' @include class-RNAmodR-args.R
NULL

#' @rdname RNAmodR-args-class
#'
#' @param x a RNAmodRargs class object
#'
#' @return a list of input files
#' @export
setMethod(
  f = "getInputFiles",
  signature = signature(x = "RNAmodRargs"),
  definition = function(x){
    return(x@files)
  }
)
#' @rdname RNAmodR-args-class
#'
#' @param x a RNAmodRargs class object
#'
#' @return the conditions of each input file. Either "Control" or "Treated"
#' @export
setMethod(
  f = "getConditions",
  signature = signature(x = "RNAmodRargs"),
  definition = function(x){
    return(x@conditions)
  }
)

#' @rdname RNAmodR-args-class
#'
#' @param x a RNAmodRargs class object
#' @param identifier an optional identifier for subsetting the arguments
#' @param param an optional param for subsetting the arguments
#' @param value a value for the chosen parameter of the given identifier
#'
#' @return a data.frame with the selected arguments
#' @export
setMethod(
  f = "getParam",
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
#' @param identifier an optional identifier for subsetting the arguments
#'
#' @return a data.frame with the selected arguments
#' @export
setMethod(
  f = "setParam",
  signature = signature(x = "RNAmodRargs",
                        identifier = "character",
                        param = "character"),
  definition = function(x,
                        identifier,
                        param,
                        value){
    # input check
    assertive::assert_is_a_string(identifier)
    assertive::assert_is_a_string(param)
    .checkValueValidity(identifier,
                        unique(x@args$Identifier))
    .checkValueValidity(param,
                        unique(x@args$Param))
    if(!missing(value)){
      stop("No value given.",
           call. = FALSE)
    }
    x@args[x@args$Identifier == identifier &
             x@args$Param == param,"Value"] <- value
    return(invisible(value))
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
    quantifier <- unique(getParam(x, param = "data_type"))
    quantifierClasses <- 
      lapply(quantifier,
             function(quant){
               className <- paste0("RNAmodRquant_",quant)
               # try to create modification detection classes
               tryCatch(
                 class <- new(className),
                 error = function(e) stop("Class for gathering '",
                                          quant,
                                          "' data does not exist (",
                                          className,
                                          ").",
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
               
             })
    names(quantifierClasses) <- quantifier
    return(quantifierClasses)
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
    identifier <- unique(x@args$Identifier)
    #
    identifierClasses <- 
      lapply(identifier,
             function(ident){
               className <- paste0("RNAmodRident_",ident)
               params <- getParam(x,
                                  identifier = ident)
               param <- params$Value
               names(param) <- params$Param
               # try to create modification detection classes
               tryCatch(
                 class <- new(className,
                              param),
                 error = function(e) stop("Class for detecting ",
                                          ident,
                                          " modification does not exist (",
                                          className,
                                          ").",
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
               
             })
    names(identifierClasses) <- identifier
    return(identifierClasses)
    
  }
)
