#' @include RNAmodR.R
NULL

#' @name PosDataFrame
#' 
#' @title PosDataFrame
#' 
#' @description 
#' title
#' 
#' 
NULL

# PosDataFrame -----------------------------------------------------------------

#' @rdname PosDataFrame
#' @export
setClass(Class = "PosDataFrame",
         contains = c("DataFrame"),
         slots = c(ranges = "GRanges",
                   sequence = "XString"))

setMethod(
  f = "initialize", 
  signature = signature(.Object = "PosDataFrame"),
  definition = function(.Object,
                        df,
                        ranges,
                        sequence){
    if(!is(df,"DataFrame")){
      stop("Invalid data object: ", class(df), " found, DataFrame expected.")
    }
    if(!is(ranges,"GRanges")){
      stop("Invalid data object: ", class(ranges), " found, GRanges expected.")
    }
    if(!is(sequence,"XString")){
      stop("Invalid data object: ", class(sequence), " found, XString expected.")
    }
    .Object@rownames <- df@rownames
    .Object@listData <- df@listData
    .Object@nrows <- df@nrows
    .Object@elementMetadata <- df@elementMetadata
    .Object@metadata <- df@metadata
    .Object@ranges <- ranges
    .Object@sequence <- sequence
    .Object
  }
)

#' @rdname PosDataFrame
#' @export
PosDataFrame <- function(df,ranges,sequence){
  new("PosDataFrame",
      df,
      ranges,
      sequence)
}

.valid_PosDataFrame <-  function(x){
  if(nrow(x) != width(ranges(x))){
    return("data length and ranges width do not match.")
  }
  if(nrow(x) != length(getSeq(x))){
    return("data length and sequence length do not match.")
  }
  S4Vectors::.valid.DataFrame(x)
  NULL
}

S4Vectors::setValidity2(Class = "PosDataFrame", .valid_PosDataFrame)

# show -------------------------------------------------------------------------
setMethod("show", "PosDataFrame",
          function(object){
            callNextMethod(object)
            cat("\ncontaining a ")
            show(object@ranges)
            cat("\nand a ")
            show(object@sequence)
          })

# accessors --------------------------------------------------------------------

#' @rdname PosDataFrame
#' @export
setMethod(
  f = "getSeq", 
  signature = signature(x = "PosDataFrame"),
  definition = function(x){x@sequence})
#' @rdname PosDataFrame
#' @export
setMethod(
  f = "ranges", 
  signature = signature(x = "PosDataFrame"),
  definition = function(x){x@ranges})

