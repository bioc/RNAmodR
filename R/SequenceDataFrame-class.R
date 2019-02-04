#' @include RNAmodR.R
NULL

#' @name SequenceDataFrame
#' 
#' @title SequenceDataFrame
#' 
#' @description 
#' The \code{SequenceDataFrame} class contains data for positions along a 
#' transcripts. It is used to describe elements from a \code{SequenceData}
#' object.
#' 
#' @param ranges a \code{GRanges} object containing all annotation elements
#' for a transcript.
#' @param sequence \code{XString} object describing the nucleotide sequence of 
#' the transcript.
#' @param conditions The condition of each column or set of columns. Either 
#' \code{control} or \code{treated}.
#' @param replicate The replicate of each column or set of columns for the 
#' individual conditions
#' 
NULL

# SequenceDataFrame ------------------------------------------------------------

#' @rdname SequenceDataFrame
#' @export
setClass(Class = "SequenceDataFrame",
         contains = c("DataFrame"),
         slots = c(ranges = "GRanges",
                   sequence = "XString",
                   conditions = "factor",
                   replicate = "factor"))

setMethod(
  f = "initialize", 
  signature = signature(.Object = "SequenceDataFrame"),
  definition = function(.Object,
                        df,
                        ranges,
                        sequence,
                        replicate,
                        condition){
    if(!is(df,"DataFrame")){
      stop("Invalid data object: ", class(df), " found, DataFrame expected.")
    }
    if(ncol(df) != length(replicate) ||
       ncol(df) != length(condition)){
      stop("Replicate and Conditions information must match the DataFrame ",
           "dimensions.")
    }
    if(!is(ranges,"GRanges")){
      stop("Invalid data object: ", class(ranges), " found, GRanges expected.")
    }
    if(!is(sequence,"XString")){
      stop("Invalid data object: ", class(sequence), " found, XString expected.")
    }
    .Object@replicate <- replicate
    .Object@condition <- condition
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

#' @rdname SequenceDataFrame
#' @export
SequenceDataFrame <- function(df,ranges,sequence,replicate,conditions){
  new("SequenceDataFrame",
      df,
      ranges,
      sequence,
      replicate,
      conditions)
}

.valid_SequenceDataFrame <-  function(x){
  if(nrow(x) != width(ranges(x))){
    return("data length and ranges width do not match.")
  }
  if(nrow(x) != length(getSeq(x))){
    return("data length and sequence length do not match.")
  }
  S4Vectors::.valid.DataFrame(x)
  NULL
}

S4Vectors::setValidity2(Class = "SequenceDataFrame", .valid_SequenceDataFrame)

# show -------------------------------------------------------------------------
setMethod("show", "SequenceDataFrame",
          function(object){
            callNextMethod(object)
            cat("\ncontaining a ")
            show(object@ranges)
            cat("\nand a ")
            show(object@sequence)
          })

# accessors --------------------------------------------------------------------

#' @rdname SequenceDataFrame
#' @export
setMethod(
  f = "sequences", 
  signature = signature(x = "SequenceDataFrame"),
  definition = function(x){x@sequence})
#' @rdname SequenceDataFrame
#' @export
setMethod(
  f = "ranges", 
  signature = signature(x = "SequenceDataFrame"),
  definition = function(x){x@ranges})

setMethod("[", "SequenceDataFrame",
          function(x, i, j, ..., drop = TRUE){
            if (!isTRUEorFALSE(drop)){
              stop("'drop' must be TRUE or FALSE")
            }
            if (length(list(...)) > 0L){
              warning("parameters in '...' not supported")
            }
            
            ## We do list-style subsetting when [ was called with no ','.
            ## NOTE: matrix-style subsetting by logical matrix not supported.
            list_style_subsetting <- (nargs() - !missing(drop)) < 3L
            if (list_style_subsetting || !missing(j)) {
              if (list_style_subsetting) {
                if (!missing(drop))
                  warning("'drop' argument ignored by list-style subsetting")
                if (missing(i))
                  return(x)
                j <- i
              }
              if (!is(j, "IntegerRanges")) {
                xstub <- setNames(seq_along(x), names(x))
                j <- normalizeSingleBracketSubscript(j, xstub)
              }
              x <- initialize(x,
                              df = as(x,"DataFrame")[, j, drop = FALSE],
                              ranges = x@ranges,
                              sequence = x@sequence,
                              replicate = x@replicate[j],
                              condition = x@condition[j])
              if (anyDuplicated(names(x))){
                names(x) <- make.unique(names(x))
              }
              if (list_style_subsetting){
                return(x)
              }
            }
            if (!missing(i)){
              x <- extractROWS(x, i)
            }
            if (missing(drop)){
              drop <- TRUE
            }  
            if (drop) {
              ## one row left
              if (nrow(x) == 1L){
                return(as(x, "list"))
              }
            }
            as(x,"DataFrame")
          }
)
