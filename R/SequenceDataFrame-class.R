#' @include RNAmodR.R
NULL

#' @name SequenceData-functions
#' @aliases show,SequenceDataFrame-method
#' 
#' @title SequenceData/SequenceDataList/SequenceDataFrame functions
#' 
#' @description 
#' The \code{SequenceData}, \code{SequenceDataList} and
#' \code{SequenceDataFrame} share functionality. Have a look at the elements
#' listed under Usage.
#' 
#' @param x,object,value a \code{SequenceData}, \code{SequenceDataList},
#' \code{SequenceDataFrame} object.
#' @param grl a \code{GRangesList} from \code{exonsBy(..., by = "tx")}
#' @param sequences a \code{XStringSet} of type \code{RNAStringSet} or 
#' \code{ModRNAStringSet}
#' @param param a \code{\link[Rsamtools:ScanBamParam-class]{ScanBamParam}} 
#' object
#' @param args a list of addition arguments
#' 
#' @return 
#' \itemize{
#' \item{\code{seqinfo}} {a \code{Seqinfo} object}
#' \item{\code{sequences}} {a \code{RNAStingSet} object}
#' \item{\code{ranges}} {a \code{GRangesList} object with each element per 
#' transcript}
#' \item{\code{bamfiles}} {a \code{BamFileList} object}
#' }
#' 
#' @examples 
#' data(e5sd,package="RNAmodR")
#' # general accessors
#' seqinfo(e5sd)
#' sequences(e5sd)
#' ranges(e5sd)
#' bamfiles(e5sd)
NULL

#' @name SequenceDataFrame-class
#' @aliases SequenceDataFrame
#' 
#' @title The SequenceDataFrame class
#' 
#' @description 
#' The \code{SequenceDataFrame} class contains data for positions along a 
#' transcripts. It is used to describe elements from a \code{SequenceData}
#' object.
#' 
#' It is derived from the \code{\link[S4Vectors:DataFrame-class]{DataFrame}} 
#' class.
#' 
#' @param df the data as a \code{DataFrame}.
#' @param ranges a \code{GRanges} object containing all annotation elements
#' for a transcript.
#' @param sequence \code{XString} object describing the nucleotide sequence of 
#' the transcript.
#' @param condition The condition of each column or set of columns. Either 
#' \code{control} or \code{treated}.
#' @param replicate The replicate of each column or set of columns for the 
#' individual conditions
#' @param x,i,j,...,drop arguments used for 
#' \code{\link[S4Vectors:DataFrame-class]{subsetting}}.
#' 
#' @return a \code{SequenceDataFrame} object
#' 
#' @examples 
#' data(e5sd,package="RNAmodR")
#' # A SequenceDataFrame can is usually constructed by subsetting from 
#' # a SequenceData object
#' sdf <- e5sd[[1]]
NULL

# SequenceDataFrame ------------------------------------------------------------

#' @rdname SequenceDataFrame-class
#' @export
setClass(Class = "SequenceDataFrame",
         contains = c("DataFrame"),
         slots = c(ranges = "GRanges",
                   sequence = "XString",
                   condition = "factor",
                   replicate = "factor"))

# constructor ------------------------------------------------------------------

.SequenceDataFrame <- function(df, ranges, sequence, replicate, condition){
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
  new2("SequenceDataFrame",
       ranges = ranges,
       sequence = sequence,
       condition = condition,
       replicate = replicate,
       rownames = df@rownames,
       nrows = df@nrows,
       listData = df@listData,
       elementMetadata = df@elementMetadata,
       metadata = df@metadata)
}

#' @rdname SequenceDataFrame-class
#' @export
SequenceDataFrame <- function(df, ranges, sequence, replicate, condition){
  .SequenceDataFrame(df, ranges, sequence, replicate, condition)
}

.valid_SequenceDataFrame <-  function(x){
  if(nrow(x) != sum(width(ranges(x)))){
    return("data length and ranges width do not match.")
  }
  if(nrow(x) != length(sequences(x))){
    return("data length and sequence length do not match.")
  }
  S4Vectors:::.valid.DataFrame(x)
  NULL
}

S4Vectors::setValidity2(Class = "SequenceDataFrame", .valid_SequenceDataFrame)

# show -------------------------------------------------------------------------

#' @rdname SequenceData-functions
setMethod("show", "SequenceDataFrame",
          function(object){
            callNextMethod(object)
            cat("\ncontaining a ")
            show(object@ranges)
            cat("\nand a ")
            show(object@sequence)
          })

# relisting --------------------------------------------------------------------

setMethod("relist", c(flesh = "SequenceDataFrame", skeleton = "ANY"),
          function(flesh, skeleton){
            stop("Relisting is not supported for 'SequenceDataFrame'")
          }
)

# accessors --------------------------------------------------------------------

#' @rdname SequenceData-functions
#' @export
setMethod(
  f = "sequences", 
  signature = signature(x = "SequenceDataFrame"),
  definition = function(x){x@sequence})
#' @rdname SequenceData-functions
#' @export
setMethod(
  f = "ranges", 
  signature = signature(x = "SequenceDataFrame"),
  definition = function(x){x@ranges})

#' @importFrom stats setNames
#' @rdname SequenceDataFrame-class
#' @export
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
                xstub <- stats::setNames(seq_along(x), names(x))
                j <- normalizeSingleBracketSubscript(j, xstub)
              }
              x <- initialize(x, df = as(x,"DataFrame")[, j, drop = FALSE],
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
            } else {
              return(x) # early exit if subset is column-only
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
