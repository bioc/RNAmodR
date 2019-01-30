#' @include RNAmodR.R
#' @include normalization.R
#' @include SequenceDataFrame-class.R
NULL

#' @name SequenceData
#' @title SequenceData
#' @description 
#' title
#' 
#' 
#' @param ... Optional arguments overwriting default values, which are
#' \itemize{
#' \item{minQuality}{class dependent}
#' }
#' 
NULL

setClass("SequenceData",
         contains = c("VIRTUAL",
                      "CompressedSplitDataFrameList"),
         slots = c(ranges = "GRangesList",
                   sequences = "XStringSet",
                   sequencesType = "character",
                   bamfiles = "BamFileList",
                   seqinfo = "Seqinfo",
                   minQuality = "integer",
                   chromosomes = "character",
                   unlistType = "character",
                   replicate = "factor",
                   conditions = "factor"),
         prototype = list(ranges = GRangesList(),
                          sequencesType = "RNAStringSet",
                          sequences = RNAStringSet(),
                          unlistType = "SequenceDataFrame"))

# class names must be compatible with this class name generation function
sequenceDataClass <- function(dataType){
  ans <- paste0(dataType,"SequenceData")
  tmp <- try(getClass(ans))
  if(is(tmp,"try-error")){
    stop("Class '",ans,"' not found: ",tmp)
  }
  
  ans
}

# for concat
#' @name SequenceData
#' @export
setMethod("parallelSlotNames", "SequenceData",
          function(x) c("ranges","sequences", callNextMethod())
)

setMethod("show", "SequenceData",
          function(object){
            k <- length(object)
            data_nc <- ncol(object@unlistData)
            ranges_mcols <- mcols(object@ranges@unlistData, use.names = FALSE)
            ranges_nmc <- if (is.null(ranges_mcols)) 0L else ncol(ranges_mcols)
            cat(classNameForDisplay(object), " with ", k, " elements ",
                "containing ",sep = "")
            cat(data_nc, ifelse(data_nc == 1L, " data column", " data columns"),
                " and ",ranges_nmc, ifelse(ranges_nmc == 1L, " metadata column",
                                           " metadata columns"),
                "\n",sep = "")
            out_data <- NULL
            # data
            if (data_nc > 0) {
              data_col_names <- colnames(object@unlistData)
              data_col_types <- 
                lapply(object@unlistData, function(x) {
                  paste0("<", classNameForDisplay(x)[1],">")
                })
              out_data <- 
                matrix(unlist(data_col_types, use.names = FALSE), nrow = 1,
                dimnames = list("", data_col_names))
            }
            cat("- Data columns:\n")
            print(out_data, quote = FALSE, right = TRUE)
            cat("-")
            cat(paste0(" ",capture.output(show(object@seqinfo)),"\n"))
          }
)
# validity ---------------------------------------------------------------------

.valid.SequenceData <- function(x){
  NULL
}
S4Vectors::setValidity2(Class = "SequenceData",.valid.SequenceData)

################################################################################
# list methods -----------------------------------------------------------------

.extractElement <- function(x, i){
  unlisted_x <- unlist(x, use.names=FALSE)
  x_partitioning <- PartitioningByEnd(x)
  window_start <- start(x_partitioning)[i]
  window_end <- end(x_partitioning)[i]
  df <- S4Vectors:::Vector_window(unlisted_x,
                                  start = window_start,
                                  end = window_end)
  new(x@unlistType,
      df,
      x@ranges[[i]],
      x@sequences[[i]],
      x@replicate,
      x@conditions)
}

# subsetting
#' @name SequenceData
#' @export
setMethod("getListElement", "SequenceData",
          function(x, i, exact=TRUE){
            i2 <- normalizeDoubleBracketSubscript(i, x, exact=exact,
                                                  allow.NA=TRUE,
                                                  allow.nomatch=TRUE)
            if (is.na(i2)){
              return(NULL)
            }
            .extractElement(x, i)
          }
)


# replacing --------------------------------------------------------------------

.remove_list_element <- function(x, i){
  stopifnot(isSingleNumberOrNA(i))
  if (is.na(i) || i < 1L || i > length(x))
    return(x)  # no-op
  ## `[<-.data.frame` does some terrible mangling of the colnames
  ## if they contain duplicates so we can't use it here.
  if (is.data.frame(x)) {
    x[[i]] <- NULL
    return(x)
  }
  x[-i]
}

.append_list_element <- function(x, value, name = NULL){
  if (is.null(name) && !is.null(names(x)))
    name <- ""
  value <- .wrap_in_length_one_list_like_object(value, name, x)
  coerce2(c(x, value), x)
}

.replace_list_element <- function(x, i, value){
  value <- .wrap_in_length_one_list_like_object(value, names(x)[[i]], x)
  ## `[<-` propagates the metadata columns from 'value' to 'x' but here
  ## we don't want that.
  if (is(x, "Vector"))
    x_mcols <- mcols(x)
  x[i] <- value
  if (is(x, "Vector"))
    mcols(x) <- x_mcols
  x
}


#' @name SequenceData
#' @export
setMethod("setListElement", "SequenceData",
          function(x, i, value){
            browser()
            if(!is(value,"SequenceDataFrame")){
              stop("invalid value. must be 'SequenceDataFrame'.")
            }
            i2 <- normalizeDoubleBracketSubscript(i, x,
                                                  allow.append=TRUE,
                                                  allow.nomatch=TRUE)
            if (is.null(value))
              return(.remove_list_element(x, i2))
            if (is.na(i2) || i2 > length(x)) {
              name <- if (is.na(i2)) as.character(i) else NULL
              return(.append_list_element(x, value, name))
            }
            .replace_list_element(x, i2, value)
          }
)

# looping ----------------------------------------------------------------------

lapply_SequenceData <- function(X, FUN, ...){
  FUN <- match.fun(FUN)
  ans <- vector(mode = "list", length = length(X))
  unlisted_X <- unlist(X, use.names = FALSE)
  X_partitioning <- PartitioningByEnd(X)
  X_elt_width <- width(X_partitioning)
  empty_idx <- which(X_elt_width == 0L)
  if (length(empty_idx) != 0L) 
    ans[empty_idx] <- list(FUN(extractROWS(unlisted_X, integer(0)), ...))
  non_empty_idx <- which(X_elt_width != 0L)
  if (length(non_empty_idx) == 0L)
    return(ans)
  ans[non_empty_idx] <-
    lapply(non_empty_idx,
           function(i){
             df <- getListElement(X,i)
             FUN(df,
                 ...)
           })
  ans
}

setMethod("lapply", "SequenceData",
          function(X, FUN, ...)
          {
            ans <- lapply_SequenceData(X, FUN, ...)
            names(ans) <- names(X)
            ans
          }
)

setMethod("extractROWS", "SequenceData",
          function(x, i){
            i <- normalizeSingleBracketSubscript(i, x, as.NSBS = TRUE)
            ans_eltNROWS <- extractROWS(width(x@partitioning), i)
            ans_breakpoints <- suppressWarnings(cumsum(ans_eltNROWS))
            nbreakpoints <- length(ans_breakpoints)
            if (nbreakpoints != 0L && is.na(ans_breakpoints[[nbreakpoints]])){
              stop("Subsetting operation on ", class(x), " object 'x' ",
                   "produces a result that is too big to be ",
                   "represented as a CompressedList object. ",
                   "This is not implemented, yet.")
            }
            idx_on_unlisted_x <- IRanges(end = extractROWS(end(x@partitioning), i),
                                         width = ans_eltNROWS)
            ans_unlistData <- extractROWS(x@unlistData, idx_on_unlisted_x)
            ans_partitioning <- new2("PartitioningByEnd",
                                     end = ans_breakpoints,
                                     NAMES = extractROWS(names(x), i),
                                     check = FALSE)
            ans_elementMetadata <- extractROWS(x@elementMetadata, i)
            ans_ranges <- extractROWS(x@ranges, i)
            ans_sequences <- extractROWS(x@sequences, i)
            initialize(x, 
                       ranges = ans_ranges,
                       sequences = ans_sequences,
                       bamfiles = x@bamfiles,
                       unlistData = ans_unlistData,
                       partitioning = ans_partitioning,
                       elementMetadata = ans_elementMetadata)
          }
)

setMethod("getListElement", "SequenceData",
          function(x, i, exact=TRUE){
            i2 <- normalizeDoubleBracketSubscript(i, x, exact = exact,
                                                  allow.NA = TRUE,
                                                  allow.nomatch = TRUE)
            if (is.na(i2))
              return(NULL)
            unlisted_x <- unlist(x, use.names=FALSE)
            x_partitioning <- PartitioningByEnd(x)
            window_start <- start(x_partitioning)[i2]
            window_end <- end(x_partitioning)[i2]
            new(x@unlistType,
                S4Vectors:::Vector_window(unlisted_x,
                                          start = window_start,
                                          end = window_end),
                x@ranges[[i2]],
                x@sequences[[i2]],
                x@replicate,
                x@conditions)
          }
)


# Concatenation ----------------------------------------------------------------

setMethod("cbind", "SequenceData",
          function(...){
            arg1 <- list(...)[[1L]]
            stop("'rbind' is not supported for ",class(arg1),".")
          }
)
setMethod("rbind", "SequenceData",
          function(...){
            arg1 <- list(...)[[1L]]
            stop("'rbind' is not supported for ",class(arg1),".")
          }
)

# object creation --------------------------------------------------------------

.norm_min_quality <- function(input,x){
  minQuality <- x@minQuality
  if(!is.null(input[["minQuality"]])){
    minQuality <- input[["minQuality"]]
    if(!is.integer(minQuality) | minQuality < 1L){
      if(!is.na(minQuality)){
        stop("'minQuality' must be integer with a value higher than 0L.",
             call. = FALSE)
      }
    }
  }
  minQuality
}

setMethod(
  f = "initialize", 
  signature = signature(.Object = "SequenceData"),
  definition = function(.Object, bamfiles, seqinfo, ...){
    # quality
    .Object@minQuality <- .norm_min_quality(list(...),.Object)
    if(!is.integer(.Object@minQuality) | 
       is.null(.Object@minQuality) | 
       .Object@minQuality == 0L){
      stop("Minimum quality is not set for '",class(.Object),"'.",
           call. = FALSE)
    }
    # set slots
    .Object@bamfiles <- bamfiles
    .Object@replicate <- factor(seq_along(bamfiles))
    .Object@conditions <- factor(names(bamfiles))
    .Object@seqinfo <- .norm_seqinfo(seqinfo)
    # additional sanity checks
    .Object <- callNextMethod(.Object,...)
    .Object
  }
)

# internal SequenceData constructor
.new_SequenceData <- function(dataType, files, annotation, sequences, seqinfo,
                              ...){
  if(is.null(dataType)){
    stop("Invalid data type.")
  }
  args <- list(...)
  className <- sequenceDataClass(dataType)
  # check bam files
  files <- .norm_bamfiles(files, className)
  # get annotation and sequence data
  txdb <- .norm_annotation(annotation, className)
  sequences <- .norm_sequences(sequences, className)
  seqinfo <- .norm_seqnames(files, txdb, sequences, seqinfo, className)
  # create the class
  ans <- new(className, files, seqinfo, args)
  # load transcript data and sequence data as well as the ScanBamParam
  grl <- .load_annotation(txdb, ans@seqinfo)
  sequences <- .load_transcript_sequences(sequences, grl)
  param <- .assemble_scanBamParam(grl, ans@minQuality, ans@seqinfo)
  # run the specific data aggregation function
  data <- .get_Data(ans, grl, sequences, param, args)
  # post process the data
  .postprocess_read_data(ans, data, grl, sequences)
}

setMethod("SequenceData",
          signature = c(annotation = "character", sequences = "character"),
          function(dataType, files, annotation, sequences, seqinfo, args, ...){
            .new_SequenceData(dataType, files, annotation, sequences, seqinfo,
                              ...)
          })
setMethod("SequenceData",
          signature = c(annotation = "character", sequences = "BSgenome"),
          function(dataType, files, annotation, sequences, seqinfo, args, ...){
            .new_SequenceData(dataType, files, annotation, sequences, seqinfo,
                              ...)
          })
setMethod("SequenceData",
          signature = c(annotation = "TxDb", sequences = "character"),
          function(dataType, files, annotation, sequences, seqinfo, args, ...){
            .new_SequenceData(dataType, files, annotation, sequences, seqinfo,
                              ...)
          })
setMethod("SequenceData",
          signature = c(annotation = "TxDb", sequences = "BSgenome"),
          function(dataType, files, annotation, sequences, seqinfo, args, ...){
            .new_SequenceData(dataType, files, annotation, sequences, seqinfo,
                              ...)
          })
setMethod("SequenceData",
          signature = c(annotation = "GFF3File", sequences = "BSgenome"),
          function(dataType, files, annotation, sequences, seqinfo, args, ...){
            .new_SequenceData(dataType, files, annotation, sequences, seqinfo,
                              ...)
          })
setMethod("SequenceData",
          signature = c(annotation = "GFF3File", sequences = "character"),
          function(dataType, files, annotation, sequences, seqinfo, args){
            .new_SequenceData(dataType, files, annotation, sequences, seqinfo,
                              ...)
          })
setMethod("SequenceData",
          signature = c(annotation = "character", sequences = "FaFile"),
          function(dataType, files, annotation, sequences, seqinfo, args, ...){
            .new_SequenceData(dataType, files, annotation, sequences, seqinfo,
                              ...)
          })
setMethod("SequenceData",
          signature = c(annotation = "GFF3File", sequences = "FaFile"),
          function(dataType, files, annotation, sequences, seqinfo, args, ...){
            .new_SequenceData(dataType, files, annotation, sequences, seqinfo,
                              ...)
          })
setMethod("SequenceData",
          signature = c(annotation = "TxDb", sequences = "FaFile"),
          function(dataType, files, annotation, sequences, seqinfo, args, ...){
            .new_SequenceData(dataType, files, annotation, sequences, seqinfo,
                              ...)
          })


################################################################################
# common utility functions -----------------------------------------------------

# load annotation as GRangesList. one element per transcript
.load_annotation <- function(txdb, seqinfo){
  ranges <- GenomicFeatures::exonsBy(txdb, by = "tx")
  ranges <- .subset_by_seqinfo(ranges, seqinfo)
  ranges
}

# load the transcript sequence per transcript aka. one sequence per GRangesList
# element
.load_transcript_sequences <- function(sequences,
                                       grl){
  seq <- getSeq(sequences, unlist(grl))
  seq <- split(seq,grl@partitioning)
  seq <- Reduce(c,lapply(seq,Biostrings::xscat))
  names(seq) <- names(grl)
  as(seq,"RNAStringSet")
}

# remove any elements, which are not in the seqinfo
.subset_by_seqinfo <- function(grl,seqinfo){
  grl <- grl[seqnames(grl) %in% seqnames(seqinfo)]
  grl <- grl[width(grl@partitioning) != 0L]
  seqlevels(grl) <- seqlevels(seqinfo)
  grl
}

################################################################################

.get_mod_data_args <- function(...){
  input <- list(...)
  max_depth <- 10000L # the default is 250, which is to small
  minLength <- NA
  maxLength <- NA
  minCoverage <- NA
  # for pileup
  if(!is.null(input[["max_depth"]])){
    max_depth <- input[["max_depth"]]
    if(!is.integer(max_depth) | max_depth <= 10L){
      stop("'max_depth' must be integer with a value higher than 10L.",
           call. = FALSE)
    }
  }
  # for protected end data
  if(!is.null(input[["minLength"]])){
    minLength <- input[["minLength"]]
    if(!is.integer(minLength) | minLength < 1L){
      if(!is.na(minLength)){
        stop("'minLength' must be integer with a value higher than 0L or NA.",
             call. = FALSE)
      }
    }
  }
  if(!is.na(minLength) && !is.na(maxLength)){
    if(minLength > maxLength){
      stop("'minLength' must be smaller or equal to 'maxLength'.",
           call. = FALSE)
    }
  }
  # for protected end data
  if(!is.null(input[["maxLength"]])){
    maxLength <- input[["maxLength"]]
    if(!is.integer(maxLength) | maxLength <= 1L){
      if(!is.na(maxLength)){
        stop("'maxLength' must be integer with a value higher than 1L or NA.",
             call. = FALSE)
      }
    }
  }
  # for norm protected end data
  if(!is.null(input[["minCoverage"]])){
    minCoverage <- input[["minCoverage"]]
    if(!is.integer(minCoverage) | minCoverage <= 1L){
      stop("'minCoverage' must be integer with a value higher than 1L or NA.",
           call. = FALSE)
    }
  }
  #
  args <- list(max_depth = max_depth,
               minLength = minLength,
               maxLength = maxLength,
               minCoverage = minCoverage)
  args
}

.norm_postprocess_read_data <- function(data){
  if(is(data[[1L]],"IntegerList") || is(data[[1L]],"NumericList")){
    data <- lapply(data, unlist)
    if(length(unique(lengths(data))) != 1L){
      stop("Data is of unequal length and cannot be coerced to a DataFrame.",
           call. = FALSE)
    }
    data <- S4Vectors::DataFrame(data)
  } else if(is(data[[1L]],"DataFrameList")) {
    data <- lapply(data, unlist, use.names = FALSE)
    ncols <- unique(unlist(lapply(data, ncol)))
    if(length(ncols) != 1L){
      stop("Something went wrong. Width of data should be constant.")
    }
    if(length(unique(vapply(data,nrow,numeric(1)))) != 1L){
      stop("Data is of unequal length and cannot be merged into a DataFrame.",
           call. = FALSE)
    }
    data <- do.call(cbind,data)
  } else {
    stop("Something went wrong.")
  }
  data
}

.postprocess_read_data <- function(x,
                                   data,
                                   grl,
                                   sequences){
  conditionsFmultiplier <- length(data)
  # work with the unlisted data and construct a CompressedSplitDataFrameList
  # from this
  data <- .norm_postprocess_read_data(data)
  conditionsFmultiplier <- ncol(data) / conditionsFmultiplier 
  # create partitioning object from ranges
  partitioning <- PartitioningByWidth(sum(width(grl)))
  if(sum(width(partitioning)) != nrow(data)){
    stop("Something went wrong. Length of data and Ranges do not match.")
  }
  # order data so that is matched the PartitioningByWidth object
  data <- IRanges::SplitDataFrameList(data)
  data@partitioning <- as(partitioning,"PartitioningByEnd")
  positions <- .seqs_rl(grl)
  rownames(data) <- IRanges::CharacterList(positions)
  # order sequences
  sequences <- sequences[match(names(grl),names(sequences))]
  # store data
  x@unlistData <- data@unlistData
  x@replicate <- rep(x@replicate, each = conditionsFmultiplier)
  x@conditions <- rep(x@conditions, each = conditionsFmultiplier)
  x@partitioning <- data@partitioning
  x@ranges <- grl
  x@sequences <- as(sequences,x@sequencesType)
  names(x) <- names(grl)
  if(any(names(x@ranges) != names(x@sequences)) || 
     any(names(x@ranges) != names(x))){
    stop("Something went wrong.")
  }
  message("OK")
  x
}

# subset to conditions
.subset_to_condition <- function(data,
                                 conditions,
                                 condition){
  if(condition != "both"){
    data <- data[,conditions == condition, drop = FALSE]
    if(ncol(data) == 0L){
      stop("No data for condition '",condition,"' found.")
    }
  }
  data
}


# accessors --------------------------------------------------------------------

#' @name SequenceData
#' @export
setMethod(f = "seqinfo", 
          signature = signature(x = "SequenceData"),
          definition = function(x){x@seqinfo})

#' @name SequenceData
#' @export
setMethod(f = "sequences", 
          signature = signature(x = "SequenceData"),
          definition = function(x){x@sequences})

#' @name SequenceData
#' @export
setMethod(f = "ranges", 
          signature = signature(x = "SequenceData"),
          definition = function(x){x@ranges})

#' @name SequenceData
#' @export
setMethod(f = "bamfiles", 
          signature = signature(x = "SequenceData"),
          definition = function(x){x@bamfiles})


# dummy functions --------------------------------------------------------------
# this needs to be implemented by each subclass

#' @name SequenceData
#' @export
setMethod(f = "aggregate", 
          signature = signature(x = "SequenceData"),
          definition = 
            function(x){
              stop("This functions needs to be implemented by '",class(x),"'.",
                   call. = FALSE)
            }
)

# data visualization -----------------------------------------------------------
# this needs to be implemented by each subclass
setMethod(
  f = ".dataTracks",
  signature = signature(x = "SequenceData",
                        data = "missing",
                        seqdata = "GRanges",
                        sequence = "XString"),
  definition = function(x, seqdata, sequence,  args) {
    stop("This functions needs to be implemented by '",class(x),"'.",
         call. = FALSE)
  }
)
