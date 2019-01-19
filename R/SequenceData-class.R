#' @include RNAmodR.R
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
                   fasta = "FaFile",
                   gff = "GFFFile",
                   minQuality = "integer",
                   chromosomes = "character",
                   unlistType = "character",
                   replicate = "factor",
                   conditions = "factor"),
         prototype = list(ranges = GRangesList(),
                          sequencesType = "RNAStringSet",
                          sequences = RNAStringSet(),
                          unlistType = "SequenceDataFrame"))

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
            out_mdata <- NULL
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
            # metadata
            if (ranges_nmc > 0) {
              mdata_col_names <- colnames(ranges_mcols)
              mdata_col_types <- 
                lapply(ranges_mcols, function(x) {
                  paste0("<", classNameForDisplay(x)[1],">")
                })
              out_mdata <- matrix(unlist(mdata_col_types,
                                         use.names = FALSE),
                                  nrow = 1,
                                  dimnames = list("", mdata_col_names))
            }
            cat("- Data columns:\n")
            print(out_data, quote = FALSE, right = TRUE)
            cat("\n- Ranges metadata columns:\n")
            print(out_mdata, quote = FALSE, right = TRUE)
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
                       fasta = x@fasta,
                       gff = x@gff,
                       unlistData = ans_unlistData,
                       partitioning = ans_partitioning,
                       elementMetadata = ans_elementMetadata)
          }
)

setMethod("getListElement", "SequenceData",
          function(x, i, exact=TRUE){
            i2 <- normalizeDoubleBracketSubscript(i, x, exact=exact,
                                                  allow.NA=TRUE,
                                                  allow.nomatch=TRUE)
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
  definition = function(.Object,
                        bamfiles,
                        fasta,
                        gff,
                        args,
                        ...){
    className <- class(.Object)
    # quality
    .Object@minQuality <- .norm_min_quality(list(...),.Object)
    if(!is.integer(.Object@minQuality) | 
       is.null(.Object@minQuality) | 
       .Object@minQuality == 0L){
      stop("Minimum quality is not set for '",class(.Object),"'.",
           call. = FALSE)
    }
    # check bam files
    bamfiles <- .norm_bamfiles(bamfiles,className)
    # check genome sequences
    fasta <- .norm_fasta(fasta,className)
    # check genome annotation
    gff <- .norm_gff(gff,className)
    # set clots
    .Object@bamfiles <- bamfiles
    .Object@replicate <- factor(seq_along(bamfiles))
    .Object@conditions <- factor(names(bamfiles))
    .Object@fasta <- fasta
    .Object@gff <- gff
    # additional sanity checks
    sequenceNames <- unique(as.character(
      GenomeInfoDb::seqnames(scanFaIndex(path(.Object@fasta)))))
    chromosomeNames <- unique(as.character(
      GenomeInfoDb::seqnames(import(.Object@gff))))
    ### add bam seqnames check
    
    ###
    .Object@chromosomes <- intersect(sequenceNames,
                                     chromosomeNames)
    .Object@chromosomes <- .Object@chromosomes[order(.Object@chromosomes)]
    if(length(.Object@chromosomes) == 0){
      stop("chromosome information in annotation input and sequence names in ",
           "fasta file do not match.",
           call. = FALSE)
    }
    .Object <- callNextMethod(.Object,...)
    .Object
  }
)

################################################################################
# common utility functions -----------------------------------------------------

# load gff annotation prepare a GRangesList per parent element
# parent element has no parent themselves
.load_annotation <- function(gfffile){
  ranges <- import(gfffile)
  if(any(is.na(ranges$ID))){
    stop("ID column of annotation must not be empty.",
         call. = FALSE)
  }
  # keep only features which can contain mmodifications
  ranges <- .subset_mod_containing_features(ranges)
  # split GRanges per parent
  ranges <- split(ranges,.get_children_factor(ranges))
  # keep parent only data as metadata. this is basically the transcript to be
  # analyzed
  metadata(ranges)[["parents"]] <- .get_parent_annotations(ranges)
  if(any(metadata(ranges)[["parents"]]$ID != names(ranges))){
    stop("Something went wrong.", call. = FALSE)
  }
  ranges
}

.load_transcript_sequences <- function(fafile,
                                       ranges){
  # apperently the FaFile object does not like to be transferred. 
  # Therefore it is recreated on-the-fly.
  fafile <- path(fafile)
  # get transcript ranges
  gr <- .get_parent_annotations(ranges)
  # get sequence per transcript
  seq <- getSeq(FaFile(fafile), gr)
  names(seq) <- gr$ID
  as(seq,"RNAStringSet")
}

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
    data <- lapply(data,
                   function(d){
                     d[order(names(d))]
                   })
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
  rownames(data) <- NULL
  data
}

.postprocess_read_data <- function(x,
                                   data,
                                   ranges,
                                   sequences){
  conditionsFmultiplier <- length(data)
  # work with the unlisted data and construct a CompressedSplitDataFrameList
  # from this
  data <- .norm_postprocess_read_data(data)
  conditionsFmultiplier <- ncol(data) / conditionsFmultiplier 
  # create partitioning object from ranges
  parentRanges <- .get_parent_annotations(ranges)
  partitioning <- PartitioningByEnd(parentRanges)
  f <- Rle(parentRanges$ID,width(parentRanges))
  if(sum(width(partitioning)) != nrow(data) ||
     length(f) != nrow(data)){
    stop("Something went wrong. Length of data and Ranges do not match.")
  }
  # order data so that is matched the PartitioningByEnd object
  data <- split(data,f)
  data <- data[match(names(data),parentRanges$ID)]
  data <- unlist(data, use.names = FALSE)
  # order sequences
  sequences <- sequences[match(names(sequences),parentRanges$ID)]
  x@chromosomes <- x@chromosomes[match(x@chromosomes,
                                       GenomeInfoDb::seqnames(parentRanges))]
  # store data
  x@unlistData <- data
  x@replicate <- rep(x@replicate, each = conditionsFmultiplier)
  x@conditions <- rep(x@conditions, each = conditionsFmultiplier)
  x@partitioning <- partitioning
  x@ranges <- ranges
  x@sequences <- as(sequences,x@sequencesType)
  names(x) <- as.character(parentRanges$ID)
  if(any(names(x@ranges) != names(x@sequences)) || 
     any(names(x@ranges) != names(x))){
    stop("Something went wrong.")
  }
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
setMethod(f = "gff", 
          signature = signature(x = "SequenceData"),
          definition = function(x){x@gff})

#' @name SequenceData
#' @export
setMethod(f = "fasta", 
          signature = signature(x = "SequenceData"),
          definition = function(x){x@fasta})

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
# these need to be implemented by each subclass

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
