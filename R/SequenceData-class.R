#' @include RNAmodR.R
#' @include normalization.R
#' @include SequenceDataFrame-class.R
NULL

#' @name SequenceData-class
#' @aliases SequenceData
#' 
#' @title The SequenceData class
#' 
#' @description 
#' The \code{SequenceData} class is implemented to contain data on each position
#' along a transcripts and holds the corresponding annotation data and
#' nucleotide sequence of these transcripts. To access this data several
#' \code{\link[=SequenceData-functions]{functions}} are available. The
#' \code{SequenceData} class is a virtual class, from which specific class can
#' be extended. Currently the following classes are implemented:
#' 
#' \itemize{
#' \item{\code{\link[=CoverageSequenceData-class]{CoverageSequenceData}}} 
#' \item{\code{\link[=EndSequenceData-class]{End5SequenceData}}, 
#' \code{\link[=EndSequenceData-class]{End3SequenceData}}, 
#' \code{\link[=EndSequenceData-class]{EndSequenceData}}}
#' \item{\code{\link[=NormEndSequenceData-class]{NormEnd5SequenceData}}, 
#' \code{\link[=NormEndSequenceData-class]{NormEnd5SequenceData}}}
#' \item{\code{\link[=PileupSequenceData-class]{PileupSequenceData}}}
#' \item{\code{\link[=ProtectedEndSequenceData-class]{ProtectedEndSequenceData}}}
#' }
#' 
#' It is derived from the
#' \code{\link[IRanges:DataFrameList-class]{CompressedSplitDataFrameList}} class
#' with additional slots for annotation and sequence data. Some functionality is
#' not inherited and not available, e.g. \code{cbind}, \code{rbind} amd
#' \code{relist}.
#' 
#' @param dataType The prefix for construction the class name of the 
#' \code{SequenceData} subclass to be constructed.
#' @param bamfiles the input which can be of the following types
#' \itemize{
#' \item{\code{BamFileList}:} {a named \code{BamFileList}}
#' \item{\code{character}:} {a \code{character} vector, which must be coercible
#' to a named \code{BamFileList} referencing existing bam files. Valid names are
#' \code{control} and \code{treated} to define conditions and replicates}
#' }
#' @param annotation annotation data, which must match the information contained
#' in the BAM files.
#' @param sequences sequences matching the target sequences the reads were 
#' mapped onto. This must match the information contained in the BAM files.
#' @param seqinfo optional \code{\link[GenomeInfoDb:Seqinfo]{Seqinfo}} to 
#' subset the transcripts analyzed on a chromosome basis.
#' @param ... Optional arguments overwriting default values, which are
#' \itemize{
#' \item{minQuality} {class dependent}
#' }
#' 
#' @slot ranges a \code{\link[GenomicRanges:GRangesList-class]{GRangesList}} 
#' object each element describing a transcript including its element. The 
#' \code{GRangesList} is constructed from the 
#' \code{\link[GenomicFeatures:transcriptsBy]{exonsBy(x, by="tx")}}.
#' @slot sequences a \code{\link[Biostrings:XStringSet-class]{XStringSet}} of 
#' type \code{sequencesType}.
#' @slot sequencesType a \code{character} value for the class name of 
#' \code{sequences}. Either \code{RNAStringSet} or \code{ModRNAStringSet}.
#' @slot bamfiles the input bam files as 
#' \code{\link[Rsamtools:BamFile-class]{BamFileList}}
#' @slot condition conditions along the 
#' \code{\link[Rsamtools:BamFile-class]{BamFileList}}: Either \code{control}
#' or \code{treated}
#' @slot replicate replicate number along the \code{BamFileList} for each of the
#' condition types.
#' @slot seqinfo a \code{\link[GenomeInfoDb:Seqinfo-class]{Seqinfo}} describing
#' the avialable chromosomes. This is an intersection of 
#' @slot minQuality a \code{integer} value describing a threshold of the minimum
#' quality of reads to be used.
NULL

setClass("SequenceData",
         contains = c("VIRTUAL", "CompressedSplitDataFrameList"),
         slots = c(ranges = "GRangesList",
                   sequences = "XStringSet",
                   sequencesType = "character",
                   bamfiles = "BamFileList",
                   condition = "factor",
                   replicate = "factor",
                   seqinfo = "Seqinfo",
                   minQuality = "integer",
                   unlistType = "character"),
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
#' @rdname RNAmodR-internals
setMethod("parallelSlotNames", "SequenceData",
          function(x) c("ranges","sequences", callNextMethod())
)

#' @rdname SequenceData-functions
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
    cat("-  ",class(object@seqinfo), " object with ", 
        summary(object@seqinfo), ":\n", sep = "")
  }
)
# validity ---------------------------------------------------------------------

.valid.SequenceData_elements <- function(x){
  nrow <- sum(unlist(width(ranges(x))))
  if(nrow != nrow(x@unlistData)){
    return("row number of data does not match position covered by annotation.")
  }
  if(nrow != sum(width(x@sequences))){
    return("Length of sequences does not match position covered by annotation.")
  }
  if(is.null(rownames(x@unlistData))){
    return("rownames of data is not set.")
  } else {
    seqs <- .seqs_rl(ranges(x))
    if(any(any(seqs != IRanges::IntegerList(rownames(x))))){
      return(paste0("Out of range rownames of data. The rownames do not match ",
                    "the ranges covered by the annotation data."))
    }
  }
  NULL
}

.valid.SequenceData <- function(x){
  c(.valid.SequenceData_elements(x),
    IRanges:::.valid.SimpleSplitDataFrameList(x))
}
S4Vectors::setValidity2(Class = "SequenceData",.valid.SequenceData)

# replacing --------------------------------------------------------------------

#' @rdname RNAmodR-internals
setReplaceMethod("[", "SequenceData",
  function(x, i, j, ..., value) {
    if (length(list(...)) > 0L){
      stop("invalid replacement")
    }
    if(!missing(j)){
      stop("replacement of columns not supported")
    }
    if(!is(value,class(x))){
      stop("replacement 'value' must be of the same class than 'x'")
    }
    if (missing(i)){
      x <- value
    } else {
      if(length(i) != length(value)){
        warning("number of items to replace is not a multiple of replacement ",
                "length")
        value <- value[seq_along(i)]
      }
      x@ranges[i] <- ranges(value)
      names(x@ranges)[i] <- names(ranges(value)) # must be set explicitly
      x@sequences[i] <- sequences(value)
      names(x@sequences)[i] <- names(sequences(value)) # must be set explicitly
      # rownames needs to be savid since a replace removes them
      rownames <- rownames(x)
      rownames[i] <- rownames(value)
      tmp <- callNextMethod(x = as(x,"SplitDataFrameList"), i = i,
                            value = as(value,"SplitDataFrameList"))
      x@unlistData <- tmp@unlistData
      x@partitioning <- tmp@partitioning
      rownames(x) <- rownames
    }
    validObject(x)
    return(x)
  }
)

#' @rdname RNAmodR-internals
setMethod("setListElement", "SequenceData",
  function(x, i, value){
    if(!is(value,"SequenceDataFrame")){
      stop("invalid value. must be 'SequenceDataFrame'.")
    }
    i2 <- S4Vectors::normalizeDoubleBracketSubscript(i, x, allow.append = TRUE,
                                                     allow.nomatch = TRUE)
    if(any(colnames(value) != colnames(unlist(x, use.names=FALSE)))){
      stop("'value' does not have matching colnames.")
    }
    x@ranges[[i2]] <- ranges(value)
    x@sequences[[i2]] <- sequences(value)
    # rownames needs to be savid since a replace removes them
    rownames <- rownames(x)
    rownames[[i2]] <- rownames(value)
    tmp <- callNextMethod(x = as(x,"SplitDataFrameList"), i = i,
                          value = as(value,"DataFrame"))
    x@unlistData <- tmp@unlistData
    x@partitioning <- tmp@partitioning
    rownames(x) <- rownames
    validObject(x)
    x
  }
)

# looping ----------------------------------------------------------------------

#' @importFrom IRanges PartitioningByEnd
lapply_SequenceData <- function(X, FUN, ...){
  FUN <- match.fun(FUN)
  ans <- vector(mode = "list", length = length(X))
  unlisted_X <- unlist(X, use.names = FALSE)
  X_partitioning <- IRanges::PartitioningByEnd(X)
  X_elt_width <- width(X_partitioning)
  empty_idx <- which(X_elt_width == 0L)
  if (length(empty_idx) != 0L){
    ans[empty_idx] <- NULL
  }
  non_empty_idx <- which(X_elt_width != 0L)
  if (length(non_empty_idx) == 0L){
    return(ans)
  }
  ans[non_empty_idx] <- 
    lapply(non_empty_idx, function(i){ FUN(getListElement(X,i), ...) })
  ans
}

setMethod("lapply", "SequenceData",
  function(X, FUN, ...){
    ans <- lapply_SequenceData(X, FUN, ...)
    names(ans) <- names(X)
    ans
  }
)

#' @importClassesFrom IRanges PartitioningByEnd
#' @importFrom IRanges PartitioningByEnd
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
    idx_on_unlisted_x <- 
      IRanges::IRanges(end = extractROWS(end(x@partitioning), i),
                       width = ans_eltNROWS)
    ans_unlistData <- extractROWS(x@unlistData, idx_on_unlisted_x)
    ans_partitioning <- new2("PartitioningByEnd", end = ans_breakpoints,
                             NAMES = extractROWS(names(x), i), check = FALSE)
    ans_elementMetadata <- extractROWS(x@elementMetadata, i)
    ans_ranges <- extractROWS(x@ranges, i)
    ans_sequences <- extractROWS(x@sequences, i)
    initialize(x, ranges = ans_ranges, sequences = ans_sequences,
               replicate = x@replicate, condition = x@condition,
               bamfiles = x@bamfiles, seqinfo = x@seqinfo, 
               minQuality = x@minQuality, unlistData = ans_unlistData,
               partitioning = ans_partitioning, 
               elementMetadata = ans_elementMetadata)
  }
)

#' @rdname RNAmodR-internals
#' @importFrom IRanges PartitioningByEnd
setMethod("getListElement", "SequenceData",
  function(x, i, exact=TRUE){
    i2 <- normalizeDoubleBracketSubscript(i, x, exact = exact,
                                          allow.NA = TRUE,
                                          allow.nomatch = TRUE)
    if (is.na(i2)){
      return(NULL)
    }
    unlisted_x <- unlist(x, use.names = FALSE)
    x_partitioning <- IRanges::PartitioningByEnd(x)
    window_start <- start(x_partitioning)[i2]
    window_end <- end(x_partitioning)[i2]
    new(x@unlistType,
        S4Vectors:::Vector_window(unlisted_x,
                                  start = window_start,
                                  end = window_end),
        x@ranges[[i2]],
        x@sequences[[i2]],
        x@replicate,
        x@condition)
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

# unlisting --------------------------------------------------------------------

setMethod("unlist", "SequenceData",
          function(x, recursive=TRUE, use.names=TRUE){
            if (!isTRUEorFALSE(use.names)){
              stop("'use.names' must be TRUE or FALSE")
            }
            unlisted_x <- x@unlistData
            if (use.names){
              unlisted_x <- S4Vectors:::set_unlisted_names(unlisted_x, x)
            }
            new(x@unlistType,
                unlisted_x,
                unlist(ranges(x), use.names = use.names),
                unlist(sequences(x)),
                x@replicate,
                x@condition)
          }
)

# accessors --------------------------------------------------------------------

setMethod("rownames", "SequenceData",
          function(x, do.NULL = TRUE, prefix = "row"){
            relist(rownames(x@unlistData),x@partitioning)
          }
)

# object creation --------------------------------------------------------------

.get_sequence_data_args <- function(input){
  minQuality <- .norm_min_quality(input, NULL)
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
  args <- list(minQuality = minQuality,
               max_depth = max_depth,
               minLength = minLength,
               maxLength = maxLength,
               minCoverage = minCoverage)
  args
}

.norm_min_quality <- function(input,minQuality){
  if(!is.null(input[["minQuality"]])){
    minQuality <- input[["minQuality"]]
    if(!is.integer(minQuality) | minQuality <= 1L){
      if(!is.na(minQuality)){
        stop("'minQuality' must be integer with a value higher than 1L.",
             call. = FALSE)
      }
    }
  }
  minQuality
}

.get_replicate_number <- function(bamfiles, conditions){
  control_rep <- seq_along(bamfiles[conditions == "control"])
  treated_rep <- seq_along(bamfiles[conditions == "treated"])
  rep <- c(control_rep,treated_rep)
  rep <- rep[c(which(conditions == "control"),
               which(conditions == "treated"))]
  factor(rep)
}

setMethod(
  f = "initialize", 
  signature = signature(.Object = "SequenceData"),
  definition = function(.Object, bamfiles, seqinfo, args, ...){
    if(!(.Object@sequencesType %in% c("RNAStringSet","ModRNAStringSet"))){
      stop("'sequencesType' must be either 'RNAStringSet' or 'ModRNAStringSet'")
    }
    if(missing(args) || !is.list(args)){
      args <- list()
    }
    # quality
    .Object@minQuality <- .norm_min_quality(args,.Object@minQuality)
    if(is.null(.Object@minQuality)){
      stop("Minimum quality is not set for '",class(.Object),"'.",
           call. = FALSE)
    }
    # set slots
    .Object@bamfiles <- bamfiles
    .Object@condition <- factor(names(bamfiles))
    .Object@replicate <- .get_replicate_number(bamfiles, .Object@condition)
    .Object@seqinfo <- .norm_seqinfo(seqinfo)
    # additional sanity checks
    .Object <- callNextMethod(.Object,...)
    .Object
  }
)

# internal SequenceData constructor
.new_SequenceData <- function(dataType, bamfiles, annotation, sequences, seqinfo,
                              ...){
  if(is.null(dataType)){
    stop("Invalid data type.")
  }
  args <- .get_sequence_data_args(list(...))
  className <- sequenceDataClass(dataType)
  # check bam files
  bamfiles <- .norm_bamfiles(bamfiles, className)
  # get annotation and sequence data
  annotation <- .norm_annotation(annotation, className)
  sequences <- .norm_sequences(sequences, className)
  seqinfo <- .norm_seqnames(bamfiles, annotation, sequences, seqinfo, className)
  # create the class
  ans <- new(className, bamfiles, seqinfo, args)
  # load transcript data and sequence data as well as the ScanBamParam
  grl <- .load_annotation(annotation)
  grl <- .subset_by_seqinfo(grl, seqinfo)
  sequences <- .load_transcript_sequences(sequences, grl)
  param <- .assemble_scanBamParam(grl, ans@minQuality, ans@seqinfo)
  # run the specific data aggregation function
  data <- getData(ans, grl, sequences, param, args)
  # post process the data
  .postprocess_read_data(ans, data, grl, sequences)
}

setMethod("SequenceData",
          signature = c(annotation = "character", sequences = "character"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
setMethod("SequenceData",
          signature = c(annotation = "character", sequences = "BSgenome"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
setMethod("SequenceData",
          signature = c(annotation = "TxDb", sequences = "character"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
setMethod("SequenceData",
          signature = c(annotation = "TxDb", sequences = "BSgenome"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
setMethod("SequenceData",
          signature = c(annotation = "GRangesList", sequences = "character"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
setMethod("SequenceData",
          signature = c(annotation = "GRangesList", sequences = "BSgenome"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
setMethod("SequenceData",
          signature = c(annotation = "GFF3File", sequences = "BSgenome"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
setMethod("SequenceData",
          signature = c(annotation = "GFF3File", sequences = "character"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
setMethod("SequenceData",
          signature = c(annotation = "character", sequences = "FaFile"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
setMethod("SequenceData",
          signature = c(annotation = "GFF3File", sequences = "FaFile"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
setMethod("SequenceData",
          signature = c(annotation = "TxDb", sequences = "FaFile"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
setMethod("SequenceData",
          signature = c(annotation = "GRangesList", sequences = "FaFile"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })


################################################################################
# common utility functions -----------------------------------------------------

# load annotation as GRangesList. one element per transcript
.load_annotation <- function(annotation){
  if(is(annotation,"TxDb")){
    ranges <- GenomicFeatures::exonsBy(annotation, by = "tx")
    rm(annotation)
    gc(FALSE)
  } else if(is(annotation,"GRangesList")) {
    ranges <- annotation
  } else {
    stop("Annotation is not a 'TxDb' or a 'GRangesList'.")
  }
  ranges
}

#' @importFrom Biostrings xscat
# load the transcript sequence per transcript aka. one sequence per GRangesList
# element
.load_transcript_sequences <- function(sequences, grl){
  seq <- Biostrings::getSeq(sequences, unlist(grl))
  seq <- split(seq,grl@partitioning)
  seq <- Reduce(c,lapply(seq,Biostrings::xscat))
  names(seq) <- names(grl)
  as(seq,"RNAStringSet")
}

# remove any elements, which are not in the seqinfo
.subset_by_seqinfo <- function(grl, seqinfo){
  grl <- grl[GenomicRanges::seqnames(grl) %in% GenomeInfoDb::seqnames(seqinfo)]
  grl <- grl[width(grl@partitioning) != 0L]
  GenomeInfoDb::seqlevels(grl) <- GenomeInfoDb::seqlevels(seqinfo)
  grl
}

################################################################################

#' @rdname SequenceData-functions
#' @export
setMethod("getData",
          signature = c(x = "SequenceData", grl = "GRangesList",
                        sequences = "XStringSet", param = "ScanBamParam"),
          definition = function(x, grl, sequences, param, args){
            stop("This functions needs to be implemented by '",class(x),"'.",
                 call. = FALSE)
          }
)

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

#' @importFrom IRanges PartitioningByWidth PartitioningByEnd
#' @importClassesFrom IRanges PartitioningByWidth PartitioningByEnd
.postprocess_read_data <- function(x, data, grl, sequences){
  conditionsFmultiplier <- length(data)
  # work with the unlisted data and construct a CompressedSplitDataFrameList
  # from this
  data <- .norm_postprocess_read_data(data)
  conditionsFmultiplier <- ncol(data) / conditionsFmultiplier 
  # create partitioning object from ranges
  partitioning <- IRanges::PartitioningByWidth(sum(width(grl)))
  if(sum(width(partitioning)) != nrow(data)){
    stop("Something went wrong. Length of data and Ranges do not match.")
  }
  # order data so that is matched the PartitioningByWidth object
  data <- relist(data, partitioning)
  positions <- .seqs_rl(grl)
  rownames(data) <- IRanges::CharacterList(positions)
  # order sequences
  sequences <- sequences[match(names(grl),names(sequences))]
  # store data
  x@unlistData <- data@unlistData
  x@partitioning <- data@partitioning
  x@replicate <- rep(x@replicate, each = conditionsFmultiplier)
  x@condition <- rep(x@condition, each = conditionsFmultiplier)
  x@ranges <- grl
  x@sequences <- as(sequences,x@sequencesType)
  names(x) <- names(grl)
  if(any(names(x@ranges) != names(x@sequences)) || 
     any(names(x@ranges) != names(x))){
    stop("Something went wrong.")
  }
  validObject(x)
  message("OK")
  x
}

# subset to conditions
.subset_to_condition <- function(conditions, condition){
  if(condition != "both"){
    f <- conditions == condition
    if(all(f == FALSE)){
      stop("No data for condition '",condition,"' found.")
    }
  } else {
    f <- rep(TRUE,length(conditions))
  }
  f
}


# accessors --------------------------------------------------------------------

#' @rdname SequenceData-functions
#' @export
setMethod(f = "seqinfo", 
          signature = signature(x = "SequenceData"),
          definition = function(x){x@seqinfo})

#' @rdname SequenceData-functions
#' @export
setMethod(f = "sequences", 
          signature = signature(x = "SequenceData"),
          definition = function(x){x@sequences})

#' @rdname SequenceData-functions
#' @export
setMethod(f = "ranges", 
          signature = signature(x = "SequenceData"),
          definition = function(x){x@ranges})

#' @rdname SequenceData-functions
#' @export
setMethod(f = "bamfiles", 
          signature = signature(x = "SequenceData"),
          definition = function(x){x@bamfiles})


# dummy functions --------------------------------------------------------------
# this needs to be implemented by each subclass

#' @rdname aggregate
#' @export
setMethod(f = "aggregate", 
          signature = signature(x = "SequenceData"),
          definition = 
            function(x, force){
              stop("This functions needs to be implemented by '",class(x),"'.",
                   call. = FALSE)
            }
)
