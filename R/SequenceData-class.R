#' @include RNAmodR.R
#' @include normalization.R
#' @include SequenceDataFrame-class.R
#' @include settings.R
NULL

#' @name SequenceData-class
#' @aliases SequenceData
#' 
#' @title The SequenceData class
#' 
#' @description 
#' The \code{SequenceData} class is implemented to contain data on each position
#' along transcripts and holds the corresponding annotation data and
#' nucleotide sequence of these transcripts. To access this data several
#' \code{\link[=SequenceData-functions]{functions}} are available. The
#' \code{SequenceData} class is a virtual class, from which specific classes can
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
#' The annotation and sequence data can be accessed through the functions
#' \code{ranges} and \code{sequences}, respectively. Beaware, that the data is
#' always provided according to genomic positions with increasing
#' \code{rownames}, but the sequence is given as the actual sequence of the
#' transcript. Therefore, it is necessary to treat the minus strand accordingly.
#' 
#' The \code{SequenceData} class is derived from the
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
#' @param ... Optional arguments overwriting default values. Not all 
#' \code{SequenceData} classes use all arguments. The arguments are:
#' \itemize{
#' \item{\code{minLength}} {single integer value setting a threshold for minimum
#' read length. Shorther reads are discarded (default: \code{minLength = NA}).}
#' \item{\code{maxLength}} {single integer value setting a threshold for maximum
#' read length. Longer reads are discarded (default: \code{maxLength = NA}).}
#' \item{\code{minQuality}} {single integer value setting a threshold for maximum
#' read quality. Reads with a lower quality are discarded (default: 
#' \code{minQuality = 5L}, but this is class dependent).}
#' \item{\code{max_depth}} {maximum depth for pileup loading (default: 
#' \code{max_depth = 10000L}).}
#' }
#' 
#' @slot ranges a \code{\link[GenomicRanges:GRangesList-class]{GRangesList}} 
#' object each element describing a transcript including its element. The 
#' \code{GRangesList} is constructed from the 
#' \code{\link[GenomicFeatures:transcriptsBy]{exonsBy(x, by="tx")}} function.
#' If during construction a \code{GRangesList} is provided instead of a 
#' character value pointing to a gff3 file or a \code{TxDb} object, it must have
#' a comparable structure. 
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
         slots = c(sequencesType = "character",
                   bamfiles = "BamFileList",
                   seqinfo = "Seqinfo",
                   minQuality = "integer",
                   unlistData = "SequenceDataFrame",
                   unlistType = "character",
                   dataDescription = "character"),
         prototype = list(sequencesType = "RNAStringSet"))

setMethod(
  f = "initialize",
  signature = signature(.Object = "SequenceData"),
  definition = function(.Object, ...){
    if(!assertive::is_a_non_empty_string(.Object@dataDescription)){
      stop("'dataDescription' must be a single non empty character value.")
    }
    if(!(.Object@sequencesType %in% c("RNAStringSet","ModRNAStringSet"))){
      stop("'sequencesType' must be either 'RNAStringSet' or 'ModRNAStringSet'")
    }
    callNextMethod(.Object, ...)
  }
)

# class names must be compatible with this class name generation function
sequenceDataClass <- function(dataType){
  ans <- paste0(dataType,"SequenceData")
  tmp <- try(getClass(ans))
  if(is(tmp,"try-error")){
    stop("Class '",ans,"' not found: ",tmp)
  }
  ans
}

#' @rdname SequenceData-functions
setMethod("show", "SequenceData",
  function(object){
    k <- length(object)
    unlisted_object <- object@unlistData
    data_nc <- ncol(unlisted_object)
    unlisted_ranges <- unlist(ranges(object),use.names = FALSE)
    ranges_mcols <- mcols(unlisted_ranges, use.names = FALSE)
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
      data_col_names <- colnames(unlisted_object)
      data_col_types <- 
        lapply(unlisted_object, function(x) {
          paste0("<", classNameForDisplay(x)[1],">")
        })
      out_data <- 
        matrix(unlist(data_col_types, use.names = FALSE), nrow = 1,
        dimnames = list("", data_col_names))
    }
    cat("- Data columns:\n")
    print(out_data, quote = FALSE, right = TRUE)
    cat("-  ",class(seqinfo(object)), " object with ", 
        summary(seqinfo(object)), ":\n", sep = "")
  }
)
# validity ---------------------------------------------------------------------

.valid.SequenceData_elements <- function(x){
  unlisted_x <- unlist(x, use.names=FALSE)
  nrow <- sum(width(ranges(unlisted_x)))
  if(nrow != nrow(unlisted_x)){
    return("row number of data does not match position covered by annotation.")
  }
  if(nrow != sum(width(sequences(x)))){
    return("Length of sequences does not match position covered by annotation.")
  }
  if(is.null(rownames(unlisted_x))){
    return("rownames of data is not set.")
  } else {
    seqs <- .seqs_rl_strand(ranges(x))
    if(any(any(seqs != IRanges::IntegerList(rownames(x))))){
      return(paste0("Out of range rownames of data. The rownames do not match ",
                    "the ranges covered by the annotation data."))
    }
  }
  NULL
}

.valid.SequenceData <- function(x){
  c(.valid.SequenceData_elements(x),
    IRanges:::.valid.CompressedList(x))
}
S4Vectors::setValidity2(Class = "SequenceData", .valid.SequenceData)

# coercion ---------------------------------------------------------------------

.as_SplitDataFrameList <- function(from){
  relist(as(unlist(from, use.names = FALSE),"DataFrame"),
         IRanges::PartitioningByWidth(from))
}
setAs("SequenceData", "SplitDataFrameList", .as_SplitDataFrameList)


# internals --------------------------------------------------------------------

#' @importClassesFrom IRanges PartitioningByEnd
#' @importFrom IRanges PartitioningByEnd
setMethod("extractROWS", "SequenceData",
  function(x, i){
    i <- normalizeSingleBracketSubscript(i, x, as.NSBS = TRUE)
    ans_eltNROWS <- extractROWS(width(IRanges::PartitioningByEnd(x)), i)
    ans_breakpoints <- suppressWarnings(cumsum(ans_eltNROWS))
    nbreakpoints <- length(ans_breakpoints)
    if (nbreakpoints != 0L && is.na(ans_breakpoints[[nbreakpoints]])){
      stop("Subsetting operation on ", class(x), " object 'x' ",
           "produces a result that is too big to be ",
           "represented as a CompressedList object. ",
           "This is not implemented, yet.")
    }
    idx_on_unlisted_x <- 
      IRanges::IRanges(end = extractROWS(end(IRanges::PartitioningByEnd(x)), i),
                       width = ans_eltNROWS)
    ans_unlistData <- extractROWS(unlist(x,use.names = FALSE),
                                  idx_on_unlisted_x)
    ans_partitioning <- new("PartitioningByEnd", end = ans_breakpoints,
                            NAMES = extractROWS(names(x), i))
    ans_elementMetadata <- extractROWS(x@elementMetadata, i)
    initialize(x, bamfiles = x@bamfiles, seqinfo = x@seqinfo, 
               minQuality = x@minQuality, unlistData = ans_unlistData,
               partitioning = ans_partitioning, 
               elementMetadata = ans_elementMetadata)
  }
)

setMethod("rownames", "SequenceData",
          function (x){
            ans <- rownames(unlist(x,use.names = FALSE), do.NULL = TRUE)
            relist(ans,x)
          }
)

# methods inherited from List and CompressedList, contain a coercion step
# x <- as(x, "List", strict = FALSE)
# 
# This does not keep the SequenceData object intact resulting in coercion
# to a CompressedSplitDataFrameList.
setMethod("[[", "SequenceData",
          function(x, i, j, ...) 
          {
            METHOD <- selectMethod("[[", "List")
            METHOD(x, i, j, ...)
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

# constructor ------------------------------------------------------------------

.quality_settings <- data.frame(
  variable = c("minQuality"),
  testFUN = c(".not_integer_bigger_equal_than_one"),
  errorValue = c(TRUE),
  errorMessage = c("'minQuality' must be integer with a value higher than 1L."),
  stringsAsFactors = FALSE)

.norm_min_quality <- function(input, minQuality){
  .norm_settings(input, .quality_settings, minQuality)[["minQuality"]]
}

.get_replicate_number <- function(bamfiles, conditions){
  control_rep <- seq_along(bamfiles[conditions == "control"])
  treated_rep <- seq_along(bamfiles[conditions == "treated"])
  rep <- c(control_rep,treated_rep)
  rep <- rep[c(which(conditions == "control"),
               which(conditions == "treated"))]
  factor(rep)
}

.check_positions_in_data <- function(data, positions){
  data_rownames <- lapply(data,rownames)
  data_rownames_non_empty <- !vapply(data_rownames,is.null,logical(1))
  # only check if data is present
  if(any(data_rownames_non_empty)){
    check <- lapply(data_rownames[data_rownames_non_empty],
                    "==", positions)
    check <- unlist(lapply(check,all))
    if(!all(check)){
      stop("rownames()/names() in data does not have the correct order. ",
           "rownames()/names() must be coercible to numeric and strictly ",
           "increasing.")
    }
  }
  NULL
}

# construct a result DataFrame by creating from Integer or NumericList
# or mergeing DataFrames into one
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
    stop("")
  }
  data
}

# reverse order of data, if originating from minus strand
.order_read_data_by_strand <- function(data, grl){
  # check if minus strand data is ordered in reverse
  strand_u <- .get_strand_u_GRangesList(grl)
  strand_minus <- strand_u == "-"
  if(any(strand_minus)){
    strand_minus_and_needs_rev <- strand_minus & 
      vapply(IRanges::IntegerList(rownames(data)),
             function(i){
               i[1] < i[2]
             },
             logical(1))
    if(any(strand_minus_and_needs_rev)){
      data[strand_minus_and_needs_rev] <- 
        lapply(data[strand_minus_and_needs_rev], rev)
    }
  }
  data
}

#' @importFrom IRanges PartitioningByWidth PartitioningByEnd
#' @importClassesFrom IRanges PartitioningByWidth PartitioningByEnd
.SequenceData <- function(className, bamfiles, ranges, sequences, seqinfo, args,
                          ...){
  ##############################################################################
  # setup additional variables
  ##############################################################################
  if(!is.list(args)){
    args <- as.list(args)
  }
  proto <- new(className)
  minQuality <- .norm_min_quality(args, proto@minQuality)
  condition <- factor(names(bamfiles))
  replicate <- .get_replicate_number(bamfiles, condition)
  if(!assertive::is_a_non_empty_string(proto@dataDescription)){
    stop("'dataDescription' must be a single non empty character value.")
  }
  if(!(proto@sequencesType %in% c("RNAStringSet","ModRNAStringSet"))){
    stop("'sequencesType' must be either 'RNAStringSet' or 'ModRNAStringSet'")
  }
  if(is.null(minQuality)){
    stop("Minimum quality is not set for '", className ,"'.",
         call. = FALSE)
  }
  param <- .assemble_scanBamParam(ranges, proto@minQuality, seqinfo)
  positions <- .seqs_rl_strand(ranges, force_continous = TRUE)
  ##############################################################################
  # run the specific data aggregation function
  ##############################################################################
  message("Loading ", proto@dataDescription, " from BAM files ... ",
          appendLF = FALSE)
  data <- getData(proto, bamfiles, ranges, sequences, param, args)
  names(data) <- paste0(names(data),".",condition,".",replicate)
  # check positions
  .check_positions_in_data(data, positions)
  ##############################################################################
  # post process the data
  ##############################################################################
  conditionsFmultiplier <- length(data)
  # work with the unlisted data and construct a CompressedSplitDataFrameList
  # from this
  data <- .norm_postprocess_read_data(data)
  # readjust replicate and condition to the normalized dimensions of the data
  conditionsFmultiplier <- ncol(data) / conditionsFmultiplier 
  replicate <- rep(replicate, each = conditionsFmultiplier)
  condition <- rep(condition, each = conditionsFmultiplier)
  # create partitioning object from ranges
  partitioning <- IRanges::PartitioningByWidth(sum(width(ranges)))
  if(sum(width(partitioning)) != nrow(data)){
    stop("Something went wrong. Length of data and Ranges do not match.")
  }
  # order data so that is matched the PartitioningByWidth object
  data <- relist(data, partitioning)
  rownames(data) <- IRanges::CharacterList(positions)
  data <- .order_read_data_by_strand(data, ranges)
  # order sequences
  sequences <- as(sequences, proto@sequencesType)
  sequences <- sequences[match(names(ranges), names(sequences))]
  # basic checks
  names(data) <- names(ranges)
  if(any(names(ranges) != names(sequences)) || 
     any(names(ranges) != names(data))){
    stop("")
  }
  ##############################################################################
  # Create SequenceData object
  ##############################################################################
  ans <- new(className, 
             bamfiles = bamfiles,
             seqinfo = seqinfo,
             minQuality = minQuality,
             unlistData = .SequenceDataFrame(gsub("SequenceData","",className),
                                             unlist(data, use.names = FALSE),
                                             unlist(ranges, use.names = FALSE),
                                             unlist(sequences, use.names = FALSE),
                                             replicate,
                                             condition),
             partitioning = IRanges::PartitioningByEnd(data),
             ...)
  message("OK")
  ans
}

.SequenceData_settings <- data.frame(
  variable = c("max_depth",
               "minLength",
               "maxLength"),
  testFUN = c(".not_integer_bigger_than_10",
              ".not_integer_bigger_equal_than_zero_nor_na",
              ".not_integer_bigger_equal_than_one_nor_na"),
  errorValue = c(TRUE,
                 TRUE,
                 TRUE),
  errorMessage = c("'max_depth' must be integer with a value higher than 10L.",
                   "'minLength' must be integer with a value higher than 0L or NA.",
                   "'maxLength' must be integer with a value higher than 1L or NA."),
  stringsAsFactors = FALSE)

.get_SequenceData_args <- function(input){
  minQuality <- .norm_min_quality(input, NULL)
  max_depth <- 10000L # the default is 250, which is to small
  minLength <- NA_integer_
  maxLength <- NA_integer_
  args <- .norm_settings(input, .SequenceData_settings, max_depth, minLength,
                         maxLength)
  if(!is.na(args[["minLength"]]) && !is.na(args[["maxLength"]])){
    if(args[["minLength"]] > args[["maxLength"]]){
      stop("'minLength' must be smaller or equal to 'maxLength'.",
           call. = FALSE)
    }
  }
  #
  args <- c(list(minQuality = minQuality),
            args)
  args
}

# internal SequenceData constructor
.new_SequenceData <- function(dataType, bamfiles, annotation, sequences, seqinfo,
                              ...){
  if(is.null(dataType)){
    stop("Invalid data type.")
  }
  args <- .get_SequenceData_args(list(...))
  className <- sequenceDataClass(dataType)
  # check bam files
  bamfiles <- .norm_bamfiles(bamfiles, className)
  # get annotation and sequence data
  annotation <- .norm_annotation(annotation, className)
  sequences <- .norm_sequences(sequences, className)
  seqinfo_missing <- missing(seqinfo)
  seqinfo <- .norm_seqnames(bamfiles, annotation, sequences, seqinfo, className)
  # load transcript data and sequence data
  grl <- .load_annotation(annotation)
  grl <- .subset_by_seqinfo(grl, seqinfo)
  if(length(grl) == 0L){
    if(seqinfo_missing){
      stop("No overlap between bamfiles and annotation.")
    } else {
      stop("No overlap between bamfiles, annotation and seqinfo.")
    }
  }
  sequences <- .load_transcript_sequences(sequences, grl)
  # create the class
  .SequenceData(className, bamfiles, grl, sequences, seqinfo, args)
}

# constructor utility functions ------------------------------------------------
# also used at other places

# check for multiple seqnames per ranges
.norm_unique_seqnames <- function(ranges){
  seqnames_ranges_u <- unique(GenomeInfoDb::seqnames(ranges))
  f <- lengths(seqnames_ranges_u) != 1L
  if(any(f)){
    message("Found transcript annotation with non unique seqnames. Removing ",
            "them ...")
    ranges <- ranges[!f]
    GenomeInfoDb::seqlevels(ranges) <- GenomeInfoDb::seqlevelsInUse(ranges)
  }
  ranges
}

# load annotation as GRangesList. one element per transcript
.load_annotation <- function(annotation){
  if(is(annotation,"TxDb")){
    ranges <- GenomicFeatures::exonsBy(annotation, by = "tx")
    ranges <- .norm_unique_seqnames(ranges)
  } else if(is(annotation,"GRangesList")) {
    ranges <- annotation
    ranges <- .norm_unique_seqnames(ranges)
    # make sure, that the elements are reverse order if on minus strand
    unlisted_ranges <- unlist(ranges)
    strand <- strand(unlisted_ranges) == "+"
    unlisted_ranges <- 
      c(unlisted_ranges[strand][order(unlisted_ranges[strand])],
        unlisted_ranges[!strand][rev(order(unlisted_ranges[!strand]))])
    ranges <- split(unname(unlisted_ranges),
                    factor(names(unlisted_ranges),
                           levels = unique(names(ranges))))
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
  seq <- relist(unlist(seq),IRanges::PartitioningByWidth(sum(width(grl))))
  names(seq) <- names(grl)
  as(seq,"RNAStringSet")
}

# remove any elements, which are not in the seqinfo
.subset_by_seqinfo <- function(grl, seqinfo){
  grl <- grl[GenomicRanges::seqnames(grl) %in% GenomeInfoDb::seqnames(seqinfo)]
  grl <- grl[width(IRanges::PartitioningByWidth(grl)) != 0L]
  GenomeInfoDb::seqlevels(grl) <- GenomeInfoDb::seqlevelsInUse(grl)
  grl
}

################################################################################

setGeneric( 
  name = "SequenceData",
  signature = c("annotation","sequences"),
  def = function(dataType, bamfiles, annotation, sequences, seqinfo, ...)
    standardGeneric("SequenceData")
) 

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

# Constructor end
################################################################################

#' @rdname SequenceData-functions
#' @export
setMethod("getData",
          signature = c(x = "SequenceData", bamfiles = "BamFileList", 
                        grl = "GRangesList", sequences = "XStringSet", 
                        param = "ScanBamParam"),
          definition = function(x, bamfiles, grl, sequences, param, args){
            stop("This functions needs to be implemented by '",class(x),"'.",
                 call. = FALSE)
          }
)

# accessors --------------------------------------------------------------------

#' @rdname SequenceData-functions
#' @export
setMethod(f = "bamfiles", 
          signature = signature(x = "SequenceData"),
          definition = function(x){x@bamfiles})
#' @rdname SequenceData-functions
#' @export
setMethod(f = "conditions", 
          signature = signature(object = "SequenceData"),
          definition = function(object){conditions(unlist(object))})
#' @rdname SequenceData-functions
#' @export
setMethod(
  f = "ranges", 
  signature = signature(x = "SequenceData"),
  definition = 
    function(x){
      partitioning <- IRanges::PartitioningByEnd(x)
      unlisted_ranges <- ranges(unlist(x))
      ends <- cumsum(width(unlisted_ranges)) == cumsum(width(partitioning))
      partitioning_relist <- IRanges::PartitioningByEnd(which(ends))
      names(partitioning_relist) <- names(x)
      if(length(x) != length(partitioning_relist)){
        stop("ranges could not be relisted.")
      }
      relist(unlisted_ranges, partitioning_relist)
    })
#' @rdname SequenceData-functions
#' @export
setMethod(f = "replicates", 
          signature = signature(x = "SequenceData"),
          definition = function(x){replicates(unlist(x))})
#' @rdname SequenceData-functions
#' @export
setMethod(f = "seqinfo", 
          signature = signature(x = "SequenceData"),
          definition = function(x){x@seqinfo})
#' @rdname SequenceData-functions
#' @export
setMethod(f = "sequences", 
          signature = signature(x = "SequenceData"),
          definition = function(x){relist(sequences(unlist(x)),x)})

# dummy functions --------------------------------------------------------------
# this needs to be implemented by each subclass

.check_aggregate_seqdata <- function(data, x){
  seqs <- .seqs_rl_strand(ranges(x))
  if(any(any(seqs != IRanges::IntegerList(rownames(data))))){
    stop("rownames() of aggregate data is not named or not named along the ",
         "genomic coordinates. Contact the maintainer of the class used.",
         call. = FALSE)
  }
  data
}

#' @rdname aggregate
#' @export
setMethod(f = "aggregate", 
          signature = signature(x = "SequenceData"),
          definition = 
            function(x, condition = c()){
              .check_aggregate_seqdata(aggregateData(x, condition), x)
            }
)

#' @rdname aggregate
#' @export
setMethod(f = "aggregateData", 
          signature = signature(x = "SequenceData"),
          definition = function(x, condition){
            stop("This functions needs to be implemented by '",class(x),"'.",
                 call. = FALSE)
          })
