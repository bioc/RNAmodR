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
#' @md
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
#' not inherited and might not available to full extend, e.g.\code{relist}.
#' 
#' **SequenceDataFrame**
#' 
#' The \code{SequenceDataFrame} class is a virtual class and  contains data for
#' positions along a single transcript. In addition to being used for returning
#' elements from a \code{SequenceData} object, the SequenceDataFrame class is
#' used to store the unlisted data within a
#' \code{\link[=SequenceData-class]{SequenceData}} object. Therefore, a matching
#' \code{SequenceData} and \code{SequenceDataFrame} class must be implemented.
#' 
#' The \code{SequenceDataFrame} class is derived from the
#' \code{\link[S4Vectors:DataFrame-class]{DataFrame}} class.
#' 
#' Subsetting of a \code{SequenceDataFrame} returns a \code{SequenceDataFrame} or 
#' \code{DataFrame}, if it is subset by a column or row, respectively. The 
#' \code{drop} argument is ignored for column subsetting.
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
#' @param deparse.level See \code{\link[base:cbind]{base::cbind}} for a 
#' description of this argument.
#' 
#' @slot sequencesType a \code{character} value for the class name of 
#' \code{sequences}. Either \code{RNAStringSet}, \code{ModRNAStringSet}, 
#' \code{DNAStringSet} or \code{ModDNAStringSet}.
#' @slot minQuality a \code{integer} value describing a threshold of the minimum
#' quality of reads to be used.
#' 
#' @return A SequenceData object
NULL

#' @name RNAmodR-development
#' @export
setClass("SequenceData",
         contains = c("VIRTUAL", "CompressedSplitDataFrameList"),
         slots = c(minQuality = "integer",
                   unlistData = "SequenceDataFrame",
                   unlistType = "character",
                   dataDescription = "character"))

setMethod(
  f = "initialize",
  signature = signature(.Object = "SequenceData"),
  definition = function(.Object, ...){
    if(!.is_non_empty_string(.Object@dataDescription)){
      stop("'dataDescription' must be a single non empty character value.")
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

setMethod("classNameForDisplay", "SequenceData",
          function(x) class(x)
)

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

coerceSequenceDataToCompressedSplitDataFrameList <- function(className){
  setMethod("coerce",
    signature = c(from = className, to = "CompressedSplitDataFrameList"),
    function(from, to, strict = TRUE){
      if (strict) {
        from <- from
        {
          value <- new("CompressedSplitDataFrameList")
          for (what in c("elementType", "elementMetadata", 
                         "metadata", "unlistData", "partitioning"
          )) slot(value, what) <- slot(from, what)
          value@unlistData <- as(value@unlistData,"DataFrame")
          value
        }
      } else from
    })
}

coerceToSequenceData <- function(className) {
  function(from) {
    if(is.list(from)) {
      classes <- unlist(lapply(from,class))
      from <- from[classes == paste0(className,"Frame")]
      if(length(from) == 0) {
        FUN <- match.fun(className)
        from <- list(FUN())
      }
    } else {
      if(is(from,className)){
        return(from)
      } else if(is(from,paste0(className,"Frame"))) {
        from <- list(from)
      } else {
        stop("Cannot coerce ",class(from)," to ",className,".")
      }
    }
    IRanges:::coerceToCompressedList(from)
  }
}

setSequenceDataCoercions <- function(type) {
  className <- sequenceDataClass(type)
  coerceSequenceDataToCompressedSplitDataFrameList(className)
  setAs("ANY", className, coerceToSequenceData(className))
  setAs("list", className, coerceToSequenceData(className))
}

# internals --------------------------------------------------------------------



# Accessors --------------------------------------------------------------------

setMethod("rownames", "SequenceData",
  function (x){
    ans <- rownames(unlist(x,use.names = FALSE), do.NULL = TRUE)
    relist(ans,x)
  }
)

# Concatenation ----------------------------------------------------------------

.check_ranges <- function(args){
  ranges <- lapply(args,ranges)
  ranges <- vapply(ranges[seq.int(2L,length(ranges))],
                   function(r){
                     all(all(r == ranges[[1L]]))
                   },
                   logical(1))
  if(!all(ranges)){
    stop("Inputs must have the same ranges.")
  }
}

.check_sequences <- function(args){
  sequences <- lapply(args,sequences)
  sequences <- vapply(sequences[seq.int(2L,length(sequences))],
                      function(s){
                        all(s == sequences[[1L]])
                      },
                      logical(1))
  if(!all(sequences)){
    stop("Inputs must have the same sequences.")
  }
}

.check_bamfiles <- function(args){
  bamfiles <- lapply(args,bamfiles)
  bamfiles <- vapply(bamfiles[seq.int(2L,length(bamfiles))],
                     function(b){
                       all(path(b) == path(bamfiles[[1L]]))
                     },
                     logical(1))
  if(!all(bamfiles)){
    stop("Inputs must be derived from the same bamfiles.")
  }
}

#' @rdname SequenceData-class
#' @export
setMethod("cbind", "SequenceData",
          function(..., deparse.level = 1) 
          {
            args <- list(...)
            if(length(args) == 1L){
              return(args[[1L]])
            }
            # input checks
            classes <- lapply(args,class)
            if(length(unique(classes)) != 1L){
              stop("Inputs must be of the same SequenceDataFrame type.")
            }
            lengths <- vapply(args,function(a){sum(lengths(a))},integer(1))
            if(length(unique(lengths)) != 1L){
              stop("Inputs must have the same lengths.")
            }
            .check_ranges(args)
            .check_sequences(args)
            callNextMethod()
          }
)

#' @rdname SequenceData-class
#' @export
setMethod("rbind", "SequenceData",
          function(..., deparse.level = 1) 
          {
            args <- list(...)
            if(length(args) == 1L){
              return(args[[1L]])
            }
            # input checks
            classes <- lapply(args,class)
            if(length(unique(classes)) != 1L){
              stop("Inputs must be of the same SequenceDataFrame type.")
            }
            lengths <- vapply(args,function(a){ncol(unlist(a))},integer(1))
            if(length(unique(lengths)) != 1L){
              stop("Inputs must have the same width.")
            }
            .check_bamfiles(args)
            callNextMethod()
          }
)

setMethod("bindROWS", "SequenceData",
          function (x, objects = list(), use.names = TRUE, ignore.mcols = FALSE, 
                    check = TRUE) 
          {
            objects <- S4Vectors:::prepare_objects_to_bind(x, objects)
            all_objects <- c(list(x), objects)
            names <- unlist(lapply(all_objects,names))
            if(any(duplicated(names))){
              stop("Input must have unique names.")
            }
            .check_bamfiles(all_objects)
            callNextMethod(x, objects, use.names = use.names, 
                           ignore.mcols = ignore.mcols, check = check)
          }
)

setMethod("unlist", "SequenceData",
          function(x, recursive = TRUE, use.names = FALSE) 
          {
            callNextMethod(x, recursive = recursive, use.names = FALSE) 
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

.get_replicate_number <- function(conditions){
  control_rep <- seq_along(conditions[conditions == "control"])
  treated_rep <- seq_along(conditions[conditions == "treated"])
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
    m_data <- lapply(lapply(data, metadata),"[[","stats")
    data <- lapply(data, unlist)
    if(length(unique(lengths(data))) != 1L){
      stop("Data is of unequal length and cannot be coerced to a DataFrame.",
           call. = FALSE)
    }
    data <- S4Vectors::DataFrame(data)
  } else if(is(data[[1L]],"DataFrameList")) {
    m_data <- lapply(lapply(data, metadata),"$","stats")
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
  metadata(data) <- list(stats = m_data) 
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
  proto@minQuality <- .norm_min_quality(args, proto@minQuality)
  condition <- factor(names(bamfiles))
  replicate <- .get_replicate_number(condition)
  if(!.is_non_empty_string(proto@dataDescription)){
    stop("'dataDescription' must be a single non empty character value.")
  }
  if(is.null(proto@minQuality)){
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
  unlist_data <- 
    .SequenceDataFrame(class = gsub("SequenceData","",className),
                       df = unlist(data, use.names = FALSE),
                       ranges = unlist(ranges, use.names = FALSE),
                       sequence = unlist(sequences, use.names = FALSE),
                       replicate = replicate,
                       condition = condition,
                       bamfiles = bamfiles,
                       seqinfo = seqinfo)
  ans <- new(className, 
             minQuality = proto@minQuality,
             unlistData = unlist_data,
             partitioning = IRanges::PartitioningByEnd(data),
             metadata = metadata(unlist_data),
             ...)
  message("OK")
  ans
}

.SequenceData_settings <- data.frame(
  variable = c("max_depth",
               "minLength",
               "maxLength",
               "seqtype"),
  testFUN = c(".not_integer_bigger_than_10",
              ".not_integer_bigger_equal_than_zero_nor_na",
              ".not_integer_bigger_equal_than_one_nor_na",
              ".is_valid_nucleotide_seqtype"),
  errorValue = c(TRUE,
                 TRUE,
                 TRUE,
                 FALSE),
  errorMessage = c("'max_depth' must be integer with a value higher than 10L.",
                   "'minLength' must be integer with a value higher than 0L or NA.",
                   "'maxLength' must be integer with a value higher than 1L or NA.",
                   paste0("'seqtype' must be either '",seqtype(RNAString()) ,"' or '",seqtype(DNAString()) ,"'.")),
  stringsAsFactors = FALSE)

.get_SequenceData_args <- function(input){
  minQuality <- .norm_min_quality(input, NULL)
  max_depth <- 10000L # the default is 250, which is to small
  minLength <- NA_integer_
  maxLength <- NA_integer_
  seqtype <- seqtype(RNAString()) 
  args <- .norm_settings(input, .SequenceData_settings, max_depth, minLength,
                         maxLength, seqtype)
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
  sequences <- .load_sequences(sequences, grl, args)
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
.load_sequences <- function(sequences, grl, args){
  seq <- Biostrings::getSeq(sequences, unlist(grl))
  seq <- relist(unlist(seq),IRanges::PartitioningByWidth(sum(width(grl))))
  names(seq) <- names(grl)
  seqtype(seq) <- args[["seqtype"]]
  seq
}

# remove any elements, which are not in the seqinfo
.subset_by_seqinfo <- function(grl, seqinfo){
  grl <- grl[GenomicRanges::seqnames(grl) %in% GenomeInfoDb::seqnames(seqinfo)]
  grl <- grl[width(IRanges::PartitioningByWidth(grl)) != 0L]
  GenomeInfoDb::seqlevels(grl) <- GenomeInfoDb::seqlevelsInUse(grl)
  grl
}

################################################################################

#' @rdname SequenceData-class
#' @export
setGeneric( 
  name = "SequenceData",
  signature = c("annotation","sequences"),
  def = function(dataType, bamfiles, annotation, sequences, seqinfo, ...)
    standardGeneric("SequenceData")
) 

#' @rdname SequenceData-class
#' @export
setMethod("SequenceData",
          signature = c(annotation = "character", sequences = "character"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
#' @rdname SequenceData-class
#' @export
setMethod("SequenceData",
          signature = c(annotation = "character", sequences = "BSgenome"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
#' @rdname SequenceData-class
#' @export
setMethod("SequenceData",
          signature = c(annotation = "TxDb", sequences = "character"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
#' @rdname SequenceData-class
#' @export
setMethod("SequenceData",
          signature = c(annotation = "TxDb", sequences = "BSgenome"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
#' @rdname SequenceData-class
#' @export
setMethod("SequenceData",
          signature = c(annotation = "GRangesList", sequences = "character"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
#' @rdname SequenceData-class
#' @export
setMethod("SequenceData",
          signature = c(annotation = "GRangesList", sequences = "BSgenome"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
#' @rdname SequenceData-class
#' @export
setMethod("SequenceData",
          signature = c(annotation = "GFF3File", sequences = "BSgenome"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
#' @rdname SequenceData-class
#' @export
setMethod("SequenceData",
          signature = c(annotation = "GFF3File", sequences = "character"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
#' @rdname SequenceData-class
#' @export
setMethod("SequenceData",
          signature = c(annotation = "character", sequences = "FaFile"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
#' @rdname SequenceData-class
#' @export
setMethod("SequenceData",
          signature = c(annotation = "GFF3File", sequences = "FaFile"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
#' @rdname SequenceData-class
#' @export
setMethod("SequenceData",
          signature = c(annotation = "TxDb", sequences = "FaFile"),
          function(dataType, bamfiles, annotation, sequences, seqinfo, ...){
            .new_SequenceData(dataType, bamfiles, annotation, sequences,
                              seqinfo, ...)
          })
#' @rdname SequenceData-class
#' @export
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
          definition = function(x){bamfiles(unlist(x))})

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
      ends <- match(cumsum(width(partitioning)),cumsum(width(unlisted_ranges)))
      partitioning_relist <- IRanges::PartitioningByEnd(ends)
      if(length(x) != length(partitioning_relist)){
        stop("ranges could not be relisted.")
      }
      names(partitioning_relist) <- names(x)
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
          definition = function(x){seqinfo(unlist(x))})

#' @rdname SequenceData-functions
#' @export
setMethod(f = "sequences", 
          signature = signature(x = "SequenceData"),
          definition = function(x){relist(sequences(unlist(x)),x)})

#' @rdname SequenceData-functions
#' @export
setMethod(f = "seqtype", 
          signature = signature(x = "SequenceData"),
          definition = function(x){seqtype(unlist(x))})

#' @rdname SequenceData-functions
#' @export
setReplaceMethod(f = "seqtype", 
                 signature = signature(x = "SequenceData"),
                 definition = function(x, value){
                   unlisted_x <- unlist(x)
                   seqtype(unlisted_x) <- value
                   relist(unlisted_x,x)
                 })

#' @rdname SequenceData-functions
#' @export
setMethod(f = "dataType",
          signature = signature(x = "SequenceData"),
          definition = function(x){dataType(unlist(x))})

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
