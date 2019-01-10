#' @include RNAmodR.R
#' @include SequenceDataFrame-class.R
#' @include SequenceData-utils.R
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
                   bamfiles = "BamFileList",
                   fasta = "FaFile",
                   gff = "GFFFile",
                   minQuality = "integer",
                   chromosomes = "character",
                   unlistType = "character",
                   sequenceDataType = "factor",
                   replicate = "factor",
                   conditions = "factor"),
         prototype = list(ranges = GRangesList(),
                          sequences = DNAStringSet(),
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
            ranges_mcols <- mcols(object@ranges@unlistData, use.names=FALSE)
            ranges_nmc <- if (is.null(ranges_mcols)) 0L else ncol(ranges_mcols)
            cat(classNameForDisplay(object), " of length ", k, " containing\n",
                sep = "")
            if (k == 0L) {
              cat("<0 elements>\n")
            }
            cat(object@unlistType,"/",class(object@ranges@unlistData),
                "/",object@sequences@elementType," elements with ",
                data_nc, ifelse(data_nc == 1L, " data column", " data columns"),
                " and ",
                ranges_nmc, ifelse(ranges_nmc == 1L, " metadata column",
                                  " metadata columns"),
                "\n",
                sep = "")
            out_data <- NULL
            out_mdata <- NULL
            sep <- NULL
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
            if(!is.null(out_data) & !is.null(out_mdata)){
              sep <- matrix("|", dimnames = list("","|"))
            }
            out <- cbind(out_data,sep,out_mdata)
            print(out, quote = FALSE, right = TRUE)
          }
)
# validity ---------------------------------------------------------------------

.valid.SequenceData <- function(x){
  NULL
}
S4Vectors::setValidity2(Class = "SequenceData",.valid.SequenceData)

################################################################################
# list methods -----------------------------------------------------------------

extractElement <- function(x, i){
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
            extractElement(x, i)
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



lsubset_List_by_List <- function(x, i, value){
  lx <- length(x)
  li <- length(i)
  if (li == 0L) {
    ## Surprisingly, in that case, `[<-` on standard vectors does not
    ## even look at 'value'. So neither do we...
    return(x)
  }
  lv <- length(value)
  if (lv == 0L){
    stop("replacement has length zero")
  }
  value <- normalizeSingleBracketReplacementValue(value, x)
  if (is.null(names(i))) {
    if (li != lx){
      stop("when list-like subscript is unnamed, it must have the ",
           "length of list-like object to subset")
    }
    if (!is(x, "SimpleList")) {
      ## We'll try to take a fast path.
      if (is(i, "List")) {
        fast_path <- .select_fast_path(i, x)
      } else {
        i2 <- as(i, "List")
        i2_elttype <- elementType(i2)
        if (length(i2) == li && all(sapply(i, is, i2_elttype))) {
          fast_path <- .select_fast_path(i2, x)
          if (!is.na(fast_path))
            i <- i2
        } else {
          fast_path <- NA_character_
        }
      }
      if (!is.na(fast_path)){
        return(.fast_lsubset_List_by_List(x, i, value))  # fast path
      }
    }
    i2x <- seq_len(li)
  } else {
    if (is.null(names(x))){
      stop("cannot subset an unnamed list-like object ",
           "by a named list-like subscript")
    }
    i2x <- match(names(i), names(x))
    if (anyMissing(i2x)){
      stop("list-like subscript has names not in ",
           "list-like object to subset")
    }
    if (anyDuplicated(i2x)){
      stop("list-like subscript has duplicated names")
    }
  }
  value <- .adjust_value_length(value, li)
  ## Slow path (loops over the list elements of 'x').
  for (k in seq_len(li)){
    x[[i2x[k]]] <- replaceROWS(x[[i2x[k]]], i[[k]], value[[k]])
  }
  return(x)
}


#' @name SequenceData
#' @export
setReplaceMethod("[", "SequenceData",
                 function(x, i, j,..., value){
                   if(!is(value,"SequenceData")){
                     stop("invalid value. must be 'SequenceData'.")
                   }
                   browser()
                   if (!missing(j) || length(list(...)) > 0L){
                     stop("invalid subsetting")
                   }
                   if (!missing(i)){
                     return(lsubset_List_by_List(x, i, value))
                   }
                   callNextMethod(x, i, value=value)
                 }
)

##################
# looping


lapply_SequenceData <- function(X, FUN, ...){
  FUN <- match.fun(FUN)
  ans <- vector(mode="list", length=length(X))
  unlisted_X <- unlist(X, use.names=FALSE)
  X_partitioning <- PartitioningByEnd(X)
  X_elt_width <- width(X_partitioning)
  empty_idx <- which(X_elt_width == 0L)
  if (length(empty_idx) != 0L) 
    ans[empty_idx] <- list(FUN(extractROWS(unlisted_X, integer(0)), ...))
  non_empty_idx <- which(X_elt_width != 0L)
  if (length(non_empty_idx) == 0L)
    return(ans)
  X_elt_start <- start(X_partitioning)
  X_elt_end <- end(X_partitioning)
  slot_names <- 
  ans[non_empty_idx] <-
    lapply(non_empty_idx,
           function(i){
             df <- new(x@unlistType,
                       extractROWS(unlisted_X,
                                   IRanges(X_elt_start[i], X_elt_end[i])),
                       X@ranges[[i]],
                       X@sequences[[i]],
                       x@replicate,
                       x@conditions)
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
            i <- normalizeSingleBracketSubscript(i, x, as.NSBS=TRUE)
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


# Concatenation ----------------------------------------------------------------

.bindROWS <- function(...){
  args <- list(...)
  if (length(dim(args[[1L]])) >= 2L)
    return(rbind(...))
  concatenateObjects(args[[1L]], args[-1L])
}

.cbind_SequenceData_objects <- function(objects){
  browser()
  
  
}
  
# .cbind_SequenceData_objects <-
#   function(x, objects=list(), use.names=TRUE, ignore.mcols=FALSE, check=TRUE)
#   {
#     browser()
#     objects <- S4Vectors:::prepare_objects_to_concatenate(x, objects)
#     all_objects <- c(list(x), objects)
#     
#     ## 1. Take care of the parallel slots
#     
#     ## Call method for Vector objects to concatenate all the parallel slots
#     ## (only "elementMetadata" in the case of CompressedList) and stick them
#     ## into 'ans'. Note that the resulting 'ans' can be an invalid object
#     ## because its "elementMetadata" slot can be longer (i.e. have more rows)
#     ## than 'ans' itself so we use 'check=FALSE' to skip validation.
#     ans <- callNextMethod(x, objects, use.names=use.names,
#                           ignore.mcols=ignore.mcols,
#                           check=FALSE)
#     
#     ## 2. Take care of the non-parallel slots
#     
#     ## Concatenate the "unlistData" slots.
#     unlistData_list <- lapply(all_objects, slot, "unlistData")
#     ans_unlistData <- do.call(.bindROWS, unlistData_list)
#     
#     ## Concatenate the "partitioning" slots.
#     ans_breakpoints <- cumsum(unlist(lapply(all_objects, elementNROWS),
#                                      use.names=use.names))
#     ans_partitioning <- PartitioningByEnd(ans_breakpoints)
#     
#     BiocGenerics:::replaceSlots(ans, unlistData=ans_unlistData,
#                                 partitioning=ans_partitioning,
#                                 check=check)
#   }

.check_able_to_cbind <- function(objects){
  obj1 <- objects[[1L]]
  chk_length <- vapply(objects,
                       function(o){
                         all(lengths(o) == lengths(obj1))
                       },
                       logical(1))
  if(!all(chk_length)){
    stop("Lengths of SequenceData elements do not match.",
         call. = FALSE)
  }
  chk_names <- vapply(objects,
                      function(o){
                        all(names(o) == names(obj1))
                      },
                      logical(1))
  if(!all(chk_names)){
    stop("Names of SequenceData elements do not match.",
         call. = FALSE)
  }
  NULL
}

setMethod("cbind", "SequenceData",
          function(...){
            objects <- list(...)
            .check_able_to_cbind(objects)
            .cbind_SequenceData_objects(objects)
          }
)

setMethod("rbind", "SequenceData",
          function(...){
            arg1 <- list(...)[[1L]]
            stop("'rbind' not supported for ",class(arg1),".")
          }
)


# object creation --------------------------------------------------------------

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
    .Object@sequenceDataType <- factor(rep(class(.Object)[1L],length(bamfiles)))
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
  seq
}

.get_mod_data_args <- function(...){
  input <- list(...)
  max_depth <- 10000L # the default is 250, which is to small
  maxLength <- NA
  # for pileup
  if(!is.null(input[["max_depth"]])){
    max_depth <- input[["max_depth"]]
    if(!is.integer(max_depth) | max_depth <= 10L){
      stop("'max_depth' must be integer with a value higher than 10L.",
           call. = FALSE)
    }
  }
  # for protected end data
  if(!is.null(input[["maxLength"]])){
    maxLength <- input[["maxLength"]]
    if(!is.integer(maxLength) | maxLength <= 1L){
      stop("'maxLength' must be integer with a value higher than 1L or NA.",
           call. = FALSE)
    }
  }
  #
  args <- list(max_depth = max_depth,
               maxLength = maxLength)
  args
}

.postprocess_read_data <- function(x,
                                   data,
                                   ranges,
                                   sequences){
  conditionsFmultiplier <- 1L
  # work with the unlisted data and construct a CompressedSplitDataFrameList
  # from this
  if(is(data[[1L]],"IntegerList")){
    data <- lapply(data,
                   function(d){
                     d[order(names(d))]
                   })
    data <- lapply(data, unlist)
    data <- DataFrame(data)
    genes_ids <- unique(rownames(df))
  } else if(is(data[[1L]],"DataFrameList")) {
    data <- lapply(data, unlist, use.names = FALSE)
    ncols <- unique(unlist(lapply(data, ncol)))
    if(length(ncols) != 1L){
      stop("Something went wrong. Width of data should be constant.")
    }
    conditionsFmultiplier <- ncols
    data <- do.call(cbind,data)
  } else {
    stop("Something went wrong.")
  }
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
  x@sequenceDataType <- rep(x@sequenceDataType, each = conditionsFmultiplier)
  x@partitioning <- partitioning
  x@ranges <- ranges
  x@sequences <- sequences
  names(x) <- as.character(parentRanges$ID)
  if(any(names(x@ranges) != names(x@sequences)) || 
     any(names(x@ranges) != names(x))){
    stop("Something went wrong.")
  }
  x
}

# returns the position data for analysis as a list of data per replicate for
# individual transcript
.subset_to_intersecting_transcripts <- function(ID,
                                                data,
                                                class){
  data
}

# detect modifications in each file
.quantify_read_data_per_file <- function(bamFile,
                                         genome,
                                         ranges,
                                         dataFUN,
                                         param = NULL,
                                         args = list()){
  # load the bamfile
  bamData <- GenomicAlignments::readGAlignments(bamFile,
                                                param = param)
  # Total counts
  totalCounts <- Rsamtools::idxstatsBam(bamFile,
                                        param = param)
  totalCounts <- sum(totalCounts$mapped)
  # process result and split into chunks based on ranges
  IDs <- .get_IDs_from_scanBamParam(param)
  if(length(IDs) == 0){
    stop("Invalid ScanBamParam object.")
  }
  ranges_subset <- ranges[.get_unique_identifiers(ranges) %in% IDs,]
  hits <- IRanges::findOverlaps(bamData,
                                ranges_subset)
  bamData <- IRanges::extractList(bamData, 
                                  S4Vectors::split(
                                    S4Vectors::from(hits),
                                    as.factor(S4Vectors::to(hits))))
  f <- as.numeric(names(bamData))
  names(bamData) <- .get_unique_identifiers(ranges_subset)[f]
  # check for data
  if(length(bamData) == 0){
    warning("No reads detected in bam file '",
            bamFile,
            "'")
    return(NULL)
  }
  # split ranges and get corresponding sequences
  rangesList <- split(ranges_subset[f],
                      names(bamData))
  sequences <- Rsamtools::getSeq(genome,ranges_subset[f])
  # get the data per transcript
  data <- mapply(dataFUN,
                 bamData,
                 sequences,
                 rangesList,
                 MoreArgs = list(totalCounts = totalCounts),
                 SIMPLIFY = FALSE)
  # remove entries for transcript for which position data is insufficient
  data <- data[!vapply(data, is.null, logical(1))]
  if(length(data) == 0) return(NULL)
  return(data)
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
            function(x,
                     ...){
              stop("This functions needs to be implemented by '",class(x),"'.",
                   call. = FALSE)
            }
)
