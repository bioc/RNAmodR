#' @include RNAmodR.R
#' @include SequenceData-class.R
NULL

#' @name EndSequenceData
#' @aliases End5SequenceData End3SequenceData EndSequenceData
#' 
#' @title End5SequenceData/End3SequenceData/EndSequenceData
#' 
#' @description
#' The \code{End5SequenceData}/\code{End3SequenceData}/\code{EndSequenceData}
#' aggregate the counts of read ends at each position along transcript. Whereas
#' the first aggregate either the 5'-end or 3'-end, the \code{EndSequenceData}
#' aggregates both.
#' 
#' \code{aggregate} calculates the mean and sd for samples in the \code{control}
#' and \code{treated} condition serparatly.
#' 
#' 
#' @examples
#' # Construct a End5SequenceData object
#' annotation <- system.file("extdata","example1.gff3",package = "RNAmodR.Data")
#' sequences <- system.file("extdata","example1.fasta",package = "RNAmodR.Data")
#' files <- c(control = system.file("extdata","example_wt_1.bam",
#'                                  package = "RNAmodR.Data"),
#'            treated = system.file("extdata","example_wt_2.bam",
#'                                  package = "RNAmodR.Data"))
#' e5sd <- End5SequenceData(files, annotation = annotation,
#'                         sequences = sequences)
#' # aggregate data
#' aggregate(e5sd)
NULL

#' @rdname EndSequenceData
#' @export
setClass(Class = "End5SequenceData",
         contains = "SequenceData",
         prototype = list(minQuality = 5L))

#' @rdname EndSequenceData
#' @export
setClass(Class = "End3SequenceData",
         contains = "SequenceData",
         prototype = list(minQuality = 5L))

#' @rdname EndSequenceData
#' @export
setClass(Class = "EndSequenceData",
         contains = "SequenceData",
         prototype = list(minQuality = 5L))

#' @rdname EndSequenceData
#' @export
End5SequenceData <- function(bamfiles, annotation, sequences, seqinfo, ...){
  SequenceData("End5", bamfiles = bamfiles, annotation = annotation,
               sequences = sequences, seqinfo = seqinfo, ...)
}
#' @rdname EndSequenceData
#' @export
End3SequenceData <- function(bamfiles, annotation, sequences, seqinfo, ...){
  SequenceData("End3", bamfiles = bamfiles, annotation = annotation,
               sequences = sequences, seqinfo = seqinfo, ...)
}
#' @rdname EndSequenceData
#' @export
EndSequenceData <- function(bamfiles, annotation, sequences, seqinfo, ...){
  SequenceData("End", bamfiles = bamfiles, annotation = annotation,
               sequences = sequences, seqinfo = seqinfo, ...)
}


# End5SequenceData ------------------------------------------------------------------

.load_bam_alignment_data <- function(bamFile, param, grl, args){
  data <- GenomicAlignments::readGAlignments(bamFile, param = param)
  if(length(data) == 0L){
    stop("No reads found in data.", call. = FALSE)
  }
  # apply length cut off if set
  if(!is.na(args[["maxLength"]])){
    data <- data[width(data) <= args[["maxLength"]],]
  }
  if(!is.na(args[["minLength"]])){
    data <- data[width(data) >= args[["minLength"]],]
  }
  if(length(data) == 0L){
    stop("No reads found in data with read length equal or between 'minLength'",
         " (",args[["minLength"]]," nt) and 'maxLength' (",args[["maxLength"]],
         " nt).", call. = FALSE)
  }
  hits <- GenomicAlignments::findOverlaps(data, grl)
  # split results per transcript
  data <- split(data[S4Vectors::queryHits(hits)],
                S4Vectors::subjectHits(hits))
  data
}

.summarize_to_position_data <- function(data, strands, type){
  # get data for lapply
  starts <- BiocGenerics::start(data)
  ends <- BiocGenerics::end(data)
  # aggregate pos of reads based on strand information
  if(type == "5prime"){
    data <- IRanges::IntegerList(lapply(seq_along(data),
                                        function(i){
                                          if(strands[i] == "+"){
                                            starts[[i]]
                                          } else {
                                            ends[[i]]
                                          }
                                        }))
  } else if(type == "3prime"){
    data <- IRanges::IntegerList(lapply(seq_along(data),
                                        function(i){
                                          if(strands[i] == "-"){
                                            starts[[i]]
                                          } else {
                                            ends[[i]]
                                          }
                                        }))
  } else if(type == "all"){
    data <- IRanges::IntegerList(lapply(seq_along(data),
                                        function(i){
                                          c(starts[[i]], ends[[i]])
                                        }))
  } else if(type == "protected_ends"){
    data <- IRanges::IntegerList(lapply(seq_along(data),
                                        function(i){
                                          # offset applied to start to
                                          # sync the data on the position
                                          c(starts[[i]] - 1L, ends[[i]])
                                        }))
  } else {
    stop("Something went wrong. Invalid type '", type, "'.")
  }
  names(data) <- names(starts)
  data
}

#'@importFrom reshape2 acast
.get_position_data_of_transcript_ends <- function(bamFile, grl, param,
                                                  type = c("5prime",
                                                           "3prime",
                                                           "all",
                                                           "protected_ends"),
                                                  args = list()){
  type <- match.arg(type)
  strands_u <- .get_strand_u_GRangesList(grl)
  data <- .load_bam_alignment_data(bamFile, param, grl, args)
  # factor for found and non found transcripts
  f <- names(data)
  f_not_found <- names(grl)[!(names(grl) %in% names(data))]
  # summarize pos of reads based on type
  data <- .summarize_to_position_data(data, strands_u[f], type)
  # calculate tables and add empty positions
  # also remove overhanging read data
  seqs <- .seqs_rl(grl)
  data <- IRanges::IntegerList(mapply(
    function(d,s){
      bg <- table(s) - 1
      d <- d[d %in% s]
      d <- table(d)
      d <- d[as.integer(names(d)) > 0L]
      d <- reshape2::acast(data.frame(pos = as.integer(c(names(bg),names(d))),
                                      count = as.integer(c(bg,d))),
                           pos ~ .,
                           value.var = "count",
                           fun.aggregate = sum)
      as.integer(d)
    },
    data,
    seqs[f],
    SIMPLIFY = FALSE))
  # get data for empty transcripts
  data_not_found <- IRanges::IntegerList(mapply(
    function(s){
      d <- table(s) - 1
      as.integer(d)
    },
    seqs[f_not_found],
    SIMPLIFY = FALSE))
  # merge and order
  data <- c(data,data_not_found)
  data <- data[match(names(grl),names(data))]
  data
}

#' @rdname RNAmodR-internals
setMethod(".getData",
          signature = c(x = "End5SequenceData",
                        grl = "GRangesList",
                        sequences = "XStringSet",
                        param = "ScanBamParam"),
          definition = function(x, grl, sequences, param, args){
            message("Loading 5'-end position data from BAM files ... ",
                    appendLF = FALSE)
            files <- bamfiles(x)
            data <- lapply(files,
                           FUN = .get_position_data_of_transcript_ends,
                           grl = grl,
                           param = param,
                           type = "5prime",
                           args = args)
            names(data) <- paste0("end5.", x@condition, ".", x@replicate)
            data
          }
)

#' @rdname RNAmodR-internals
setMethod(".getData",
          signature = c(x = "End3SequenceData",
                        grl = "GRangesList",
                        sequences = "XStringSet",
                        param = "ScanBamParam"),
          definition = function(x, grl, sequences, param, args){
            message("Loading 3'-end position data from BAM files ... ",
                    appendLF = FALSE)
            files <- bamfiles(x)
            data <- lapply(files,
                           FUN = .get_position_data_of_transcript_ends,
                           grl = grl,
                           param = param,
                           type = "3prime",
                           args = args)
            names(data) <- paste0("end3.", x@condition, ".", x@replicate)
            data
          }
)

#' @rdname RNAmodR-internals
setMethod(".getData",
          signature = c(x = "EndSequenceData",
                        grl = "GRangesList",
                        sequences = "XStringSet",
                        param = "ScanBamParam"),
          definition = function(x, grl, sequences, param, args){
            message("Loading read end position data (5' and 3') from BAM ",
                    "files ... ", appendLF = FALSE)
            files <- bamfiles(x)
            data <- lapply(files,
                           FUN = .get_position_data_of_transcript_ends,
                           grl = grl,
                           param = param,
                           type = "all",
                           args = args)
            names(data) <- paste0("end.", x@condition, ".", x@replicate)
            data
          }
)

# aggregation ------------------------------------------------------------------

#' @importFrom matrixStats rowSds
.aggregate_list_data_mean_sd <- function(x, condition){
  f <- .subset_to_condition(x@condition, condition)
  df <- x@unlistData[f]
  conditions <- unique(x@condition[f])
  # set up some base values. replicates is here the same as the number of 
  # columns, since a list per replicate is assumed
  # get means
  means <- IRanges::NumericList(
    lapply(conditions,
           function(con){
             rowMeans(as.data.frame(df[,x@condition[f] == con]),
                      na.rm = TRUE)
           }))
  names(means) <- paste0("means.", conditions)
  # get sds
  sds <- IRanges::NumericList(
    lapply(conditions,
           function(con){
             matrixStats::rowSds(as.matrix(df[,x@condition[f] == con]),
                                 na.rm = TRUE)
           }))
  names(sds) <- paste0("sds.", conditions)
  # assemble data
  ans <- cbind(do.call(DataFrame, means),
               do.call(DataFrame, sds))
  ans <- SplitDataFrameList(ans)
  ans@partitioning <- x@partitioning
  ans
}

#' @rdname EndSequenceData
#' @export
setMethod("aggregate",
          signature = c(x = "End5SequenceData"),
          function(x, condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_list_data_mean_sd(x,condition)
          }
)

#' @rdname EndSequenceData
#' @export
setMethod("aggregate",
          signature = c(x = "End3SequenceData"),
          function(x, condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_list_data_mean_sd(x,condition)
          }
)

#' @rdname EndSequenceData
#' @export
setMethod("aggregate",
          signature = c(x = "EndSequenceData"),
          function(x, condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_list_data_mean_sd(x,condition)
          }
)

# data visualization -----------------------------------------------------------

#' @rdname RNAmodR-internals
setMethod(
  f = ".dataTracks",
  signature = signature(x = "EndSequenceData",
                        data = "missing",
                        seqdata = "GRanges",
                        sequence = "XString"),
  definition = function(x, seqdata, sequence,  args) {
    requireNamespace("Gviz")
    browser()
  }
)

#' @rdname RNAmodR-internals
setMethod(
  f = ".dataTracks",
  signature = signature(x = "End5SequenceData",
                        data = "missing",
                        seqdata = "GRanges",
                        sequence = "XString"),
  definition = function(x, seqdata, sequence,  args) {
    requireNamespace("Gviz")
    browser()
  }
)

#' @rdname RNAmodR-internals
setMethod(
  f = ".dataTracks",
  signature = signature(x = "End3SequenceData",
                        data = "missing",
                        seqdata = "GRanges",
                        sequence = "XString"),
  definition = function(x, seqdata, sequence,  args) {
    requireNamespace("Gviz")
    browser()
  }
)
