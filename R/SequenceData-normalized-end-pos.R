#' @include RNAmodR.R
#' @include SequenceData-class.R
NULL

#' @name NormEndSequenceData
#' @aliases NormEnd5SequenceData NormEnd3SequenceData
#' 
#' @title NormEndSequenceData
#' 
#' @description
#' title
NULL

#' @rdname NormEndSequenceData
#' @export
setClass(Class = "NormEnd5SequenceData",
         contains = "SequenceData",
         prototype = list(minQuality = 5L))

#' @rdname NormEndSequenceData
#' @export
setClass(Class = "NormEnd3SequenceData",
         contains = "SequenceData",
         prototype = list(minQuality = 5L))

# End5SequenceData ------------------------------------------------------------------

#'@importFrom reshape2 acast
.get_position_data_of_transcript_ends_norm <- function(bamFile, grl, param,
                                                       type = c("5prime",
                                                                "3prime"),
                                                       args = list()){
  type <- match.arg(type)
  strands_u <- .get_strand_u_GRangesList(grl)
  data <- .load_bam_alignment_data(bamFile, param, grl, args)
  # factor for found and non found transcripts
  f <- names(data)
  f_not_found <- names(grl)[!(names(grl) %in% names(data))]
  # summarize pos of reads based on type
  enddata <- .summarize_to_position_data(data, strands_u[f], type)
  # tabulate the counts per position
  seqs <- .seqs_rl(grl)
  enddata <- IRanges::IntegerList(mapply(
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
    enddata,
    seqs[f],
    SIMPLIFY = FALSE))
  # noralize against total number transcript or against the overlap per position
  normTranscript <- (enddata / BiocGenerics::lengths(data)) * 1000
  normTranscript <- IRanges::NumericList(lapply(normTranscript,unname))
  normOverlap <- IRanges::NumericList(mapply(
    function(d,end){
      gr <- GenomicRanges::GRanges(seqnames = as.character(unique(seqnames(d))),
                    ranges = ir <- IRanges(seq_along(end),width = 1),
                    strand = as.character(unique(strand(d))))
      end / GenomicRanges::countOverlaps(gr,d)
    },
    data,
    enddata,
    SIMPLIFY = FALSE))
  # calculate tables and add empty positions
  data_not_found <- IRanges::IntegerList(mapply(
    function(s){
      d <- table(s) - 1
      as.integer(d)
    },
    seqs[f_not_found],
    SIMPLIFY = FALSE))
  # merge data with empty data and order based on factor numbers
  enddata <- c(enddata,data_not_found)
  enddata@unlistData[is.na(enddata@unlistData)] <- 0L
  enddata <- enddata[match(names(grl),names(enddata))]
  normTranscript <- c(normTranscript,data_not_found)
  normTranscript@unlistData[is.na(normTranscript@unlistData)] <- 0
  normTranscript <- normTranscript[match(names(grl),names(normTranscript))]
  normOverlap <- c(normOverlap,data_not_found)
  normOverlap@unlistData[is.na(normOverlap@unlistData)] <- 0
  normOverlap <- normOverlap[match(names(grl),names(normOverlap))]
  # name results based on transcript ID
  data <- IRanges::SplitDataFrameList(
    S4Vectors::DataFrame(ends = unlist(enddata),
                         norm.tx = unlist(normTranscript),
                         norm.ol = unlist(normOverlap)))
  data@partitioning <- enddata@partitioning
  data
}

setMethod(".get_Data",
          signature = c(x = "NormEnd5SequenceData",
                        grl = "GRangesList",
                        sequences = "XStringSet",
                        param = "ScanBamParam"),
          definition = function(x, grl, sequences, param, args){
            message("Loading normalized 5'-end position data from BAM files ",
                    "... ", appendLF = FALSE)
            files <- bamfiles(x)
            data <- lapply(files,
                           FUN = .get_position_data_of_transcript_ends_norm,
                           grl = grl,
                           param = param,
                           type = "5prime",
                           args = args)
            names(data) <- paste0("norm.end5.",
                                  names(files),
                                  ".",
                                  seq_along(files))
            data
          }
)

setMethod(".get_Data",
          signature = c(x = "NormEnd3SequenceData",
                        grl = "GRangesList",
                        sequences = "XStringSet",
                        param = "ScanBamParam"),
          definition = function(x, grl, sequences, param, args){
            message("Loading normalized 3'-end position data from BAM files ",
                    "... ", appendLF = FALSE)
            files <- bamfiles(x)
            data <- lapply(files,
                           FUN = .get_position_data_of_transcript_ends_norm,
                           grl = grl,
                           param = param,
                           type = "3prime",
                           args = args)
            names(data) <- paste0("norm.end3.",
                                  names(files),
                                  ".",
                                  seq_along(files))
            data
          }
)

#' @rdname NormEndSequenceData
#' @export
NormEnd5SequenceData <- function(bamfiles, annotation, sequences, seqinfo, ...){
  SequenceData("NormEnd5", bamfiles = bamfiles, annotation = annotation,
               sequences = sequences, seqinfo = seqinfo, ...)
}
#' @rdname NormEndSequenceData
#' @export
NormEnd3SequenceData <- function(bamfiles, annotation, sequences, seqinfo, ...){
  SequenceData("NormEnd3", bamfiles = bamfiles, annotation = annotation,
               sequences = sequences, seqinfo = seqinfo, ...)
}

# aggregation ------------------------------------------------------------------

# aggregate
# - calculate mean per observation
# - calculate sd per observation
#' @importFrom matrixStats rowSds
.aggregate_data_frame_mean_sd <- function(x,
                                           condition){
  df <- .subset_to_condition(x@unlistData,
                             x@conditions,
                             condition)
  # set up some base values
  ncol <- ncol(df[,x@replicate == 1L,drop = FALSE])
  seqAdd <- seq.int(from = 0, to = ncol(df) - 1, by = ncol)
  colNames <- strsplit(colnames(df)[seq_len(ncol)],"\\.")
  colNames <- IRanges::CharacterList(colNames)[as.list(lengths(colNames))]
  # get means
  means <- IRanges::NumericList(lapply(seq_len(ncol),
                              function(i){
                                unname(rowMeans(as.data.frame(df[,i + seqAdd]),
                                         na.rm = TRUE))
                              }))
  names(means) <- paste0("means.",colNames)
  # get sds
  sds <- IRanges::NumericList(lapply(seq_len(ncol),
                            function(i){
                              unname(matrixStats::rowSds(as.matrix(df[,i + seqAdd]),
                                                  na.rm = TRUE))
                            }))
  names(sds) <- paste0("sds.",colNames)
  # merge data
  ans <- cbind(do.call(S4Vectors::DataFrame, means),
               do.call(S4Vectors::DataFrame, sds))
  ans <- IRanges::SplitDataFrameList(ans)
  ans@partitioning <- x@partitioning
  ans
}

#' @name NormEndSequenceData
#' @export
setMethod("aggregate",
          signature = c(x = "NormEnd5SequenceData"),
          function(x, condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_data_frame_mean_sd(x,condition)
          }
)

#' @name NormEndSequenceData
#' @export
setMethod("aggregate",
          signature = c(x = "NormEnd3SequenceData"),
          function(x, condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_data_frame_mean_sd(x,condition)
          }
)

# data visualization -----------------------------------------------------------
setMethod(
  f = ".dataTracks",
  signature = signature(x = "NormEnd5SequenceData",
                        data = "missing",
                        seqdata = "GRanges",
                        sequence = "XString"),
  definition = function(x, seqdata, sequence,  args) {
    requireNamespace("Gviz")
    browser()
  }
)

setMethod(
  f = ".dataTracks",
  signature = signature(x = "NormEnd3SequenceData",
                        data = "missing",
                        seqdata = "GRanges",
                        sequence = "XString"),
  definition = function(x, seqdata, sequence,  args) {
    requireNamespace("Gviz")
    browser()
  }
)
