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
.get_position_data_of_transcript_ends_norm <- function(bamFile,
                                                       ranges,
                                                       param,
                                                       type = c("5prime",
                                                                "3prime"),
                                                       args = list()){
  type <- match.arg(type)
  parentRanges <- RNAmodR:::.get_parent_annotations(ranges)
  data <- .load_bam_alignment_data(bamFile,param,parentRanges,args)
  # factor for found and non found transcripts
  f <- as.integer(names(data))
  f_not_found <- as.integer(
    seq_along(parentRanges)[!(seq_along(parentRanges) %in% 
                                unique(names(data)))])
  # summarize pos of reads based on type
  strands <- as.character(BiocGenerics::strand(parentRanges)[f])
  enddata <- .summarize_to_position_data(data,strands,type)
  # tabulate the counts per position
  rl <- split(ranges(parentRanges),seq_along(parentRanges))
  enddata <- IRanges::IntegerList(mapply(
    function(d,r){
      bg <- table(start(r):end(r)) - 1
      d <- d[d >= start(r) & d <= end(r)]
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
    rl[f],
    SIMPLIFY = FALSE))
  names(enddata) <- f
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
    function(r){
      d <- table(start(r):end(r)) - 1
      as.integer(d)
    },
    rl[f_not_found],SIMPLIFY = FALSE))
  names(data_not_found) <- f_not_found
  # merge data with empty data and order based on factor numbers
  enddata <- c(enddata,data_not_found)
  enddata[is.na(enddata)] <- 0
  enddata <- enddata[order(as.integer(names(enddata))),]
  normTranscript <- c(normTranscript,data_not_found)
  normTranscript[is.na(normTranscript)] <- 0
  normTranscript <- normTranscript[order(as.integer(names(normTranscript))),]
  normOverlap <- c(normOverlap,data_not_found)
  normOverlap[is.na(normOverlap)] <- 0
  normOverlap <- normOverlap[order(as.integer(names(normOverlap))),]
  # name results based on transcript ID
  names(enddata) <- parentRanges$ID
  names(normTranscript) <- parentRanges$ID
  names(normOverlap) <- parentRanges$ID
  data <- IRanges::SplitDataFrameList(
    S4Vectors::DataFrame(ends = unlist(enddata),
                         norm.tx = unlist(normTranscript),
                         norm.ol = unlist(normOverlap)))
  data@partitioning <- normTranscript@partitioning
  data
}

#' @rdname NormEndSequenceData
#' @export
NormEnd5SequenceData <- function(bamfiles,
                                 fasta,
                                 gff,
                                 ...){
  args <- .get_mod_data_args(...)
  ans <- new("NormEnd5SequenceData",
             bamfiles,
             fasta,
             gff,
             args)
  ranges <- .load_annotation(ans@gff)
  sequences <- .load_transcript_sequences(ans@fasta,
                                          ranges)
  param <- .assemble_scanBamParam(ranges,
                                  ans@minQuality,
                                  ans@chromosomes)
  message("Loading normalized 5'-end position data from BAM files...")
  data <- lapply(ans@bamfiles,
                 FUN = .get_position_data_of_transcript_ends_norm,
                 ranges = ranges,
                 param = param,
                 type = "5prime",
                 args = args)
  
  
  names(data) <- paste0("norm.end5.",
                        names(ans@bamfiles),
                        ".",
                        seq_along(ans@bamfiles))
  .postprocess_read_data(ans,
                         data,
                         ranges,
                         sequences)
}

#' @rdname NormEndSequenceData
#' @export
NormEnd3SequenceData <- function(bamfiles,
                                 fasta,
                                 gff,
                                 ...){
  args <- .get_mod_data_args(...)
  ans <- new("NormEnd3SequenceData",
             bamfiles,
             fasta,
             gff,
             args)
  ranges <- .load_annotation(ans@gff)
  sequences <- .load_transcript_sequences(ans@fasta,
                                          ranges)
  param <- .assemble_scanBamParam(ranges,
                                  ans@minQuality,
                                  ans@chromosomes)
  message("Loading normalized 5'-end position data from BAM files...")
  data <- lapply(ans@bamfiles,
                 FUN = .get_position_data_of_transcript_ends_norm,
                 ranges = ranges,
                 param = param,
                 type = "3prime",
                 args = args)
  names(data) <- paste0("norm.end3.",
                        names(ans@bamfiles),
                        ".",
                        seq_along(ans@bamfiles))
  .postprocess_read_data(ans,
                         data,
                         ranges,
                         sequences)
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
          function(x,
                   condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_data_frame_mean_sd(x,condition)
          }
)

#' @name NormEndSequenceData
#' @export
setMethod("aggregate",
          signature = c(x = "NormEnd3SequenceData"),
          function(x,
                   condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_data_frame_mean_sd(x,condition)
          }
)