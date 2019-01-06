#' @include RNAmodR.R
#' @include SequenceData-class.R
NULL

#' @name EndSequenceData
#' 
#' @title EndSequenceData
#' 
#' @description
#' title
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



# End5SequenceData ------------------------------------------------------------------

#'@importFrom reshape2 acast
.get_position_data_of_transcript_ends <- function(bamFile,
                                                  ranges,
                                                  param,
                                                  type = "5prime",
                                                  args = list()){
  parentRanges <- .get_parent_annotations(ranges)
  # add some flags
  # bamWhat(param) <- c("flag","mapq")
  data <- GenomicAlignments::readGAlignments(bamFile, param = param)
  hits <- findOverlaps(data,parentRanges)
  # split results per transcript
  data <- split(subsetByOverlaps(data, parentRanges),
                subjectHits(hits))
  # factor for found and non found transcripts
  f <- as.integer(names(data))
  f_not_found <- as.integer(
    seq_along(parentRanges)[!(seq_along(parentRanges) %in% 
                                unique(subjectHits(hits)))])
  # get data for lapply
  starts <- start(data)
  ends <- end(data)
  widths <- width(parentRanges)
  strands <- as.character(strand(parentRanges)[f])
  # aggregate 5'-pos of reads based on strand information
  if(type == "5prime"){
    data <- IntegerList(lapply(seq_along(data),
                               function(i){
                                 if(strands[i] == "+"){
                                   starts[[i]]
                                 } else {
                                   ends[[i]]
                                 }
                               }))
  } else if(type == "3prime"){
    data <- IntegerList(lapply(seq_along(data),
                               function(i){
                                 if(strands[i] == "-"){
                                   starts[[i]]
                                 } else {
                                   ends[[i]]
                                 }
                               }))
  } else if(type == "all"){
    data <- IntegerList(lapply(seq_along(data),
                               function(i){
                                 c(starts[[i]],ends[[i]])
                               }))
  } else if(type == "protected_ends"){
    data <- IntegerList(lapply(seq_along(data),
                               function(i){
                                 c(starts[[i]]-1,ends[[i]])
                               }))
  } else {
    stop("Something went wrong. Invalid type '", type, "'.")
  }
  # calculate tables and add empty positions
  data <- IntegerList(mapply(
    function(d,w){
      bg <- table(seq_len(w)) - 1
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
    widths[f],SIMPLIFY = FALSE))
  names(data) <- f
  # get data for empty transcripts
  data_not_found <- IntegerList(mapply(
    function(w){
      d <- table(seq_len(w)) - 1
      as.integer(d)
    },
    widths[f_not_found],SIMPLIFY = FALSE))
  names(data_not_found) <- f_not_found
  # merge and order based on factor numbers
  data <- c(data,data_not_found)
  data <- data[order(as.integer(names(data))),]
  # name results based on transcript ID
  names(data) <- parentRanges$ID
  data
}

#' @rdname EndSequenceData
#' @export
End5SequenceData <- function(bamfiles,
                             fasta,
                             gff,
                             ...){
  ans <- new("End5SequenceData",
             bamfiles,
             fasta,
             gff,
             ...)
  args <- .get_mod_data_args(...)
  ranges <- .load_annotation(ans@gff)
  sequences <- .load_transcript_sequences(ans@fasta,
                                          ranges)
  param <- .assemble_scanBamParam(ranges,
                                  ans@minQuality,
                                  ans@chromosomes)
  message("Loading 5'-end position data from BAM files...")
  data <- lapply(ans@bamfiles,
                 FUN = .get_position_data_of_transcript_ends,
                 ranges = ranges,
                 param = param,
                 type = "5prime",
                 args = args)
  names(data) <- paste0("end5.",
                        names(ans@bamfiles),
                        ".",
                        seq_along(ans@bamfiles))
  .postprocess_read_data(ans,
                         data,
                         ranges,
                         sequences)
}

#' @rdname EndSequenceData
#' @export
End3SequenceData <- function(bamfiles,
                             fasta,
                             gff,
                             ...){
  ans <- new("End3SequenceData",
             bamfiles,
             fasta,
             gff,
             ...)
  args <- .get_mod_data_args(...)
  ranges <- .load_annotation(ans@gff)
  sequences <- .load_transcript_sequences(ans@fasta,
                                          ranges)
  param <- .assemble_scanBamParam(ranges,
                                  ans@minQuality,
                                  ans@chromosomes)
  message("Loading 3'-end position data from BAM files...")
  data <- lapply(ans@bamfiles,
                 FUN = .get_position_data_of_transcript_ends,
                 ranges = ranges,
                 param = param,
                 type = "3prime",
                 args = args)
  names(data) <- paste0("end3.",
                        names(ans@bamfiles),
                        ".",
                        seq_along(ans@bamfiles))
  .postprocess_read_data(ans,
                         data,
                         ranges,
                         sequences)
}

#' @rdname EndSequenceData
#' @export
EndSequenceData <- function(bamfiles,
                            fasta,
                            gff,
                            ...){
  ans <- new("EndSequenceData",
             bamfiles,
             fasta,
             gff,
             ...)
  args <- .get_mod_data_args(...)
  ranges <- .load_annotation(ans@gff)
  sequences <- .load_transcript_sequences(ans@fasta,
                                          ranges)
  param <- .assemble_scanBamParam(ranges,
                                  ans@minQuality,
                                  ans@chromosomes)
  message("Loading read end position data (5' and 3') from BAM files...")
  data <- lapply(ans@bamfiles,
                 FUN = .get_position_data_of_transcript_ends,
                 ranges = ranges,
                 param = param,
                 type = "all",
                 args = args)
  names(data) <- paste0("end.",
                        names(ans@bamfiles),
                        ".",
                        seq_along(ans@bamfiles))
  .postprocess_read_data(ans,
                         data,
                         ranges,
                         sequences)
}

# aggregation ------------------------------------------------------------------

.aggregate_end_data_mean_sd <- function(x,
                                        condition){
  browser()
  condition <- tolower(match.arg(condition))
  if(is(x,"CompressedSplitDataFrameList")){
    df <- x@unlistData
  } else {
    df <- x
  }
  if(condition != "both"){
    df <- df[,x@conditions == condition, drop = FALSE]
    if(ncol(df) == 0L){
      stop("No data for condition '",condition,"' found.")
    }
  }
  replicates <- unique(x@replicate)
  means <- NumericList(lapply(replicates,
                              function(rep){
                                rowMeans(as.data.frame(df[,x@replicate == rep]),
                                         na.rm = TRUE)
                              }))
  names(means) <- paste0("means.",replicates)
  sds <- NumericList(lapply(replicates,
                            function(rep){
                              matrixStats::rowSds(as.matrix(df[,x@replicate == rep]),
                                                  na.rm = TRUE)
                            }))
  names(sds) <- paste0("sds.",replicates)
  ans <- cbind(do.call(DataFrame, means),
               do.call(DataFrame, sds))
  if(is(x,"CompressedSplitDataFrameList")){
    ans <- SplitDataFrameList(ans)
    ans@partitioning <- x@partitioning
  }
  ans
}

#' @name EndSequenceData
#' @export
setMethod("aggregate",
          signature = c(x = "EndSequenceData"),
          function(x,
                   condition = c("Both","Treated","Control")){
            .aggregate_end_data_mean_sd(x,condition)
          }
)

#' @name EndSequenceData
#' @export
setMethod("aggregate",
          signature = c(x = "End5SequenceData"),
          function(x,
                   condition = c("Both","Treated","Control")){
            .aggregate_end_data_mean_sd(x,condition)
          }
)


#' @name EndSequenceData
#' @export
setMethod("aggregate",
          signature = c(x = "End3SequenceData"),
          function(x,
                   condition = c("Both","Treated","Control")){
            .aggregate_end_data_mean_sd(x,condition)
          }
)
