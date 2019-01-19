#' @include RNAmodR.R
#' @include SequenceData-class.R
NULL

#' @name EndSequenceData
#' @aliases End5SequenceData End3SequenceData EndSequenceData
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

.load_bam_alignment_data <- function(bamFile,param,ranges,args){
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
  hits <- GenomicAlignments::findOverlaps(data,ranges)
  # split results per transcript
  data <- split(IRanges::subsetByOverlaps(data, ranges),
                S4Vectors::subjectHits(hits))
}

.summarize_to_position_data <- function(data,strands,type){
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
                                          c(starts[[i]],ends[[i]])
                                        }))
  } else if(type == "protected_ends"){
    data <- IRanges::IntegerList(lapply(seq_along(data),
                                        function(i){
                                          # offset applied to start to
                                          # sync the data on the position
                                          c(starts[[i]]-1,ends[[i]])
                                        }))
  } else {
    stop("Something went wrong. Invalid type '", type, "'.")
  }
  data
}

#'@importFrom reshape2 acast
.get_position_data_of_transcript_ends <- function(bamFile,
                                                  ranges,
                                                  param,
                                                  type = c("5prime",
                                                           "3prime",
                                                           "all",
                                                           "protected_ends"),
                                                  args = list()){
  type <- match.arg(type)
  parentRanges <- .get_parent_annotations(ranges)
  data <- .load_bam_alignment_data(bamFile,param,parentRanges,args)
  # factor for found and non found transcripts
  f <- as.integer(names(data))
  f_not_found <- as.integer(
    seq_along(parentRanges)[!(seq_along(parentRanges) %in% 
                                unique(names(data)))])
  # summarize pos of reads based on type
  strands <- as.character(BiocGenerics::strand(parentRanges)[f])
  data <- .summarize_to_position_data(data,strands,type)
  # calculate tables and add empty positions
  # also remove overhanging read data
  rl <- split(ranges(parentRanges),seq_along(parentRanges))
  data <- IRanges::IntegerList(mapply(
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
    data,
    rl[f],
    SIMPLIFY = FALSE))
  names(data) <- f
  # get data for empty transcripts
  data_not_found <- IRanges::IntegerList(mapply(
    function(r){
      d <- table(start(r):end(r)) - 1
      as.integer(d)
    },
    rl[f_not_found],
    SIMPLIFY = FALSE))
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
  args <- .get_mod_data_args(...)
  ans <- new("End5SequenceData",
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
  args <- .get_mod_data_args(...)
  ans <- new("End3SequenceData",
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
  args <- .get_mod_data_args(...)
  ans <- new("EndSequenceData",
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

#' @importFrom matrixStats rowSds
.aggregate_list_data_mean_sd <- function(x,
                                        condition){
  df <- .subset_to_condition(x@unlistData,
                             x@conditions,
                             condition)
  # set up some base values. replicates is here the same as the number of 
  # columns, since a list per replicate is assumed
  replicates <- unique(x@replicate)
  # get means
  means <- NumericList(lapply(replicates,
                              function(rep){
                                rowMeans(as.data.frame(df[,x@replicate == rep]),
                                         na.rm = TRUE)
                              }))
  names(means) <- paste0("means.",replicates)
  # get sds
  sds <- NumericList(lapply(replicates,
                            function(rep){
                              matrixStats::rowSds(as.matrix(df[,x@replicate == rep]),
                                                  na.rm = TRUE)
                            }))
  names(sds) <- paste0("sds.",replicates)
  # assemble data
  ans <- cbind(do.call(DataFrame, means),
               do.call(DataFrame, sds))
  ans <- SplitDataFrameList(ans)
  ans@partitioning <- x@partitioning
  ans
}

#' @name EndSequenceData
#' @export
setMethod("aggregate",
          signature = c(x = "EndSequenceData"),
          function(x,
                   condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_list_data_mean_sd(x,condition)
          }
)

#' @name EndSequenceData
#' @export
setMethod("aggregate",
          signature = c(x = "End5SequenceData"),
          function(x,
                   condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_list_data_mean_sd(x,condition)
          }
)


#' @name EndSequenceData
#' @export
setMethod("aggregate",
          signature = c(x = "End3SequenceData"),
          function(x,
                   condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_list_data_mean_sd(x,condition)
          }
)
