#' @include RNAmodR.R
#' @include SequenceData-class.R
NULL

#' @name EndSequenceData-class
#' @aliases End5SequenceData End3SequenceData EndSequenceData
#'
#' @title End5SequenceData/End3SequenceData/EndSequenceData
#'
#' @description
#' The \code{End5SequenceData}/\code{End3SequenceData}/\code{EndSequenceData}
#' classes aggregate the counts of read ends at each position along a
#' transcript. \code{End5SequenceData}/\code{End3SequenceData} classes aggregate
#' either the 5'-end or 3'-end, the \code{EndSequenceData} aggregates both.
#'
#' All three classes contain one column per data file named using the following
#' naming convention \code{(end5/end3/end).condition.replicate}.
#'
#' \code{aggregate} calculates the mean and sd for samples in the \code{control}
#' and \code{treated} condition separatly.
#'
#' @param bamfiles,annotation,seqinfo,grl,sequences,param,args,... See
#' \code{\link[=SequenceData-class]{SequenceData}} and
#' \code{\link[=SequenceData-functions]{SequenceData-functions}}
#' @param x a \code{End5SequenceData}, \code{End3SequenceData} or
#' \code{EndSequenceData} object
#' @param name For \code{\link[=plotDataByCoord]{getDataTrack}}: a valid
#' transcript name. Must be a name of \code{ranges(x).}
#' @param condition For \code{\link{aggregate}}: condition for which the data
#' should be aggregated.
#' @param df,ranges,sequence,replicate inputs for creating a 
#' \code{SequenceDataFrame}. See 
#' \code{\link[=SequenceDataFrame-class]{SequenceDataFrame}}.
#'
#' @return a \code{End5SequenceData}, a \code{End3SequenceData} or a
#' \code{EndSequenceData} object
#'
#' @examples
#' # Construction of a End5SequenceData object
#' library(RNAmodR.Data)
#' library(rtracklayer)
#' annotation <- GFF3File(RNAmodR.Data.example.man.gff3())
#' sequences <- RNAmodR.Data.example.man.fasta()
#' files <- c(treated = RNAmodR.Data.example.wt.1())
#' e5sd <- End5SequenceData(files, annotation = annotation,
#'                         sequences = sequences)
NULL

#' @rdname EndSequenceData-class
#' @export
setClass(Class = "End5SequenceDataFrame",
         contains = "SequenceDFrame")
#' @rdname EndSequenceData-class
#' @export
End5SequenceDataFrame <- function(df, ranges, sequence, replicate,
                                  condition, bamfiles, seqinfo){
  .SequenceDataFrame("End5",df, ranges, sequence, replicate, condition,
                     bamfiles, seqinfo)
}
#' @rdname EndSequenceData-class
#' @export
setClass(Class = "End5SequenceData",
         contains = "SequenceData",
         slots = c(unlistData = "End5SequenceDataFrame"),
         prototype = list(unlistData = End5SequenceDataFrame(),
                          unlistType = "End5SequenceDataFrame",
                          minQuality = 5L,
                          dataDescription = "5'-end position data"))

#' @rdname EndSequenceData-class
#' @export
setClass(Class = "End3SequenceDataFrame",
         contains = "SequenceDFrame")
#' @rdname EndSequenceData-class
#' @export
End3SequenceDataFrame <- function(df, ranges, sequence, replicate, condition,
                                  bamfiles, seqinfo){
  .SequenceDataFrame("End3",df, ranges, sequence, replicate, condition,
                     bamfiles, seqinfo)
}
#' @rdname EndSequenceData-class
#' @export
setClass(Class = "End3SequenceData",
         contains = "SequenceData",
         slots = c(unlistData = "End3SequenceDataFrame"),
         prototype = list(unlistData = End3SequenceDataFrame(),
                          unlistType = "End3SequenceDataFrame",
                          minQuality = 5L,
                          dataDescription = "3'-end position data"))

#' @rdname EndSequenceData-class
#' @export
EndSequenceDataFrame <- function(df, ranges, sequence, replicate, condition,
                                 bamfiles, seqinfo){
  .SequenceDataFrame("End",df, ranges, sequence, replicate, condition,
                     bamfiles, seqinfo)
}
#' @rdname EndSequenceData-class
#' @export
setClass(Class = "EndSequenceDataFrame",
         contains = "SequenceDFrame")
#' @rdname EndSequenceData-class
#' @export
setClass(Class = "EndSequenceData",
         contains = "SequenceData",
         slots = c(unlistData = "EndSequenceDataFrame"),
         prototype = list(unlistData = EndSequenceDataFrame(),
                          unlistType = "EndSequenceDataFrame",
                          minQuality = 5L,
                          dataDescription = "read end position data (5' and 3')"))

#' @rdname EndSequenceData-class
#' @export
End5SequenceData <- function(bamfiles, annotation, sequences, seqinfo, ...){
  .new_SequenceData("End5", bamfiles = bamfiles, annotation = annotation,
                    sequences = sequences, seqinfo = seqinfo, ...)
}
#' @rdname EndSequenceData-class
#' @export
End3SequenceData <- function(bamfiles, annotation, sequences, seqinfo, ...){
  .new_SequenceData("End3", bamfiles = bamfiles, annotation = annotation,
                    sequences = sequences, seqinfo = seqinfo, ...)
}
#' @rdname EndSequenceData-class
#' @export
EndSequenceData <- function(bamfiles, annotation, sequences, seqinfo, ...){
  .new_SequenceData("End", bamfiles = bamfiles, annotation = annotation,
                    sequences = sequences, seqinfo = seqinfo, ...)
}

setSequenceDataCoercions("End5")
setSequenceDataCoercions("End3")
setSequenceDataCoercions("End")

# End5SequenceData ------------------------------------------------------------------

.summarize_to_position_data <- function(data, hits, names, strands, type){
  # get data for lapply
  starts <- BiocGenerics::start(data)
  ends <- BiocGenerics::end(data)
  rm(data)
  starts <- IRanges::IntegerList(split(starts[S4Vectors::queryHits(hits)],
                                       S4Vectors::subjectHits(hits)))
  ends <- IRanges::IntegerList(split(ends[S4Vectors::queryHits(hits)],
                                     S4Vectors::subjectHits(hits)))
  f <- as.integer(names(starts))
  strands <- strands[f]
  names(starts) <- names[f]
  names(ends) <- names[f]
  # aggregate pos of reads based on strand information
  if(type == "5prime"){
    data <- c(starts[strands == "+"],ends[strands == "-"])
    data <- data[names[f]]
  } else if(type == "3prime"){
    data <- c(starts[strands == "-"],ends[strands == "+"])
    data <- data[names[f]]
  } else if(type == "all"){
    data <- pc(starts, ends)
  } else if(type == "protected_ends"){
    data <- c(pc(ends[strands == "-"] - 1L,starts[strands == "-"]),
              pc(starts[strands == "+"] - 1L,ends[strands == "+"]))
    data <- data[names[f]]
  } else {
    stop("Something went wrong. Invalid type '", type, "'.")
  }
  data
}

#' @importFrom GenomicAlignments findOverlaps
.get_position_data_of_transcript_ends <- function(bamFile, grl, param,
                                                  type = c("5prime",
                                                           "3prime",
                                                           "all",
                                                           "protected_ends"),
                                                  args = list()){
  type <- match.arg(type)
  strands_u <- .get_strand_u_GRangesList(grl)
  data <- .load_bam_alignment_data(bamFile, param, grl, args)
  # get hits
  hits <- GenomicAlignments::findOverlaps(data, grl)
  # summarize pos of reads based on type
  data <- .summarize_to_position_data(data, hits, names(grl), strands_u, type)
  rm(hits)
  # get subsetting vectors
  f <- names(data)
  f_not_found <- names(grl)[!(names(grl) %in% f)]
  # calculate tables and add empty positions
  # also remove overhanging read data
  seqs <- .seqs_rl(grl)
  background <- relist(rep(0L,sum(lengths(seqs))),seqs)
  data <- data[data %in% seqs[f]]
  data <- data[lengths(data) > 0L]
  # update subsetting vectors
  f <- names(data)
  f_not_found <- names(grl)[!(names(grl) %in% f)]
  # calculated table
  data <- pc(data,seqs[f])
  # vectorized attempt does work for human transcriptome, since 
  # prod(dimensions) is bigger than .Machine$integer.max
  data <- IRanges::IntegerList(lapply(data,function(d){unname(table(d))}))
  data <- data - 1L
  # combine with empty positions
  data <- c(data,background[f_not_found])
  data <- data[match(names(grl),names(data))]
  data
}

#' @rdname EndSequenceData-class
#' @export
setMethod("getData",
          signature = c(x = "End5SequenceData",
                        bamfiles = "BamFileList",
                        grl = "GRangesList",
                        sequences = "XStringSet",
                        param = "ScanBamParam"),
          definition = function(x, bamfiles, grl, sequences, param, args){
            data <- lapply(bamfiles,
                           FUN = .get_position_data_of_transcript_ends,
                           grl = grl,
                           param = param,
                           type = "5prime",
                           args = args)
            names(data) <- rep("end5",length(data))
            data
          }
)

#' @rdname EndSequenceData-class
#' @export
setMethod("getData",
          signature = c(x = "End3SequenceData",
                        bamfiles = "BamFileList",
                        grl = "GRangesList",
                        sequences = "XStringSet",
                        param = "ScanBamParam"),
          definition = function(x, bamfiles, grl, sequences, param, args){
            data <- lapply(bamfiles,
                           FUN = .get_position_data_of_transcript_ends,
                           grl = grl,
                           param = param,
                           type = "3prime",
                           args = args)
            names(data) <- rep("end3",length(data))
            data
          }
)

#' @rdname EndSequenceData-class
#' @export
setMethod("getData",
          signature = c(x = "EndSequenceData",
                        bamfiles = "BamFileList",
                        grl = "GRangesList",
                        sequences = "XStringSet",
                        param = "ScanBamParam"),
          definition = function(x, bamfiles, grl, sequences, param, args){
            data <- lapply(bamfiles,
                           FUN = .get_position_data_of_transcript_ends,
                           grl = grl,
                           param = param,
                           type = "all",
                           args = args)
            names(data) <- rep("end",length(data))
            data
          }
)

# aggregation ------------------------------------------------------------------

#' @importFrom matrixStats rowSds
.aggregate_list_data_mean_sd <- function(x, condition){
  conditions <- conditions(x)
  f <- .subset_to_condition(conditions, condition)
  df <- as(unlist(x,use.names=FALSE),"DataFrame")[,f,drop=FALSE]
  conditions_u <- unique(conditions[f])
  # set up some base values. replicates is here the same as the number of
  # columns, since a list per replicate is assumed
  # get means
  means <- IRanges::NumericList(
    lapply(conditions_u,
           function(con){
             rowMeans(as.data.frame(df[,conditions[f] == con]),
                      na.rm = TRUE)
           }))
  names(means) <- paste0("means.", conditions_u)
  # get sds
  sds <- IRanges::NumericList(
    lapply(conditions_u,
           function(con){
             matrixStats::rowSds(as.matrix(df[,conditions[f] == con]),
                                 na.rm = TRUE)
           }))
  names(sds) <- paste0("sds.", conditions_u)
  # assemble data
  ans <- cbind(do.call(DataFrame, means),
               do.call(DataFrame, sds))
  ans <- relist(ans, IRanges::PartitioningByEnd(x))
  positions <- .seqs_rl_strand(ranges(x))
  rownames(ans) <- IRanges::CharacterList(positions)
  ans
}

#' @rdname EndSequenceData-class
#' @export
setMethod("aggregateData",
          signature = c(x = "End5SequenceData"),
          function(x, condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_list_data_mean_sd(x,condition)
          }
)

#' @rdname EndSequenceData-class
#' @export
setMethod("aggregateData",
          signature = c(x = "End3SequenceData"),
          function(x, condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_list_data_mean_sd(x,condition)
          }
)

#' @rdname EndSequenceData-class
#' @export
setMethod("aggregateData",
          signature = c(x = "EndSequenceData"),
          function(x, condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_list_data_mean_sd(x,condition)
          }
)

# data visualization -----------------------------------------------------------

RNAMODR_PLOT_SEQ_END_NAMES <- c("end5" = "mean(5'-ends)",
                                "end3" = "mean(3'-ends)",
                                "end" = "mean(ends)")

.clean_mcols_end <- function(seqdata){
  d <- mcols(seqdata@unlistData)
  d <- d[,grepl("means",colnames(d)),drop=FALSE]
  mcols(seqdata@unlistData) <- d
  seqdata
}

#' @rdname EndSequenceData-class
#' @export
setMethod(
  f = "getDataTrack",
  signature = signature(x = "EndSequenceData"),
  definition = function(x, name, ...) {
    args <- list(...)
    # DataTrack for sequence data
    seqdata <- .get_data_for_visualization(x, name)
    # clean meta data columns
    seqdata <- .clean_mcols_end(seqdata)
    seqdata <- unlist(seqdata)
    conditions <- unique(conditions(x))
    if("control" %in% conditions){
      d <- seqdata[,grepl("control",colnames(mcols(seqdata)))]
      colnames(mcols(d)) <- gsub(".control","",colnames(mcols(d)))
      dt.control <- Gviz::DataTrack(range = d,
                                    group = "means",
                                    name = paste0(RNAMODR_PLOT_SEQ_END_NAMES["end"],
                                                  "\ncontrol"),
                                    type = "histogram")
      Gviz::displayPars(dt.control)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.control)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.control)$col.axis <- "#000000"
      Gviz::displayPars(dt.control) <- args
      track <- list("End" = dt.control)
    }
    if("treated" %in% conditions){
      d <- seqdata[,grepl("treated",colnames(mcols(seqdata)))]
      colnames(mcols(d)) <- gsub(".treated","",colnames(mcols(d)))
      dt.treated <- Gviz::DataTrack(range = d,
                                    group = "means",
                                    name = paste0(RNAMODR_PLOT_SEQ_END_NAMES["end"],
                                                  "\ntreated"),
                                    type = "histogram")
      Gviz::displayPars(dt.treated)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.treated)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.treated)$col.axis <- "#000000"
      Gviz::displayPars(dt.treated) <- args
      track <- list("End" = dt.treated)
    }
    if(length(conditions) == 2L){
      track <- list("End" = dt.control,
                    "End" = dt.treated)
    }
    track
  }
)

#' @rdname EndSequenceData-class
#' @export
setMethod(
  f = "getDataTrack",
  signature = signature(x = "End5SequenceData"),
  definition = function(x, name, ...) {
    args <- list(...)
    # DataTrack for sequence data
    seqdata <- .get_data_for_visualization(x, name)
    # clean meta data columns
    seqdata <- .clean_mcols_end(seqdata)
    seqdata <- unlist(seqdata)
    conditions <- unique(conditions(x))
    if("control" %in% conditions){
      d <- seqdata[,grepl("control",colnames(mcols(seqdata)))]
      colnames(mcols(d)) <- gsub(".control","",colnames(mcols(d)))
      dt.control <- Gviz::DataTrack(range = d,
                                    group = "means",
                                    name = paste0(RNAMODR_PLOT_SEQ_END_NAMES["end5"],
                                                  "\ncontrol"),
                                    type = "histogram")
      Gviz::displayPars(dt.control)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.control)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.control)$col.axis <- "#000000"
      Gviz::displayPars(dt.control) <- args
      track <- list("End5" = dt.control)
    }
    if("treated" %in% conditions){
      d <- seqdata[,grepl("treated",colnames(mcols(seqdata)))]
      colnames(mcols(d)) <- gsub(".treated","",colnames(mcols(d)))
      dt.treated <- Gviz::DataTrack(range = d,
                                    group = "means",
                                    name = paste0(RNAMODR_PLOT_SEQ_END_NAMES["end5"],
                                                  "\ntreated"),
                                    type = "histogram")
      Gviz::displayPars(dt.treated)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.treated)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.treated)$col.axis <- "#000000"
      Gviz::displayPars(dt.treated) <- args
      track <- list("End5" = dt.treated)
    }
    if(length(conditions) == 2L){
      track <- list("End5" = dt.control,
                    "End5" = dt.treated)
    }
    track
  }
)

#' @rdname EndSequenceData-class
#' @export
setMethod(
  f = "getDataTrack",
  signature = signature(x = "End3SequenceData"),
  definition = function(x, name, ...) {
    args <- list(...)
    # DataTrack for sequence data
    seqdata <- .get_data_for_visualization(x, name)
    # clean meta data columns
    seqdata <- .clean_mcols_end(seqdata)
    seqdata <- unlist(seqdata)
    conditions <- unique(conditions(x))
    if("control" %in% conditions){
      d <- seqdata[,grepl("control",colnames(mcols(seqdata)))]
      colnames(mcols(d)) <- gsub(".control","",colnames(mcols(d)))
      dt.control <- Gviz::DataTrack(range = d,
                                    group = "means",
                                    name = paste0(RNAMODR_PLOT_SEQ_END_NAMES["end3"],
                                                  "\ncontrol"),
                                    type = "histogram")
      Gviz::displayPars(dt.control)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.control)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.control)$col.axis <- "#000000"
      Gviz::displayPars(dt.control) <- args
      track <- list("End3" = dt.control)
    }
    if("treated" %in% conditions){
      d <- seqdata[,grepl("treated",colnames(mcols(seqdata)))]
      colnames(mcols(d)) <- gsub(".treated","",colnames(mcols(d)))
      dt.treated <- Gviz::DataTrack(range = d,
                                    group = "means",
                                    name = paste0(RNAMODR_PLOT_SEQ_END_NAMES["end3"],
                                                  "\ntreated"),
                                    type = "histogram")
      Gviz::displayPars(dt.treated)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.treated)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.treated)$col.axis <- "#000000"
      Gviz::displayPars(dt.treated) <- args
      track <- list("End3" = dt.treated)
    }
    if(length(conditions) == 2L){
      track <- list("End3" = dt.control,
                    "End3" = dt.treated)
    }
    track
  }
)
