#' @include RNAmodR.R
#' @include SequenceData-class.R
NULL

#' @name NormEndSequenceData-class
#' @aliases NormEnd5SequenceData NormEnd3SequenceData
#' 
#' @title NormEnd5SequenceData/NormEnd3SequenceData
#' 
#' @description
#' The \code{NormEnd5SequenceData}/\code{NormEnd3SequenceData}
#' aggregate the counts of read ends (Either 5' or 3') at each position along a
#' transcript. In addition, the number of counts are then normalized to the
#' length of the transcript and to the overlapping reads.
#' 
#' Both classes contain three columns per data file named using the
#' following naming convention \code{(normend5/normend3).condition.replicate}.
#' The three columns are distinguished by additional identifiers \code{ends},
#' \code{norm.tx} and \code{norm.ol}.
#' 
#' \code{aggregate} calculates the mean and sd for samples in the \code{control}
#' and \code{treated} condition separatly. Similar to the stored results for 
#' each of the two conditions six columns are returned (three for mean and sd 
#' each) ending in \code{ends}, \code{tx} and \code{ol}.
#' 
#' @param bamfiles,annotation,seqinfo,grl,sequences,param,args,... See 
#' \code{\link[=SequenceData-class]{SequenceData}} and
#' \code{\link[=SequenceData-functions]{SequenceData-functions}}
#' @param x a \code{CoverageSequenceData}
#' @param name For \code{\link[=plotDataByCoord]{getDataTrack}}: a valid 
#' transcript name. Must be a name of \code{ranges(x)}
#' @param condition For \code{\link{aggregate}}: condition for which the data 
#' should be aggregated.
#' @param df,ranges,sequence,replicate inputs for creating a 
#' \code{SequenceDataFrame}. See 
#' \code{\link[=SequenceDataFrame-class]{SequenceDataFrame}}.
#' 
#' @return a \code{NormEnd5SequenceData} or \code{NormEnd3SequenceData} object
#' 
#' @examples
#' # Construction of a NormEnd5SequenceData object
#' \dontrun{
#' library(RNAmodR.Data)
#' library(rtracklayer)
#' annotation <- GFF3File(RNAmodR.Data.example.man.gff3())
#' sequences <- RNAmodR.Data.example.man.fasta()
#' files <- c(treated = RNAmodR.Data.example.wt.1())
#' ne5sd <- NormEnd5SequenceData(files, annotation = annotation,
#'                               sequences = sequences)
#' }
NULL

#' @rdname NormEndSequenceData-class
#' @export
setClass(Class = "NormEnd5SequenceDataFrame",
         contains = "SequenceDataFrame")
#' @rdname NormEndSequenceData-class
#' @export
NormEnd5SequenceDataFrame <- function(df, ranges, sequence, replicate,
                                      condition, bamfiles, seqinfo){
  .SequenceDataFrame("NormEnd5",df, ranges, sequence, replicate, condition,
                     bamfiles, seqinfo)
}
#' @rdname NormEndSequenceData-class
#' @export
setClass(Class = "NormEnd5SequenceData",
         contains = "SequenceData",
         slots = c(unlistData = "NormEnd5SequenceDataFrame"),
         prototype = list(unlistData = NormEnd5SequenceDataFrame(),
                          unlistType = "NormEnd5SequenceDataFrame",
                          minQuality = 5L,
                          dataDescription = "normalized 5'-end position data"))

#' @rdname NormEndSequenceData-class
#' @export
setClass(Class = "NormEnd3SequenceDataFrame",
         contains = "SequenceDataFrame")
#' @rdname NormEndSequenceData-class
#' @export
NormEnd3SequenceDataFrame <- function(df, ranges, sequence, replicate,
                                      condition, bamfiles, seqinfo){
  .SequenceDataFrame("NormEnd3",df, ranges, sequence, replicate, condition,
                     bamfiles, seqinfo)
}
#' @rdname NormEndSequenceData-class
#' @export
setClass(Class = "NormEnd3SequenceData",
         contains = "SequenceData",
         slots = c(unlistData = "NormEnd3SequenceDataFrame"),
         prototype = list(unlistData = NormEnd3SequenceDataFrame(),
                          unlistType = "NormEnd3SequenceDataFrame",
                          minQuality = 5L,
                          dataDescription = "normalized 3'-end position data"))

#' @rdname NormEndSequenceData-class
#' @export
NormEnd5SequenceData <- function(bamfiles, annotation, sequences, seqinfo, ...){
  .new_SequenceData("NormEnd5", bamfiles = bamfiles, annotation = annotation,
                    sequences = sequences, seqinfo = seqinfo, ...)
}
#' @rdname NormEndSequenceData-class
#' @export
NormEnd3SequenceData <- function(bamfiles, annotation, sequences, seqinfo, ...){
  .new_SequenceData("NormEnd3", bamfiles = bamfiles, annotation = annotation,
                    sequences = sequences, seqinfo = seqinfo, ...)
}

setSequenceDataCoercions("NormEnd5")
setSequenceDataCoercions("NormEnd3")

# summary ----------------------------------------------------------------------

.get_summary_MultiColSequenceData <- function(object){
  bamfilesstats <- .get_bamfiles_stats(object)
  n <- ncol(bamfilesstats)
  datastats <- .get_data_stats(object)
  ndt <- ncol(datastats)
  ncols <- ndt / n
  datastats <- lapply(seq_len(ncols),
                      function(i){
                        d <- datastats[,seq.int(i,ndt,ncols),drop=FALSE]
                        coln <- colnames(object@unlistData)[i]
                        coln <- strsplit(coln,"\\.")[[1]]
                        coln <- coln[length(coln)]
                        rownames(d) <- paste0("Data.",coln,
                                              gsub("Data","",rownames(d)))
                        colnames(d) <- seq_len(ncol(d))
                        d
                      })
  .merge_summary_data(bamfilesstats,Reduce(rbind,datastats))
}

setMethod("summary",
          signature = "NormEnd5SequenceData",
          function(object){
            .get_summary_MultiColSequenceData(object)
          })

setMethod("summary",
          signature = "NormEnd3SequenceData",
          function(object){
            .get_summary_MultiColSequenceData(object)
          })


# End5SequenceData ------------------------------------------------------------------

#' @importFrom GenomicAlignments findOverlaps
.get_position_data_of_transcript_ends_norm <- function(bamFile, grl, param,
                                                       type = c("5prime",
                                                                "3prime"),
                                                       args = list()){
  type <- match.arg(type)
  strands_u <- .get_strand_u_GRangesList(grl)
  data <- .load_bam_alignment_data(bamFile, param, grl, args)
  # get hits
  hits <- GenomicAlignments::findOverlaps(data, grl)
  # summarize pos of reads based on type
  enddata <- .summarize_to_position_data(data, hits, names(grl), strands_u,
                                         type)
  rm(hits)
  # get subsetting vectors
  f <- names(enddata)
  f_not_found <- names(grl)[!(names(grl) %in% f)]
  # tabulate the counts per position
  seqs <- .seqs_rl_strand(grl, force_continous = TRUE)
  background <- relist(rep(0L,sum(lengths(seqs))),seqs)
  enddata <- enddata[enddata %in% seqs[f]]
  enddata <- enddata[lengths(enddata) > 0L]
  # update subsetting vectors
  f <- names(enddata)
  f_not_found <- names(grl)[!(names(grl) %in% f)]
  # calculated table
  enddata <- pc(enddata,seqs[f])
  # vectorized attempt does work for human transcriptome, since 
  # prod(dimensions) is bigger than .Machine$integer.max
  enddata <- IRanges::IntegerList(lapply(enddata,function(d){unname(table(d))}))
  enddata <- enddata - 1L
  # noralize against total number transcript or against the overlap per position
  normTranscript <- (enddata / unname(sum(enddata))) * 1000
  coverage <- .get_coverage_from_GA(data, grl[f])
  if(any(names(coverage) != names(enddata))){
    stop("")
  }
  normOverlap <- enddata / coverage
  # merge data with empty data and order based on factor numbers
  enddata <- c(enddata,background[f_not_found])
  enddata@unlistData[is.na(enddata@unlistData)] <- 0L
  enddata <- enddata[match(names(grl),names(enddata))]
  normTranscript <- c(normTranscript,background[f_not_found])
  normTranscript@unlistData[is.na(normTranscript@unlistData)] <- 0
  normTranscript <- normTranscript[match(names(grl),names(normTranscript))]
  normOverlap <- c(normOverlap,background[f_not_found])
  normOverlap@unlistData[is.na(normOverlap@unlistData)] <- 0
  normOverlap@unlistData[is.infinite(normOverlap@unlistData)] <- 0
  normOverlap <- normOverlap[match(names(grl),names(normOverlap))]
  rm(background)
  # name results based on transcript ID
  df <- S4Vectors::DataFrame(ends = unlist(enddata),
                             norm.tx = unlist(normTranscript),
                             norm.ol = unlist(normOverlap))
  ans <- relist(df, IRanges::PartitioningByWidth(enddata))
  rownames(ans) <- IRanges::CharacterList(seqs)
  ans
}

#' @rdname NormEndSequenceData-class
#' @export
setMethod("getData",
          signature = c(x = "NormEnd5SequenceData",
                        bamfiles = "BamFileList",
                        grl = "GRangesList",
                        sequences = "XStringSet",
                        param = "ScanBamParam"),
          definition = function(x, bamfiles, grl, sequences, param, args){
            data <- lapply(bamfiles,
                           FUN = .get_position_data_of_transcript_ends_norm,
                           grl = grl,
                           param = param,
                           type = "5prime",
                           args = args)
            names(data) <- rep("normend5",length(data))
            data
          }
)

#' @rdname NormEndSequenceData-class
#' @export
setMethod("getData",
          signature = c(x = "NormEnd3SequenceData",
                        bamfiles = "BamFileList",
                        grl = "GRangesList",
                        sequences = "XStringSet",
                        param = "ScanBamParam"),
          definition = function(x, bamfiles, grl, sequences, param, args){
            data <- lapply(bamfiles,
                           FUN = .get_position_data_of_transcript_ends_norm,
                           grl = grl,
                           param = param,
                           type = "3prime",
                           args = args)
            names(data) <- rep("normend3",length(data))
            data
          }
)


# aggregation ------------------------------------------------------------------

# aggregate
# - calculate mean per observation
# - calculate sd per observation
#' @importFrom matrixStats rowSds
.aggregate_data_frame_mean_sd <- function(x, condition){
  conditions <- conditions(x)
  f <- .subset_to_condition(conditions, condition)
  df <- as(unlist(x,use.names=FALSE),"DataFrame")[,f,drop=FALSE]
  conditions_u <- unique(conditions[f])
  replicates <- replicates(x)[f]
  # set up some base values
  sample_width <- length(replicates[conditions[f] == conditions_u[1] & 
                                      replicates == unique(replicates)[1]])
  colNames <- strsplit(colnames(df)[seq_len(sample_width)],"\\.")
  colNames <- IRanges::CharacterList(colNames)[as.list(lengths(colNames))]
  # set up some base values. replicates is here the same as the number of 
  # columns, since a list per replicate is assumed
  # get means
  means <- do.call(
    c,
    lapply(conditions_u,
           function(con){
             ff <- conditions[f] == con
             ncol <- ncol(df[,ff,drop = FALSE]
                          [,replicates[ff] == unique(replicates[ff])[1],
                            drop = FALSE])
             seqAdd <- seq.int(from = 0, 
                               to = ncol(df[,ff,drop=FALSE]) - 1, 
                               by = sample_width)
             means <- IRanges::NumericList(
               lapply(seq_len(ncol),
                      function(i){
                        unname(rowMeans(as.data.frame(df[,ff,drop=FALSE][,i + seqAdd,drop=FALSE]),
                                        na.rm = TRUE))
                      }))
             names(means) <- paste0("means.", con, ".", colNames)
             means
           }))
  # get sds
  sds <- do.call(
    c,
    lapply(conditions_u,
           function(con){
             ff <- conditions[f] == con
             ncol <- ncol(df[,ff,drop = FALSE]
                          [,replicates[ff] == unique(replicates[ff])[1],
                            drop = FALSE])
             seqAdd <- seq.int(from = 0, 
                               to = ncol(df[,ff,drop=FALSE]) - 1, 
                               by = sample_width)
             means <- IRanges::NumericList(
               lapply(seq_len(ncol),
                      function(i){
                        unname(matrixStats::rowSds(as.matrix(df[,ff,drop=FALSE][,i + seqAdd,drop=FALSE]),
                                                   na.rm = TRUE))
                      }))
             names(means) <- paste0("sds.", con, ".", colNames)
             means
           }))
  # merge data
  ans <- cbind(do.call(S4Vectors::DataFrame, means),
               do.call(S4Vectors::DataFrame, sds))
  ans <- relist(ans, IRanges::PartitioningByEnd(x))
  positions <- .seqs_rl_strand(ranges(x))
  rownames(ans) <- IRanges::CharacterList(positions)
  ans
}

#' @rdname NormEndSequenceData-class
#' @export
setMethod("aggregateData",
          signature = c(x = "NormEnd5SequenceData"),
          function(x, condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_data_frame_mean_sd(x,condition)
          }
)

#' @rdname NormEndSequenceData-class
#' @export
setMethod("aggregateData",
          signature = c(x = "NormEnd3SequenceData"),
          function(x, condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_data_frame_mean_sd(x,condition)
          }
)

# data visualization -----------------------------------------------------------

RNAMODR_PLOT_SEQ_NORMEND_NAMES <- 
  c("normend5tx" = "mean(5'-ends transcript normalized)",
    "normend3tx" = "mean(3'-ends transcript normalized)",
    "normend5ol" = "mean(5'-ends overlap normalized)",
    "normend3ol" = "mean(3'-ends overlap normalized)")

.clean_mcols_normend <- function(seqdata){
  d <- mcols(seqdata@unlistData)
  d <- d[,stringr::str_detect(colnames(d),"means"),drop=FALSE]
  
  mcols(seqdata@unlistData) <- d
  seqdata
}

#' @rdname NormEndSequenceData-class
#' @export
setMethod(
  f = "getDataTrack",
  signature = signature(x = "NormEnd5SequenceData"),
  definition = function(x, name, ...) {
    args <- list(...)
    # DataTrack for sequence data
    seqdata <- .get_data_for_visualization(x, name)
    # clean meta data columns
    seqdata <- .clean_mcols_normend(seqdata)
    seqdata <- unlist(seqdata)
    conditions_u <- unique(conditions(x))
    if("control" %in% conditions_u){
      d <- seqdata[,stringr::str_detect(colnames(mcols(seqdata)),"control")]
      colnames(mcols(d)) <- gsub(".control","",colnames(mcols(d)))
      dt.control.tx <- Gviz::DataTrack(
        range = d[,"means.tx"],
        group = "means",
        name = paste0(RNAMODR_PLOT_SEQ_NORMEND_NAMES["normend5tx"],
                      "\ncontrol"),
        type = "histogram")
      Gviz::displayPars(dt.control.tx)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.control.tx)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.control.tx)$col.axis <- "#000000"
      Gviz::displayPars(dt.control.tx) <- args
      d <- seqdata[,stringr::str_detect(colnames(mcols(seqdata)),"control")]
      colnames(mcols(d)) <- gsub(".control","",colnames(mcols(d)))
      dt.control.ol <- Gviz::DataTrack(
        range = d[,"means.ol"],
        group = "means",
        name = paste0(RNAMODR_PLOT_SEQ_NORMEND_NAMES["normend5ol"],
                      "\ncontrol"),
        type = "histogram")
      Gviz::displayPars(dt.control.ol)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.control.ol)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.control.ol)$col.axis <- "#000000"
      Gviz::displayPars(dt.control.ol) <- args
      tracks <- list("NormEnd5tx" = dt.control.tx,
                     "NormEnd5ol" = dt.control.ol)
    }
    if("treated" %in% conditions_u){
      d <- seqdata[,stringr::str_detect(colnames(mcols(seqdata)),"treated")]
      colnames(mcols(d)) <- gsub(".treated","",colnames(mcols(d)))
      dt.treated.tx <- Gviz::DataTrack(
        range = d[,"means.tx"],
        group = "means.tx",
        name = paste0(RNAMODR_PLOT_SEQ_NORMEND_NAMES["normend5tx"],
                      "\ntreated"),
        type = "histogram")
      Gviz::displayPars(dt.treated.tx)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.treated.tx)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.treated.tx)$col.axis <- "#000000"
      Gviz::displayPars(dt.treated.tx) <- args
      dt.treated.ol <- Gviz::DataTrack(
        range = d[,"means.ol"],
        group = "means.ol",
        name = paste0(RNAMODR_PLOT_SEQ_NORMEND_NAMES["normend5ol"],
                      "\ntreated"),
        type = "histogram")
      Gviz::displayPars(dt.treated.ol)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.treated.ol)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.treated.ol)$col.axis <- "#000000"
      Gviz::displayPars(dt.treated.ol) <- args
      tracks <- list("NormEnd5tx" = dt.treated.tx,
                     "NormEnd5ol" = dt.treated.ol)
    }
    if(length(conditions_u) == 2L){
      tracks <- list("NormEnd5tx" = dt.control.tx,
                     "NormEnd5ol" = dt.control.ol,
                     "NormEnd5tx" = dt.treated.tx,
                     "NormEnd5ol" = dt.treated.ol)
    }
    tracks
  }
)

#' @rdname NormEndSequenceData-class
#' @export
setMethod(
  f = "getDataTrack",
  signature = signature(x = "NormEnd3SequenceData"),
  definition = function(x, name, ...) {
    args <- list(...)
    # DataTrack for sequence data
    seqdata <- .get_data_for_visualization(x, name)
    # clean meta data columns
    seqdata <- .clean_mcols_normend(seqdata)
    seqdata <- unlist(seqdata)
    conditions_u <- unique(conditions(x))
    if("control" %in% conditions_u){
      d <- seqdata[,stringr::str_detect(colnames(mcols(seqdata)),"control")]
      colnames(mcols(d)) <- gsub(".control","",colnames(mcols(d)))
      dt.control.tx <- Gviz::DataTrack(
        range = d[,"means.tx"],
        group = "means",
        name = paste0(RNAMODR_PLOT_SEQ_NORMEND_NAMES["normend3tx"],
                      "\ncontrol"),
        type = "histogram")
      Gviz::displayPars(dt.control.tx)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.control.tx)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.control.tx)$col.axis <- "#000000"
      Gviz::displayPars(dt.control.tx) <- args
      d <- seqdata[,stringr::str_detect(colnames(mcols(seqdata)),"control")]
      colnames(mcols(d)) <- gsub(".control","",colnames(mcols(d)))
      dt.control.ol <- Gviz::DataTrack(
        range = d[,"means.ol"],
        group = "means",
        name = paste0(RNAMODR_PLOT_SEQ_NORMEND_NAMES["normend3ol"],
                      "\ncontrol"),
        type = "histogram")
      Gviz::displayPars(dt.control.ol)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.control.ol)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.control.ol)$col.axis <- "#000000"
      Gviz::displayPars(dt.control.ol) <- args
      tracks <- list("NormEnd3tx" = dt.control.tx,
                     "NormEnd3ol" = dt.control.ol)
    }
    if("treated" %in% conditions_u){
      d <- seqdata[,stringr::str_detect(colnames(mcols(seqdata)),"treated")]
      colnames(mcols(d)) <- gsub(".treated","",colnames(mcols(d)))
      dt.treated.tx <- Gviz::DataTrack(
        range = d[,"means.tx"],
        group = "means.tx",
        name = paste0(RNAMODR_PLOT_SEQ_NORMEND_NAMES["normend3tx"],
                      "\ntreated"),
        type = "histogram")
      Gviz::displayPars(dt.treated.tx)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.treated.tx)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.treated.tx)$col.axis <- "#000000"
      Gviz::displayPars(dt.treated.tx) <- args
      dt.treated.ol <- Gviz::DataTrack(
        range = d[,"means.ol"],
        group = "means.ol",
        name = paste0(RNAMODR_PLOT_SEQ_NORMEND_NAMES["normend3ol"],
                      "\ntreated"),
        type = "histogram")
      Gviz::displayPars(dt.treated.ol)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.treated.ol)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.treated.ol)$col.axis <- "#000000"
      Gviz::displayPars(dt.treated.ol) <- args
      tracks <- list("NormEnd3tx" = dt.treated.tx,
                     "NormEnd3ol" = dt.treated.ol)
    }
    if(length(conditions_u) == 2L){
      tracks <- list("NormEnd3tx" = dt.control.tx,
                     "NormEnd3ol" = dt.control.ol,
                     "NormEnd3tx" = dt.treated.tx,
                     "NormEnd3ol" = dt.treated.ol)
    }
    tracks
  }
)
