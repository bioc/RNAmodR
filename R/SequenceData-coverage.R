#' @include RNAmodR.R
#' @include SequenceData-class.R
NULL

#' @name CoverageSequenceData
#' 
#' @title CoverageSequenceData
#' 
#' @description
#' \code{CoverageSequenceData} implements
#' \code{\link[=SequenceData-class]{SequenceData}} to contain and aggregate the
#' coverage of reads per position along the transcripts.
#' 
#' \code{aggregate} calculates the mean and sd for samples in the \code{control}
#' and \code{treated} condition serparatly.
#' 
#' @param bamfiles,annotation,sequences,seqinfo,... See 
#' \code{\link[=SequenceData-class]{SequenceData}}
#' @param x a \code{CoverageSequenceData}
#' @param condition For \code{\link{aggregate}}: condition for which the data 
#' should be aggregated.
#' 
#' @examples
#' # Construct a CoverageSequenceData object
#' annotation <- system.file("extdata","example1.gff3",package = "RNAmodR.Data")
#' sequences <- system.file("extdata","example1.fasta",package = "RNAmodR.Data")
#' files <- c(control = system.file("extdata","example_wt_1.bam",
#'                                  package = "RNAmodR.Data"),
#'            treated = system.file("extdata","example_wt_2.bam",
#'                                  package = "RNAmodR.Data"))
#' csd <- CoverageSequenceData(files, annotation = annotation,
#'                             sequences = sequences)
#' # aggregate data
#' aggregate(csd)
NULL

#' @rdname CoverageSequenceData
#' @export
setClass(Class = "CoverageSequenceData",
         contains = "SequenceData",
         prototype = list(minQuality = 5L))

#' @rdname CoverageSequenceData
#' @export
CoverageSequenceData <- function(bamfiles, annotation, sequences, seqinfo, ...){
  SequenceData("Coverage", bamfiles = bamfiles, annotation = annotation,
               sequences = sequences, seqinfo = seqinfo, ...)
}

# CoverageSequenceData ---------------------------------------------------------
.get_position_data_of_transcript_coverage <- function(bamFile, grl, param,
                                                      args = list()){
  # get data per chromosome
  coverage <- GenomicAlignments::coverage(bamFile, param = param)
  coverage <- as(coverage,"IntegerList")
  # subset per transcript
  seqs <- .seqs_rl(grl)
  coverage <- IRanges::IntegerList(mapply(
    function(gr,s){
      coverage[[GenomicRanges::seqnames(gr)]][s]
    },
    grl,
    seqs,
    SIMPLIFY = FALSE))
  coverage
}

#' @rdname RNAmodR-internals
setMethod(".getData",
          signature = c(x = "CoverageSequenceData",
                        grl = "GRangesList",
                        sequences = "XStringSet",
                        param = "ScanBamParam"),
          definition = function(x, grl, sequences, param, args){
            message("Loading Coverage data from BAM files ... ",
                    appendLF = FALSE)
            files <- bamfiles(x)
            data <- lapply(files,
                           FUN = .get_position_data_of_transcript_coverage,
                           grl = grl,
                           param = param,
                           args = args)
            names(data) <- paste0("coverage.", x@condition, ".", x@replicate)
            data
          })

# aggregation ------------------------------------------------------------------

#' @rdname CoverageSequenceData
#' @export
setMethod("aggregate",
          signature = c(x = "CoverageSequenceData"),
          function(x, condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_list_data_mean_sd(x, condition)
          }
)

# data visualization -----------------------------------------------------------

RNAMODR_PLOT_SEQ_COVERAGE_NAMES <- c("means" = "mean(coverage)")

.clean_mcols_coverage <- function(seqdata){
  d <- mcols(seqdata@unlistData)
  d <- d[,!stringr::str_detect(colnames(d),"sds."),drop=FALSE]
  mcols(seqdata@unlistData) <- d
  seqdata
}

#' @rdname CoverageSequenceData
setMethod(
  f = "getDataTrack",
  signature = signature(x = "CoverageSequenceData"),
  definition = function(x, ...) {
    args <- list(...)
    name <- .norm_viz_name(args[["name"]])
    # DataTrack for sequence data
    seqdata <- .get_data_for_visualization(x, name)
    # clean meta data columns
    seqdata <- .clean_mcols_coverage(seqdata)
    seqdata <- unlist(seqdata)
    conditions <- unique(x@condition)
    if("control" %in% conditions){
      d <- seqdata[,stringr::str_detect(colnames(mcols(seqdata)),"control")]
      colnames(mcols(d)) <- gsub(".control","",colnames(mcols(d)))
      dt.control <- Gviz::DataTrack(range = d,
                                    group = "means",
                                    name = paste0(RNAMODR_PLOT_SEQ_COVERAGE_NAMES["means"],
                                                  "\ncontrol"),
                                    type = "histogram")
      Gviz::displayPars(dt.control)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.control)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.control)$col.axis <- "#000000"
      Gviz::displayPars(dt.control) <- args
      track <- dt.control
    }
    if("treated" %in% conditions){
      d <- seqdata[,stringr::str_detect(colnames(mcols(seqdata)),"treated")]
      colnames(mcols(d)) <- gsub(".treated","",colnames(mcols(d)))
      dt.treated <- Gviz::DataTrack(range = d,
                                    group = "means",
                                    name = paste0(RNAMODR_PLOT_SEQ_COVERAGE_NAMES["means"],
                                                  "\ntreated"),
                                    type = "histogram")
      Gviz::displayPars(dt.treated)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.treated)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.treated)$col.axis <- "#000000"
      Gviz::displayPars(dt.treated) <- args
      track <- dt.treated
    }
    if(length(conditions) == 2L){
      track <- list("1" = dt.control,
                    "1" = dt.treated)
    }
    track
  }
)
