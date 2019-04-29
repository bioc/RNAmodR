#' @include RNAmodR.R
#' @include SequenceData-class.R
NULL

#' @name CoverageSequenceData-class
#' @aliases CoverageSequenceData
#' 
#' @title CoverageSequenceData
#' 
#' @description
#' \code{CoverageSequenceData} implements
#' \code{\link[=SequenceData-class]{SequenceData}} to contain and aggregate the
#' coverage of reads per position along the transcripts.
#' 
#' \code{CoverageSequenceData} contains one column per data file named using the
#' following naming convention \code{coverage.condition.replicate}.
#' 
#' \code{aggregate} calculates the mean and sd for samples in the \code{control}
#' and \code{treated} condition separatly.
#' 
#' @param bamfiles,annotation,seqinfo,grl,sequences,param,args,... See 
#' \code{\link[=SequenceData-class]{SequenceData}}
#' @param x a \code{CoverageSequenceData}
#' @param name For \code{getDataTrack}: a valid transcript name. Must be a name
#' of \code{ranges(x)}
#' @param condition For \code{\link{aggregate}}: condition for which the data 
#' should be aggregated.
#' 
#' @return a \code{CoverageSequenceData} object
#' 
#' @examples
#' # Construction of a CoverageSequenceData objectobject
#' library(RNAmodR.Data)
#' library(rtracklayer)
#' annotation <- GFF3File(RNAmodR.Data.example.man.gff3())
#' sequences <- RNAmodR.Data.example.man.fasta()
#' files <- c(treated = RNAmodR.Data.example.wt.1())
#' csd <- CoverageSequenceData(files, annotation = annotation,
#'                             sequences = sequences)
NULL

#' @rdname CoverageSequenceData-class
#' @export
setClass(Class = "CoverageSequenceData",
         contains = "SequenceData",
         prototype = list(minQuality = 5L,
                          dataDescription = "Coverage data"))

#' @rdname CoverageSequenceData-class
#' @export
CoverageSequenceData <- function(bamfiles, annotation, sequences, seqinfo, ...){
  SequenceData("Coverage", bamfiles = bamfiles, annotation = annotation,
               sequences = sequences, seqinfo = seqinfo, ...)
}

# CoverageSequenceData ---------------------------------------------------------
#' @importFrom GenomicAlignments coverage
.get_position_data_of_transcript_coverage <- function(bamFile, grl, param,
                                                      args = list()){
  # get data per chromosome
  coverage <- GenomicAlignments::coverage(bamFile, param = param)
  coverage <- as(coverage,"IntegerList")
  # subset per transcript
  seqs <- .seqs_rl_strand(grl, force_continous = TRUE)
  coverage <- IRanges::IntegerList(mapply(
    function(gr,s){
      coverage[[unique(GenomicRanges::seqnames(gr))]][s]
    },
    grl,
    seqs,
    SIMPLIFY = FALSE))
  coverage
}

#' @rdname CoverageSequenceData-class
#' @export
setMethod("getData",
          signature = c(x = "CoverageSequenceData",
                        bamfiles = "BamFileList",
                        grl = "GRangesList",
                        sequences = "XStringSet",
                        param = "ScanBamParam"),
          definition = function(x, bamfiles, grl, sequences, param, args){
            data <- lapply(bamfiles,
                           FUN = .get_position_data_of_transcript_coverage,
                           grl = grl,
                           param = param,
                           args = args)
            names(data) <- rep("coverage",length(data))
            data
          })

# aggregation ------------------------------------------------------------------

#' @rdname CoverageSequenceData-class
#' @export
setMethod("aggregateData",
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

#' @rdname CoverageSequenceData-class
setMethod(
  f = "getDataTrack",
  signature = signature(x = "CoverageSequenceData"),
  definition = function(x, name, ...) {
    args <- list(...)
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
      track <- list("Coverage" = dt.control)
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
      track <- list("Coverage" = dt.treated)
    }
    if(length(conditions) == 2L){
      track <- list("Coverage" = dt.control,
                    "Coverage" = dt.treated)
    }
    track
  }
)
