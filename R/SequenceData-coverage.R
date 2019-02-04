#' @include RNAmodR.R
#' @include SequenceData-class.R
NULL

#' @name CoverageSequenceData
#' 
#' @title CoverageSequenceData
#' 
#' @description
#' \code{CoverageSequenceData} implements \code{\link{SequenceData}} to contain
#' and aggregate the coverage of reads per position along the transcripts.
#' 
#' \code{aggregate} calculates the mean and sd for samples in the \code{control}
#' and \code{treated} condition serparatly.
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

#' @rdname RNAmodR-internals
setMethod(
  f = ".dataTracks",
  signature = signature(x = "CoverageSequenceData",
                        data = "missing",
                        seqdata = "GRanges",
                        sequence = "XString"),
  definition = function(x, seqdata, sequence,  args) {
    requireNamespace("Gviz")
    browser()
  }
)
