#' @include RNAmodR.R
#' @include SequenceData-class.R
NULL

#' @name CoverageSequenceData
#' 
#' @title CoverageSequenceData
#' 
#' @description
#' title
NULL

#' @rdname CoverageSequenceData
#' @export
setClass(Class = "CoverageSequenceData",
         contains = "SequenceData",
         prototype = list(minQuality = 5L))

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

setMethod(".get_Data",
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
            names(data) <- paste0("coverage.",
                                  names(files),
                                  ".",
                                  seq_along(files))
            data
          })

#' @rdname CoverageSequenceData
#' @export
CoverageSequenceData <- function(bamfiles, annotation, sequences, seqinfo, ...){
  SequenceData("Coverage", bamfiles = bamfiles, annotation = annotation,
               sequences = sequences, seqinfo = seqinfo, ...)
}

# aggregation ------------------------------------------------------------------

#' @name CoverageSequenceData
#' @export
setMethod("aggregate",
          signature = c(x = "CoverageSequenceData"),
          function(x, condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_list_data_mean_sd(x, condition)
          }
)

# data visualization -----------------------------------------------------------
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
