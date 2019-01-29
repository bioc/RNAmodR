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
.get_position_data_of_transcript_coverage <- function(bamFile,
                                                      ranges,
                                                      param,
                                                      args = list()){
  ranges <- .get_parent_annotations(ranges)
  # get data per chromosome
  coverage <- GenomicAlignments::coverage(bamFile, param = param)
  coverage <- 
    coverage[names(coverage) %in% as.character(seqnames(ranges))]
  coverage <- as(coverage,"IntegerList")
  coverage <- IRanges::IntegerList(
    lapply(split(ranges,
                 seq_along(ranges)),
           function(r){
             coverage[[seqnames(r)]][start(r):end(r)]
           }))
  names(coverage) <- ranges$ID
  coverage
}

setMethod(".get_Data",
          signature = c(x = "CoverageSequenceData",
                        grl = "GRangesList",
                        sequences = "XStringSet",
                        param = "ScanBamParam"),
          definition = function(x,
                                grl,
                                sequences,
                                param){
            data <- lapply(x@bamfiles,
                           FUN = .get_position_data_of_transcript_coverage,
                           ranges = grl,
                           param = param,
                           args = settings(x))
            names(data) <- paste0("coverage.",
                                  names(x@bamfiles),
                                  ".",
                                  seq_along(x@bamfiles))
            data
          })

#' @rdname CoverageSequenceData
#' @export
CoverageSequenceData <- function(bamfiles, annotation, sequences, seqinfo, ...){
  browser()
  SequenceData("Coverage", files = bamfiles, annotation = annotation,
               sequences = sequences, seqinfo = seqinfo, ...)
}
