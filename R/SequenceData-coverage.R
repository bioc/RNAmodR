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

#' @rdname CoverageSequenceData
#' @export
CoverageSequenceData <- function(bamfiles,
                                 annotation,
                                 sequences,
                                 seqinfo,
                                 ...){
  args <- .get_mod_data_args(...)
  txdb <- .norm_annotation(annotation)
  sequences <- .norm_sequences(sequences)
  seqinfo <- .norm_seqnames(bamfiles, annotation, sequences, seqinfo)
  ans <- new("CoverageSequenceData",
             bamfiles,
             seqinfo,
             args)
  ranges <- .load_annotation(txdb, ans@seqinfo)
  sequences <- .load_transcript_sequences(sequences, ranges)
  param <- .assemble_scanBamParam(ranges, ans@minQuality, ans@seqinfo)
  message("Loading Coverage data from BAM files...")
  data <- lapply(ans@bamfiles,
                 FUN = .get_position_data_of_transcript_coverage,
                 ranges = ranges,
                 param = param,
                 args = args)
  names(data) <- paste0("coverage.",
                        names(ans@bamfiles),
                        ".",
                        seq_along(ans@bamfiles))
  .postprocess_read_data(ans,
                         data,
                         ranges,
                         sequences)
}
