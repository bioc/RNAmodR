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
  # get data per chromosome
  coverage <- coverage(bamFile, param = param)
  coverage <- coverage[names(coverage) %in% as.character(seqnames(ranges))]
  coverage <- as(coverage,"IntegerList")
  coverage <- coverage[order(names(coverage))]
  # split into per gene data
  ranges <- .get_parent_annotations(ranges)
  ranges <- split(ranges,GenomeInfoDb::seqnames(ranges))
  ranges <- ranges[order(names(ranges))]
  if(any(names(ranges) != names(coverage))){
    stop("Something went wrong.")
  }
  coverage <- mapply(
    function(co,r){
      # can this be streamlined ?
      starts <- start(r)
      ends <- end(r)
      ans <- lapply(seq_along(r),
             function(i){
               co[starts[i]:ends[i]]
             })
      names(ans) <- r$ID
      ans
    },
    coverage,
    ranges,
    SIMPLIFY = FALSE)
  as(do.call(c,unname(coverage)),"IntegerList")
}

#' @rdname CoverageSequenceData
#' @export
CoverageSequenceData <- function(bamfiles,
                            fasta,
                            gff,
                            ...){
  ans <- new("CoverageSequenceData",
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
