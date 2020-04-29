#' @include RNAmodR.R
NULL

# reading BAM input ------------------------------------------------------------

# return parameters to be used by scanBam
# ranges a GRanges object containing the ranges for search in the BAM file
# quality quality argument used for scanBamParam
#' @importFrom GenomeInfoDb seqnames
#' @importFrom Rsamtools ScanBamParam
.assemble_scanBamParam <- function(grl,
                                   quality,
                                   seqinfo,
                                   what = character(0)){
  if(!is(grl,"GRangesList")){
    stop("GRangesList expected.")
  }
  if(is.null(what)){
    what <- character(0)
  }
  # make sure the GRangesList is split by seqlevels
  grl@unlistData <- unname(grl@unlistData)
  gr <- unlist(grl)
  grl <- split(gr,
               GenomicRanges::seqnames(gr))
  # assemble param
  param <- Rsamtools::ScanBamParam(which = grl, what = what, mapqFilter = quality)
  return(param)
}

#' @importFrom GenomicAlignments readGAlignments
.load_bam_alignment_data <- function(bamFile, param, args){
  data <- GenomicAlignments::readGAlignments(bamFile, param = param)
  if(length(data) == 0L){
    stop("No reads found in data.", call. = FALSE)
  }
  # apply length cut off if set
  if(!is.na(args[["maxLength"]])){
    data <- data[GenomicAlignments::qwidth(data) <= args[["maxLength"]],]
  }
  if(!is.na(args[["minLength"]])){
    data <- data[GenomicAlignments::qwidth(data) >= args[["minLength"]],]
  }
  if(length(data) == 0L){
    stop("No reads found in data with read length equal or between 'minLength'",
         " (",args[["minLength"]]," nt) and 'maxLength' (",args[["maxLength"]],
         " nt).", call. = FALSE)
  }
  data
}
