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

#' @importFrom Rsamtools idxstatsBam
# Extracts sequence names aka. chromosome identifier from list of bam files
.get_acceptable_chrom_ident <- function(bamFiles){
  bf_seqnames <- lapply(bamFiles, function(file){
    res <- Rsamtools::idxstatsBam(file)
    return(as.character(res$seqnames))
  })
  return(unique(unlist(bf_seqnames)))
}
