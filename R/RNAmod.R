#' @title 
#' RNAmod
#' 
#' @author Felix G M Ernst [aut]
#' 
#' @description
#' todo
#' 
#'
#' @docType package
#' @name RNAmod
NULL

#' @import methods
#' @import assertive
#' @import BiocParallel
#' @import Biostrings
#' @import rtracklayer
#' @import GenomicRanges
NULL

# Load it directly to allow user interacion with SE object right away
requireNamespace("BiocParallel")
requireNamespace("Biostrings")
requireNamespace("rtracklayer")
requireNamespace("GenomicRanges")
requireNamespace("assertive")
