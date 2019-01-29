#' @title RNAmodR
#' 
#' @author Felix G M Ernst [aut]
#' 
#' @description
#' todo xczCZXcZZcc cCX 
#'
#' @docType package
#' @name RNAmodR
NULL

#' @import methods
#' @import assertive
#' @import S4Vectors
#' @import IRanges
#' @import rtracklayer
#' @import GenomicRanges
#' @import Rsamtools
#' @import GenomicAlignments
#' @import Biostrings
#' @import Modstrings
#' @import BiocParallel
NULL
requireNamespace("S4Vectors")
requireNamespace("IRanges")
requireNamespace("GenomicRanges")
requireNamespace("GenomicAlignments")
requireNamespace("Modstrings")
requireNamespace("Rsamtools")
requireNamespace("rtracklayer")
requireNamespace("BiocParallel")

# constants --------------------------------------------------------------------

SAMPLE_TYPES <- c("treated","control")

# constants for annotation -----------------------------------------------------

#' @name RNAmodR-internals
#' @aliases .dataTrack
#' 
#' @title RNAmodR internal functions
#' 
#' @description
#' These functions are not intended for general use, but maybe exported for 
#' additional package development.
NULL
