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
#' @import Rsamtools
#' @import GenomicRanges
#' @import ggplot2
NULL

requireNamespace("BiocParallel")
requireNamespace("Biostrings")
requireNamespace("rtracklayer")
requireNamespace("GenomicRanges")
requireNamespace("assertive")


RNAMOD_GFF_FEATURE_TYPES <- c("gene","mRNA","transcript","CDS",
                                "biological_region","five_prime_UTR",
                                "three_prime_UTR")
RNAMOD_GFF_COLNAMES <- c("source","type","score","phase","ID","Name","Parent",
                           "Alias","gene")

RNAMOD_MOD_CONTAINING_FEATURES <- c("transcript",
                                    "mRNA",
                                    "rRNA_gene",
                                    "tRNA_gene",
                                    "ncRNA_gene",
                                    "snoRNA",
                                    "snRNA_gene")


.onLoad <- function(libname,pkgname){
  options("RNAmod_sample_transcripts" = c("RDN18-1"))
  options("RNAmod_dpi" = 600)
  options("RNAmod_use_p" = TRUE)
}