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

RNAMOD_MOD_CONTAINING_FEATURES <- c("transcript","mRNA","rRNA_gene","tRNA_gene",
                                    "internal_transcribed_spacer_region",
                                    "external_transcribed_spacer_region",
                                    "ncRNA_gene","snoRNA","snRNA_gene")