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
                                    "gene",
                                    "mRNA",
                                    "rRNA",
                                    "rRNA_gene",
                                    "tRNA",
                                    "tRNA_gene",
                                    "ncRNA",
                                    "ncRNA_gene",
                                    "snoRNA",
                                    "snoRNA_gene",
                                    "snRNA",
                                    "snRNA_gene",
                                    "SRP_RNA")
RNAMOD_MOD_SEQ_FEATURES <- c("five_prime_UTR",
                             "three_prime_UTR",
                             "CDS",
                             "noncoding_exon",
                             "exon")
RNAMOD_MOD_DIRECT_SEQ_FEATURES <- c("rRNA",
                                    "rRNA_gene")
# RNAMOD_MOD_SEQ_FEATURES <- c("CDS",
#                              "noncoding_exon",
#                              "exon",
#                              "rRNA",
#                              "rRNA_gene",
#                              "tRNA",
#                              "tRNA_gene",
#                              "ncRNA",
#                              "ncRNA_gene",
#                              "snoRNA",
#                              "snoRNA_gene",
#                              "snRNA",
#                              "snRNA_gene",
#                              "SRP_RNA")

# Settings
RNAMOD_DEFAULT_PALETTE <- "Set1"
RNAMOD_DEFAULT_TRANSCRIPT_MAX_ITERATIONS <- 5

.onLoad <- function(libname,pkgname){
  options("RNAmod_sample_transcripts" = c("RDN18-1"))
  options("RNAmod_transcript_max_iteration" = 
            RNAMOD_DEFAULT_TRANSCRIPT_MAX_ITERATIONS)
  options("RNAmod_use_p" = TRUE)
  options("RNAmod_palette" = RNAMOD_DEFAULT_PALETTE)
  options("RNAmod_dpi" = 600)
  options("RNAmod_use_cairo" = TRUE)
  options("RNAmod_debug" = FALSE)
  options("RNAmod_debug_transcripts" = c("tG(UCC)G"))
}