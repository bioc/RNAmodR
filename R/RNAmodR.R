#' @title RNAmodR
#' 
#' @author Felix G M Ernst [aut]
#' 
#' @description
#' todo
#'
#' @docType package
#' @name RNAmodR
NULL

#' @import methods
#' @import assertive
#' @import BiocParallel
#' @import Biostrings
#' @import rtracklayer
#' @import Rsamtools
#' @import GenomicRanges
#' @import GenomicAlignments
#' @import ggplot2
NULL

requireNamespace("BiocParallel")
requireNamespace("Biostrings")
requireNamespace("rtracklayer")
requireNamespace("GenomicRanges")
requireNamespace("GenomicAlignments")
requireNamespace("assertive")


# RNAMODR_GFF_FEATURE_TYPES <- c("gene",
#                               "mRNA",
#                               "transcript",
#                               "CDS",
#                               "biological_region",
#                               "five_prime_UTR",
#                               "three_prime_UTR")
RNAMODR_GFF_COLNAMES <- c("source",
                         "type",
                         "score",
                         "phase",
                         "ID",
                         "Name",
                         "Parent",
                         "Alias",
                         "gene",
                         "protein_id",
                         "transcript_id")
RNAMODR_MOD_CONTAINING_FEATURES <- c("transcript",
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
RNAMODR_MOD_INVALID_FEATURES <- c("chromosome",
                                 "telomere",
                                 "telomeric_repeat")
RNAMODR_MOD_TRANSCRIPT_FEATURES <- c("gene",
                                    "mRNA",
                                    "five_prime_UTR",
                                    "three_prime_UTR",
                                    "CDS",
                                    "noncoding_exon",
                                    "exon",
                                    "intron",
                                    "five_prime_UTR_intron",
                                    "three_prime_UTR_intron")


# Settings
RNAMODR_DEFAULT_PALETTE <- "Set1"
RNAMODR_DEFAULT_TRANSCRIPT_MAX_ITERATIONS <- 5
RNAMODR_DEFAULT_MAPQ <- 5

.onLoad <- function(libname,pkgname){
  options("RNAmodR_map_quality" = RNAMODR_DEFAULT_MAPQ)
  options("RNAmodR_sample_transcripts" = c("RDN18-1"))
  options("RNAmodR_transcript_max_iteration" = 
            RNAMODR_DEFAULT_TRANSCRIPT_MAX_ITERATIONS)
  options("RNAmodR_use_p" = TRUE)
  options("RNAmodR_palette" = RNAMODR_DEFAULT_PALETTE)
  options("RNAmodR_dpi" = 600)
  options("RNAmodR_use_cairo" = FALSE)
  options("RNAmodR_debug" = FALSE)
  options("RNAmodR_debug_transcripts" = c("tG(UCC)G"))
}