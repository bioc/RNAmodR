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
#' @import BiocParallel
#' @import Biostrings
#' @import rtracklayer
#' @import Rsamtools
#' @import GenomicRanges
#' @import GenomicFeatures
#' @import GenomicAlignments
#' @import ggplot2
#' @importFrom stats setNames wilcox.test sd median
NULL

requireNamespace("BiocParallel")
requireNamespace("Biostrings")
requireNamespace("rtracklayer")
requireNamespace("GenomicRanges")
requireNamespace("GenomicFeatures")
requireNamespace("GenomicAlignments")
requireNamespace("assertive")

# RNAMODR_GFF_FEATURE_TYPES <- c("gene",
#                               "mRNA",
#                               "transcript",
#                               "CDS",
#                               "biological_region",
#                               "five_prime_UTR",
#                               "three_prime_UTR")


# GFF/GRanges contants ---------------------------------------------------------
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

# Import constants -------------------------------------------------------------
RNAMODR_MIN_MAPQUALITY <- 5
RNAMODR_PARAM_COL <- c("Identifier",
                       "Param",
                       "Value")
RNAMODR_SAMPLE_TYPES <- c("Treated",
                          "Control")
RNAMODR_DEFAULT_COLNAMES <- c("ExperimentNo",
                              "SampleName",
                              "ShortName",
                              "Replicate",
                              "Conditions",
                              "BamFile")

# Visualization settings -------------------------------------------------------
RNAMODR_DEFAULT_PALETTE <- "Set1"

# Settings ---------------------------------------------------------------------
.onLoad <- function(libname,pkgname){
  options("RNAmodR_palette" = RNAMODR_DEFAULT_PALETTE)
  options("RNAmodR_dpi" = 600)
  options("RNAmodR_debug" = FALSE)
  options("RNAmodR_min_mapquality" = RNAMODR_MIN_MAPQUALITY)
}

