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
requireNamespace("Rsamtools")
requireNamespace("IRanges")
requireNamespace("GenomicRanges")
requireNamespace("GenomicAlignments")
requireNamespace("Biostrings")
requireNamespace("Modstrings")
requireNamespace("rtracklayer")
requireNamespace("BiocParallel")


# constants --------------------------------------------------------------------

SAMPLE_TYPES <- c("treated","control")

# constants for annotation -----------------------------------------------------

GFF_COLNAMES <- c("source",
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
MOD_CONTAINING_FEATURES <- c("transcript",
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
MOD_INVALID_FEATURES <- c("chromosome",
                          "telomere",
                          "telomeric_repeat")
MOD_TRANSCRIPT_FEATURES <- c("gene",
                             "mRNA",
                             "five_prime_UTR",
                             "three_prime_UTR",
                             "CDS",
                             "noncoding_exon",
                             "exon",
                             "intron",
                             "five_prime_UTR_intron",
                             "three_prime_UTR_intron")


#' @name RNAmodR-internals
#' @aliases .dataTrack
#' 
#' @title RNAmodR internal functions
#' 
#' @description
#' These functions are not intended for general use, but maybe exported for 
#' additional package development.
NULL