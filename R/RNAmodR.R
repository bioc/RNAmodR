#' @title RNAmodR
#' 
#' @author Felix G M Ernst [aut]
#' 
#' @description
#' \code{RNAmodR} implements class and a workflows to detect
#' post-transcriptional RNA modifications in high throughput sequencing data.
#' 
#' From the a basic \code{SequenceData} speficic subclasses are derived
#' containing specific aspects of aligned reads, e.g. 5'-end positions or pileup
#' data. From this \code{Modifier} class can be used to detect specific patterns
#' for individual types of modifications. The \code{SequenceData} classes can be
#' shared by different \code{Modifier} classes allowing easy adapattion to new
#' methods.
#' 
#' @seealso The \code{RNAmodR.RiboMethSeq} and \code{RNAmodR.AlkAnilineSeq}
#' package.
#'
#' @docType package
#' @name RNAmodR
NULL

#' @import methods
#' @import assertive
#' @import S4Vectors
#' @import IRanges
#' @import BiocGenerics
#' @import BiocParallel
#' @import GenomicRanges
#' @import GenomicFeatures
#' @import GenomicAlignments
#' @import Biostrings
#' @import Modstrings
#' @import rtracklayer
#' @import Rsamtools
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
#' @aliases .getData
#' 
#' @title RNAmodR internal functions
#' 
#' @description
#' These functions are not intended for general use, but are used for 
#' additional package development.
#' 
#' @param object,range,data,modType,scoreFun,source,type,x,grl,sequences,param,args,i,exact,value
#' Internal arguments
NULL
