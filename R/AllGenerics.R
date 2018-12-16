#' @include RNAmodR.R
NULL


# experiment data---------------------------------------------------------------

#' @title .saveGffResult
#' 
#' @description
#' Saves GRanges object as gff3 file
setGeneric ( 
  name= ".saveGffResult", 
  def=function(.Object, 
               gff, 
               experiment, 
               modification ) standardGeneric(".saveGffResult")
) 

# analysis type functions - modifications parsing ------------------------------

#' @name parseMod
#' 
#' @title parseMod
#' 
#' @param object mod object 
#' @param gff a GRanges object for the genome
#' @param fafile a FaFile object for the genome
#' @param modClasses list of modification classes used for parsing
#' 
#' @export
setGeneric( 
  name = "parseMod",
  def = function(object,
                 gff,
                 fafile,
                 modClasses) standardGeneric("parseMod")
) 

#' @name mergeDataOfReplicates
setGeneric( 
  name = "mergeDataOfReplicates",
  def = function(x) standardGeneric("mergeDataOfReplicates")
) 


