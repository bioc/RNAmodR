#' @include RNAmodR.R
NULL


# experiment data---------------------------------------------------------------

#' @title .loadSummarizedExperiment
#' 
#' @description
#' Loads saved SummarizedExperiment from file
setGeneric ( 
  name= ".loadSummarizedExperiment", 
  def=function(.Object, 
               experiment, 
               modification,
               failOnNonExist = TRUE){
    standardGeneric(".loadSummarizedExperiment")
  } 
) 

#' @title .saveSummarizedExperiments
#' 
#' @description
#' Saves SummarizedExperiment to file
setGeneric ( 
  name= ".saveSummarizedExperiments", 
  def=function(.Object, 
               se, 
               experiment, 
               modification ) standardGeneric(".saveSummarizedExperiments")
) 

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

#' @name convertReadsToData
setGeneric( 
  name = "convertReadsToData",
  def = function(x,
                 files,
                 conditions,
                 gff,
                 param) standardGeneric("convertReadsToData")
) 

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

#' @name .getDataOfTranscript
#' @title Returns the data for each transcript
setGeneric ( 
  name= ".getDataOfTranscript", 
  def=function(x,
               bamData,
               counts,
               gff) standardGeneric(".getDataOfTranscript")
) 


# mod type functions -----------------------------------------------------------

#' @name checkForModification
#' 
#' @title checking for modifications
#' 
#' @param x a mod class
#' @param location a local position on a transcript
#' @param locations all local positions on the transcript
#' @param data the data to be analyzed
setGeneric( 
  name = "checkForModification",
  def = function(x,
                 location,
                 locations,
                 data) standardGeneric("checkForModification")
) 