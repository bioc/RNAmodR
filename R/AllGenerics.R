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
               modification ){standardGeneric(".saveSummarizedExperiments")} 
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
               modification ){standardGeneric(".saveGffResult")} 
) 