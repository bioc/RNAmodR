#' @include RNAmod.R
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