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

#' #' @title .loadGff
#' #' 
#' #' @description
#' #' Loads saved gff3 as GRanges object from file
#' setGeneric ( 
#'   name= ".loadGff", 
#'   def=function(.Object, 
#'                experiment, 
#'                modification,
#'                failOnNonExist = TRUE){
#'     standardGeneric(".loadGff")
#'   } 
#' ) 

#' @title .saveGff
#' 
#' @description
#' Saves GRanges object as gff3 file
setGeneric ( 
  name= ".saveGff", 
  def=function(.Object, 
               gff, 
               experiment, 
               modification ){standardGeneric(".saveGff")} 
) 