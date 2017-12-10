#' @include AllGenerics.R
NULL


# RNAmod ================================================================
#' @name RNAmod-Accessors 
#' 
#' @title Accessors for RNAmod Object
#' 
#' @description
#' Accessors for RNAmod Object
#'
#' @param .Object a RNAmod object
#'
#' @return 
#' *Folder: the requested folder path
#' 
#' *Name: the title of the experiment collection. This is also the path 
#' identifier, which is expected to contain the data and source folder.
NULL

#' @rdname RNAmod-Accessors
#' 
#' @export
setGeneric ( 
  name= "getInputFolder",
  def=function(.Object ){standardGeneric("getInputFolder")} 
) 
#' @rdname RNAmod-Accessors
#' 
#' @export
setGeneric ( 
  name= "getOutputFolder",
  def=function(.Object ){standardGeneric("getOutputFolder")} 
) 
#' @rdname RNAmod-Accessors
#' 
#' @export
setGeneric ( 
  name= "getExperimentName",
  def=function(.Object ){standardGeneric("getExperimentName")} 
) 


#' @title setupWorkEnvir
#' 
#' @description
#' conveniance function for setting up folders for an experiment
#' 
#' @export
setGeneric ( 
  name= "setupWorkEnvir", 
  def=function(experimentName){standardGeneric("setupWorkEnvir")} 
) 


# experiment data---------------------------------------------------------------

#' @name getExperimentData
#' 
#' @param .Object 
#' @param number 
#' 
#' @export
setGeneric ( 
  name= "getExperimentData",
  def=function(.Object,
               number){standardGeneric("getExperimentData")} 
) 

# parsing ----------------------------------------------------------------------

setGeneric ( 
  name= "parse",
  def=function(.Object,
               number){standardGeneric("parse")} 
) 
