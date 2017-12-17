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
#' @export
setGeneric ( 
  name= "getExperimentData",
  def=function(.Object,
               number){standardGeneric("getExperimentData")} 
) 

# parsing ----------------------------------------------------------------------

#' @name analyzeModifications
#'
#' @export
setGeneric ( 
  name= "analyzeModifications",
  def=function(.Object,
               number,
               modifications){standardGeneric("analyzeModifications")} 
) 


# modification parsing ---------------------------------------------------------


#' @name parseMod
#' 
#' @title parseMod
#' 
#' @param object mod object 
#' @param counts total counts in the BAM file containing the data
#' @param gff a GRanges object for a single gene
#' @param seq a DNAString object for a single gene
#' @param data the read data as DataFrame for a single gene
#' 
#' @export
setGeneric ( 
  name= "parseMod",
  def=function(object,
               counts,
               gff,
               seq,
               data){standardGeneric("parseMod")} 
) 


#' @name mergeModsOfReplicates
#' 
#' @title mergeModsOfReplicates
#' 
#' @param object mod object 
#' @param gff a GRanges object
#' @param seq a DNAStringSet object
#' @param res the results for individual replicates
#' 
#' @export
setGeneric ( 
  name= "mergeModsOfReplicates",
  def=function(object,
               gff,
               seq,
               data){standardGeneric("mergeModsOfReplicates")} 
) 


#' @name convertReadsToPositions
#' 
#' @title convertReadsToPositions
#' 
#' @param object mod object 
#' @param counts total counts in the BAM file containing the data
#' @param gff a GRanges object for a single gene
#' @param data the read data as DataFrame for a single gene
#' 
#' @export
setGeneric ( 
  name= "convertReadsToPositions",
  def=function(object,
               counts,
               gff,
               seq,
               data){standardGeneric("convertReadsToPositions")} 
) 

#' @name mergePositionsOfReplicates
#' 
#' @title mergePositionsOfReplicates
#' 
#' @param object mod object 
#' @param gff a GRanges object for a single gene
#' @param seq a DNAString object for a single gene
#' @param data the read data as DataFrame for a single gene
#' 
#' @export
setGeneric ( 
  name= "mergePositionsOfReplicates",
  def=function(object,
               gff,
               seq,
               data){standardGeneric("mergePositionsOfReplicates")} 
) 


# modification visualization ---------------------------------------------------


#' @name visualizeModifications
#'
#' @title visualizeModifications
#' 
#' @export
setGeneric ( 
  name= "visualizeModifications",
  def=function(.Object,
               number,
               modifications,
               genes){standardGeneric("visualizeModifications")} 
)


# comparison -------------------------------------------------------------------

#' @name compareSample
#'
#' @title compareSample
#' 
#' @export
setGeneric ( 
  name= "compareSample",
  def=function(.Object,
               number,
               modifications){standardGeneric("compareSample")} 
)

#' @name compareAsHeatmap
#'
#' @title compareAsHeatmap
#' 
#' @export
setGeneric ( 
  name= "compareAsHeatmap",
  def=function(.Object,
               number,
               modifications,
               genes){standardGeneric("compareAsHeatmap")} 
)

