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
#' \code{getInputFolder}, \code{getOutputFolder}: the requested folder path
#' 
#' \code{getExperimentName}: the title of the experiment collection. This is 
#' also the path identifier, which is expected to contain the data and source 
#' folder.
#' 
#' \code{getGFF}: the information of the input gff file as GRanges object
#' 
#' \code{getFasta}: the link to the fasta file as FaFile class
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
#' @rdname RNAmod-Accessors
#' 
#' @export
setGeneric ( 
  name= "getGFF",
  def=function(.Object ){standardGeneric("getGFF")} 
) 
#' @rdname RNAmod-Accessors
#' 
#' @export
setGeneric ( 
  name= "getFasta",
  def=function(.Object ){standardGeneric("getFasta")} 
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

#' @title matchFastaToGff
#' 
#' @description
#' title
#' 
#' @export
setGeneric ( 
  name= "matchFastaToGff", 
  def=function(inputFasta,
               inputGFF){standardGeneric("matchFastaToGff")} 
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


#' @name getSummarizedExperiment
#' 
#' @title returns one or more SummarizedExperiment
#' 
#' @description
#' Global access to SummarizedExperiments stored by RpfExperiment. 
#' \code{getSummarizedExperiment()} returns the result of experiment, whereas 
#' \code{getSummarizedExperiments()}, is the vectorized version returning a
#' list of experiment results.
#' 
#' @export
setGeneric ( 
  name= "getSummarizedExperiment",
  def=function(.Object,
               number, 
               modification){standardGeneric("getSummarizedExperiment")} 
) 

#' @rdname  getSummarizedExperiment
#' 
#' @export
setGeneric ( 
  name= "getSummarizedExperiments",
  def=function(.Object,
               number, 
               modification){standardGeneric("getSummarizedExperiments")} 
) 

#' @name setSummarizedExperiment
#' 
#' @title sets a SummarizedExperiment object
#' 
#' @description
#' Saves/overwrites a SummarizedExperiment object for certain experiment and 
#' passes the SummarizedExperiment object on to be saved to disk as .RData file
#' in the \code{results\\SE} folder
#' 
#' @export
setGeneric ( 
  name= "setSummarizedExperiment",
  def=function(.Object,
               se,
               number, 
               modification ){standardGeneric("setSummarizedExperiment")} 
) 

#' @name getGff
#' 
#' @description getGff
#' 
#' @export
setGeneric ( 
  name= "getGff",
  def=function(.Object,
               number, 
               modification){standardGeneric("getGff")} 
) 

#' @name setGff
#' 
#' @description setGff
#' 
#' @export
setGeneric ( 
  name= "setGff",
  def=function(.Object,
               gff,
               number, 
               modification){standardGeneric("setGff")} 
) 

# parsing ----------------------------------------------------------------------

#' @name parseForModifications
#'
#' @export
setGeneric ( 
  name= "parseForModifications",
  def=function(.Object,
               number,
               modifications){standardGeneric("parseForModifications")} 
) 

# modtype accessors ------------------------------------------------------------

#' @name mod-accessors
#' 
#' @title Accessor for 'mod' class
#' 
#' @param object mod object 
#' 
#' @export
setGeneric ( 
  name= "getPlotType",
  def=function(object){standardGeneric("getPlotType")} 
) 

#' @rdname mod-accessors
#' 
#' @export
setGeneric ( 
  name= "getModType",
  def=function(object){standardGeneric("getModType")} 
) 


# modification parsing ---------------------------------------------------------

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
               data){standardGeneric("convertReadsToPositions")} 
) 

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
               gff,
               seq,
               data){standardGeneric("parseMod")} 
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

#' @name getModPlot
#'
#' @title getModPlot
#' 
#' @export
setGeneric ( 
  name= "getModPlot",
  def=function(.Object,
               se,
               modifications,
               gene,
               focus = FALSE){standardGeneric("getModPlot")} 
)

#' @rdname getModPlot
#' 
#' @export
setGeneric ( 
  name= "saveModPlot",
  def=function(.Object,
               se,
               modifications,
               genes,
               focus = FALSE,
               filetype = "pdf"){standardGeneric("saveModPlot")} 
)

#' @name heatmapModifications
#'
#' @export
setGeneric (
  name= "heatmapModifications",
  def=function(.Object,
               se,
               modifications,
               genes){standardGeneric("heatmapModifications")}
)

#' @rdname heatmapModifications
#'
#' @export
setGeneric (
  name= "saveHeatmapModifications",
  def=function(.Object,
               se,
               modifications,
               genes){standardGeneric("saveHeatmapModifications")}
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

