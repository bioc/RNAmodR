#' @include AllGenerics.R
NULL

# RNAmodR ================================================================
#' @name RNAmodR-Accessors 
#' 
#' @title Accessors for RNAmodR Object
#' 
#' @description
#' Accessors for RNAmodR Object
#'
#' @param .Object a RNAmodR object
#'
#' @return 
#' \code{getInputFolder}, \code{getOutputFolder}: the requested folder path
#' 
#' \code{getExperimentName}: the title of the experiment collection. This is 
#' also the path identifier, which is expected to contain the data and source 
#' folder.
#' 
#' \code{getGFFFile}: the information of the input gff file as GRanges object
#' 
#' \code{getFastaFile}: the link to the fasta file as FaFile class
#' 
#' @examples
#' \donttest{
#' unzip(system.file("extdata", package = "RNAmodR", file = "RNAmodR.zip" ))
#' mod <- RNAmodR("test",
#'                "test_layout.csv",
#'                "test_gff.gff3",
#'                "test_masked.fasta")
#' getInputFolder(mod)  
#' getOutputFolder(mod)  
#' getExperimentName(mod)  
#' getGFFFile(mod)  
#' getFastaFile(mod)
#' }
NULL

#' @rdname RNAmodR-Accessors
#' @export
setGeneric( 
  name = "getInputFolder",
  def = function(.Object ){standardGeneric("getInputFolder")} 
) 
#' @rdname RNAmodR-Accessors
#' @export
setGeneric( 
  name = "getOutputFolder",
  def = function(.Object ){standardGeneric("getOutputFolder")} 
) 
#' @rdname RNAmodR-Accessors
#' @export
setGeneric( 
  name = "getExperimentName",
  def = function(.Object ){standardGeneric("getExperimentName")} 
) 
#' @rdname RNAmodR-Accessors
#' @export
setGeneric( 
  name = "getGFFFile",
  def = function(.Object ){standardGeneric("getGFFFile")} 
)
#' @rdname RNAmodR-Accessors
#' @export
setGeneric( 
  name = "getFastaFile",
  def = function(.Object ){standardGeneric("getFastaFile")} 
) 

# setup helper functions -------------------------------------------------------

#' @title setupWorkEnvir
#' 
#' @description
#' conveniance function for setting up folders for an experiment
#' 
#' @export
setGeneric( 
  name = "setupWorkEnvir", 
  def = function(experimentName){standardGeneric("setupWorkEnvir")} 
) 


# experiment data---------------------------------------------------------------

#' @name getExperimentData
#' @export
setGeneric( 
  name = "getExperimentData",
  def = function(.Object,
                 number){standardGeneric("getExperimentData")} 
) 


#' @name getSummarizedExperiment
#' @export
setGeneric( 
  name = "getSummarizedExperiment",
  def = function(.Object,
                 number, 
                 modification){standardGeneric("getSummarizedExperiment")} 
) 

#' @rdname  getSummarizedExperiment
#' @export
setGeneric( 
  name = "getSummarizedExperiments",
  def = function(.Object,
                 number, 
                 modification){standardGeneric("getSummarizedExperiments")} 
) 

#' @name setSummarizedExperiment
#' @export
setGeneric( 
  name = "setSummarizedExperiment",
  def = function(.Object,
                 se,
                 number, 
                 modification ){standardGeneric("setSummarizedExperiment")} 
) 

#' @name getGffResult
#' @export
setGeneric( 
  name = "getGffResult",
  def = function(.Object,
                 number, 
                 modification,
                 genomicCoordinates = FALSE){standardGeneric("getGffResult")} 
) 

#' @name setGffResult
#' 
#' @description setGff
#' 
#' @export
setGeneric( 
  name = "setGffResult",
  def = function(.Object,
                 gff,
                 number, 
                 modification){standardGeneric("setGffResult")} 
) 

# parsing ----------------------------------------------------------------------

#' @name parseForModifications
#' @export
setGeneric( 
  name = "parseForModifications",
  def = function(.Object,
                 number,
                 name,
                 gff,
                 fasta,
                 files,
                 mapQuality,
                 modifications){standardGeneric("parseForModifications")} 
) 

# analysis type accessors ------------------------------------------------------------

#' @name analysis-accessors
#' @export
setGeneric( 
  name = "getPlotType",
  def = function(object){standardGeneric("getPlotType")} 
) 

#' @rdname analysis-accessors
#' @export
setGeneric( 
  name = "getPositions",
  def = function(object){standardGeneric("getPositions")} 
) 

#' @rdname analysis-accessors
#' @export
setGeneric( 
  name = "getModifications",
  def = function(object){standardGeneric("getModifications")} 
) 

#' @rdname analysis-accessors
#' @export
setGeneric( 
  name = "getDataLabel",
  def = function(object){standardGeneric("getDataLabel")} 
) 

#' @rdname analysis-accessors
#' @export
setGeneric( 
  name = "getDataFormat",
  def = function(object){standardGeneric("getDataFormat")} 
) 

# mod type accessors -----------------------------------------------------------

#' @name mod-accessors
#' @aliases getAnalysisType getModType
#' 
#' @title accessors for \code{mod} class objects
#' 
#' @description 
#' title 
#' 
#' @export
setGeneric( 
  name = "getModType",
  def = function(object){standardGeneric("getModType")} 
) 
#' @rdname mod-accessors
#' @export
setGeneric( 
  name = "getAnalysisType",
  def = function(object){standardGeneric("getAnalysisType")} 
) 
#' @rdname mod-accessors
#' @export
setGeneric( 
  name = "getPositionOffset",
  def = function(object){standardGeneric("getPositionOffset")} 
) 


# analysis type functions - modifications parsing ------------------------------

#' @name convertReadsToPositions
#' 
#' @title convertReadsToPositions
#' 
#' @param object analysis object 
#' @param files files which to analyze
#' @param param ScanBamParam to use for loading the bam files
#' @param gff a GRanges object for the genome
#'
#' @export
setGeneric( 
  name = "convertReadsToPositions",
  def = function(object,
                 files,
                 gff,
                 param){standardGeneric("convertReadsToPositions")} 
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
                 modClasses){standardGeneric("parseMod")} 
) 

#' @name mergePositionsOfReplicates
#' 
#' @title mergePositionsOfReplicates
#' 
#' @param object a mod object 
#' 
#' @export
setGeneric( 
  name = "mergePositionsOfReplicates",
  def = function(object){standardGeneric("mergePositionsOfReplicates")} 
) 


# mod type functions -----------------------------------------------------------

#' @name preTest
#' 
#' @title Testing data for compatibility with modification detection
#' 
#' 
#' @param object a mod class
#' @param location a local position on a transcript
#' @param locations all local positions on the transcript
#' @param data the data to be analyzed
#' 
#' @export
setGeneric( 
  name = "preTest",
  def = function(object,
                 location,
                 locations,
                 data){standardGeneric("preTest")} 
) 

#' @name checkForModification
#' 
#' @title checking for modifications
#' 
#' @param object a mod class
#' @param location a local position on a transcript
#' @param locations all local positions on the transcript
#' @param data the data to be analyzed
#' 
#' @export
setGeneric( 
  name = "checkForModification",
  def = function(object,
                 location,
                 locations,
                 data){standardGeneric("checkForModification")} 
) 


# modification visualization ---------------------------------------------------

#' @name getModPlot
#' @export
setGeneric( 
  name = "getModPlot",
  def = function(x,
                 number,
                 modifications,
                 gene,
                 gff,
                 fasta,
                 focus = FALSE){standardGeneric("getModPlot")} 
)

#' @rdname getModPlot
#' @export
setGeneric( 
  name = "saveModPlot",
  def = function(x,
                 number,
                 modifications,
                 genes,
                 gff,
                 fasta,
                 focus = FALSE,
                 folder,
                 filetype = "pdf"){standardGeneric("saveModPlot")} 
)

#' @name heatmapModifications
#' @export
setGeneric(
  name = "heatmapModifications",
  def = function(x,
                 numbers,
                 sampleNames,
                 modifications,
                 genes){standardGeneric("heatmapModifications")}
)
#' @rdname heatmapModifications
#' @export
setGeneric(
  name = "saveHeatmapModifications",
  def = function(x,
                 numbers,
                 sampleNames,
                 modifications,
                 genes,
                 folder){standardGeneric("saveHeatmapModifications")}
)

# comparison -------------------------------------------------------------------

#' @name compareSample
#'
#' @title compareSample
#' 
#' @param .Object 
#' @param number 
#' @param modifications 
#' 
#' @export
setGeneric( 
  name = "compareSample",
  def = function(.Object,
                 number,
                 modifications){standardGeneric("compareSample")} 
)

#' @name compareAsHeatmap
#'
#' @title compareAsHeatmap
#' 
#' @param .Object 
#' @param number 
#' @param modifications 
#' @param genes 
#' 
#' @export
setGeneric( 
  name = "compareAsHeatmap",
  def = function(.Object,
                 number,
                 modifications,
                 genes){standardGeneric("compareAsHeatmap")} 
)


# SummarizedExperiment ---------------------------------------------------------


#' @name seGetModifications
#' @export
setGeneric( 
  name = "seGetModifications",
  def = function(se,
                 genes,
                 modifications){standardGeneric("seGetModifications")} 
)
#' @rdname seGetModifications
#' @export
setGeneric( 
  name = "seGetPositions",
  def = function(se,
                 genes){standardGeneric("seGetPositions")} 
)