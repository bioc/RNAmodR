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
  def = function(.Object ) standardGeneric("getInputFolder")
) 
#' @rdname RNAmodR-Accessors
#' @export
setGeneric( 
  name = "getOutputFolder",
  def = function(.Object ) standardGeneric("getOutputFolder")
) 
#' @rdname RNAmodR-Accessors
#' @export
setGeneric( 
  name = "getExperimentName",
  def = function(.Object ) standardGeneric("getExperimentName")
) 
#' @rdname RNAmodR-Accessors
#' @export
setGeneric( 
  name = "getGFFFile",
  def = function(.Object ) standardGeneric("getGFFFile")
)
#' @rdname RNAmodR-Accessors
#' @export
setGeneric( 
  name = "getFastaFile",
  def = function(.Object ) standardGeneric("getFastaFile")
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
  def = function(experimentName) standardGeneric("setupWorkEnvir")
) 

# RNAModRargs ------------------------------------------------------------------

#' @rdname RNAmodR-args-class
#' @export
setGeneric(
  name = "getArgs",
  def = function(x,
                 identifier,
                 param) standardGeneric("getArgs")
)
#' @rdname RNAmodR-args-class
#' @export
setGeneric(
  name = "loadQuantifier",
  def = function(x) standardGeneric("loadQuantifier")
)
#' @rdname RNAmodR-args-class
#' @export
setGeneric(
  name = "loadIdentifier",
  def = function(x) standardGeneric("loadIdentifier")
)



# experiment data---------------------------------------------------------------

#' @name getExperimentData
#' @export
setGeneric( 
  name = "getExperimentData",
  def = function(.Object,
                 number) standardGeneric("getExperimentData")
) 
#' @name getGffResult
#' @export
setGeneric( 
  name = "getGffResult",
  def = function(.Object,
                 number, 
                 genomicCoordinates = FALSE) standardGeneric("getGffResult")
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
                 number) standardGeneric("setGffResult")
) 




# parsing ----------------------------------------------------------------------

#' @name parseForModifications
#' @export
setGeneric( 
  name = "parseForModifications",
  def = function(x,
                 number,
                 param) standardGeneric("parseForModifications")
) 
#' @name parseForModifications
#' @export
setGeneric( 
  name = "parseForModificationsWithArgs",
  def = function(x,
                 name,
                 gff,
                 fasta) standardGeneric("parseForModificationsWithArgs")
)

# quantifier class accessors ---------------------------------------------------

#' @name RNAmodRquant-accessors
#' @export
setGeneric( 
  name = "getPlotType",
  def = function(x) standardGeneric("getPlotType")
) 

#' @rdname RNAmodRquant-accessors
#' @export
setGeneric( 
  name = "getPositions",
  def = function(x) standardGeneric("getPositions")
) 

#' @rdname RNAmodRquant-accessors
#' @export
setGeneric( 
  name = "getModifications",
  def = function(x) standardGeneric("getModifications")
) 

#' @rdname RNAmodRquant-accessors
#' @export
setGeneric( 
  name = "getDataLabel",
  def = function(x) standardGeneric("getDataLabel")
) 

#' @rdname RNAmodRquant-accessors
#' @export
setGeneric( 
  name = "getDataFormat",
  def = function(x) standardGeneric("getDataFormat")
) 

# quantifier -------------------------------------------------------------------

#' @name quantifyReadData
#' @export
setGeneric( 
  name = "quantifyReadData",
  def = function(x,
                 args,
                 gff,
                 param) standardGeneric("quantifyReadData")
) 

#' @name quantifiyReadDataPerTranscript
#' @title Returns the data for each transcript
setGeneric ( 
  name = "quantifiyReadDataPerTranscript", 
  def = function(x,
                 bamData,
                 args,
                 counts,
                 gff) standardGeneric("quantifiyReadDataPerTranscript")
) 


# ident class accessors --------------------------------------------------------

#' @name RNAmodRident-accessors
#' @export
setGeneric( 
  name = "getModType",
  def = function(x) standardGeneric("getModType")
) 
#' @rdname mod-accessors
#' @export
setGeneric( 
  name = "getDataType",
  def = function(x) standardGeneric("getDataType")
) 
#' @rdname mod-accessors
#' @export
setGeneric( 
  name = "getPositionOffset",
  def = function(x) standardGeneric("getPositionOffset")
) 

# identifier -------------------------------------------------------------------

#' @name identifyModification
#' @export
setGeneric( 
  name = "identifyModifications",
  def = function(x,
                 args,
                 data,
                 gff,
                 fafile) standardGeneric("identifyModifications")
) 

#' @name identifyModificationPerTranscript
#' @title Returns the data for each transcript
setGeneric ( 
  name = "identifyModificationsPerTranscript", 
  def = function(x,
                 location,
                 locations,
                 data) standardGeneric("identifyModificationsPerTranscript")
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
                 focus = FALSE) standardGeneric("getModPlot")
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
                 filetype = "pdf") standardGeneric("saveModPlot")
)

#' @name heatmapModifications
#' @export
setGeneric(
  name = "heatmapModifications",
  def = function(x,
                 numbers,
                 sampleNames,
                 modifications,
                 genes) standardGeneric("heatmapModifications")
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
                 folder) standardGeneric("saveHeatmapModifications")
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
                 modifications) standardGeneric("compareSample")
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
                 genes) standardGeneric("compareAsHeatmap")
)


# SummarizedExperiment ---------------------------------------------------------


#' @name seGetModifications
#' @export
setGeneric( 
  name = "seGetModifications",
  def = function(se,
                 genes,
                 modifications) standardGeneric("seGetModifications")
)
#' @rdname seGetModifications
#' @export
setGeneric( 
  name = "seGetPositions",
  def = function(se,
                 genes) standardGeneric("seGetPositions")
)