#' @include RNAmodR.R
#' @include PosData-class.R
NULL

#' @name ModExperiment
#' @title ModExperiment
#' @description 
#' title
NULL

#' @rdname ModExperiment
#' @export
setClass("ModExperiment",
         contains = c("VIRTUAL"),
         slots = c(mod = "character",
                   files = "BamFileList",
                   conditions = "CharacterList",
                   genome = "FaFile",
                   genomeFeatures = "character",
                   mods = "GRanges"))

setMethod(
  f = "initialize", 
  signature = signature(.Object = "ModExperiment"),
  definition = function(.Object,
                        files,
                        genome,
                        ranges,
                        mods = NULL) {
    .Object@files <- files
    .Object@conditions <- names(files)
    .Object@genome <- genome
    .Object@genomeFeatures <- ranges
    return(.Object)
  }
)

#' @rdname ModExperiment
#' @export
setMethod(
  f = "show", 
  signature = signature(object = "ModExperiment"),
  definition = function(object) {
    cat("test")
  }
)
