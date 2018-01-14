#' @include RNAmod.R
NULL

#' @rdname heatmapModifications
#'
#' @aliases heatmapModifications saveHeatmapModifications
#'
#' @title Comparative visualization of modifications
#'
#' @param .Object a RNAmod object.
#' @param se a SummarizedExperiment containg the experimental data.
#' @param modifications name of modification to be used for analysis.
#' @param gene a single gene name
#'
#' @return a heatmap plot for given experiments and modifications positions
#' @export
#'
#' @import ggplot2
#'
#' @examples
#' \donttest{
#' }
setMethod(
  f = "heatmapModifications",
  signature = signature(.Object = "RNAmod",
                        se = "SummarizedExperiment",
                        modifications = "character",
                        gene = "character"),
  definition = function(.Object,
                        se,
                        modifications,
                        genes){
    NULL
  }
)

#' @rdname heatmapModifications
#'
#' @export
setMethod(
  f = "saveHeatmapModifications",
  signature = signature(.Object = "RNAmod",
                        se = "SummarizedExperiment",
                        modifications = "character",
                        gene = "character"),
  definition = function(.Object,
                        se,
                        modifications,
                        genes){
    NULL
  }
)