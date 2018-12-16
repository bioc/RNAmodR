#' @include RNAmodR.R
#' @include ModExperiment-class.R
NULL


setClass("ModExperimentInosine",
         contains = c("ModExperiment"),
         prototype = list(mod = "I",
                          data = "PileupModData"))

ModExperimentInosine <- function(){
  new("ModExperimentInosine")
}

