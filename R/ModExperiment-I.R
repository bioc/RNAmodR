#' @include RNAmodR.R
#' @include ModExperiment-class.R
NULL


setClass("ModExperimentInosine",
         contains = c("ModExperiment"),
         prototype = list(mod = "I",
                          dataClass = "PileupModData"))

setMethod(
  f = "initialize", 
  signature = signature(.Object = "ModExperiment"),
  definition = function(.Object,
                        bamfiles,
                        genome,
                        ranges) {
    .Object <- callNextMethod(.Object,
                              bamfiles,
                              genome,
                              ranges)
    return(.Object)
  }
)

ModExperimentInosine <- function(bamfiles,
                                 genome,
                                 ranges){
  ans <- new("ModExperimentInosine",
             bamfiles,
             genome,
             ranges)
  
  message("Starting to search for ")
}

