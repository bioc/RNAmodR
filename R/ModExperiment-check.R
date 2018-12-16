#' @include RNAmodR.R
NULL


setMethod(
  f = "match.genome", 
  signature = signature(x = "ModExperiment"),
  definition = function(x,
                        ...){
    ModExperiments <- .get_ModExperiment_objects(...)
    genomes <- lapply(c(list(x),
                        ModExperiments),
                      genome)
    paths <- unique(unlist(lapply(genomes,path)))
    if(length(paths) != 1){
      stop("")
    }
    return(paths)
  }
)

setMethod(
  f = "match.annotation", 
  signature = signature(x = "ModExperiment"),
  definition = function(x,
                        ...){
    ModExperiments <- .get_ModExperiment_objects(...)
    genomeFeatureFiles <- lapply(c(list(x),
                                   ModExperiments),
                                 genomeFeatures)
    genomeFeatureFiles <- unique(unlist(genomeFeatureFiles))
    if(length(genomeFeatureFiles) != 1){
      stop("")
    }
    return(genomeFeatureFiles)
  }
)
