#' @include RNAmodR.R
#' @include Modifier-class.R
NULL


setMethod(
  f = "match.genome", 
  signature = signature(x = "Modifier"),
  definition = function(x,
                        ...){
    Modifiers <- .get_Modifier_objects(...)
    genomes <- lapply(c(list(x),
                        Modifiers),
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
  signature = signature(x = "Modifier"),
  definition = function(x,
                        ...){
    Modifiers <- .get_Modifier_objects(...)
    genomeFeatureFiles <- lapply(c(list(x),
                                   Modifiers),
                                 genomeFeatures)
    genomeFeatureFiles <- unique(unlist(genomeFeatureFiles))
    if(length(genomeFeatureFiles) != 1){
      stop("")
    }
    return(genomeFeatureFiles)
  }
)
