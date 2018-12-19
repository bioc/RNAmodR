#' @include RNAmodR.R
#' @include Modifier-class.R
NULL


# 
# setMethod(
#   f = "getMod", 
#   signature = signature(x = "ModExperiment"),
#   definition = function(x,
#                         ...){
#     # get ModExperiment objects from input
#     ModExperiments <- c(list(x),.get_ModExperiment_objects(...))
#     # check if sequence and annotation are equivalent
#     genome <- match.genome(x, ModExperiments)
#     annotation <- match.annotation(x, ModExperiments)
#     # subset to relevant annotations 
#     annotation_features <- .get_parent_annotations(
#       .subset_rnamod_containing_features(annotation)
#     )
#     # retrieve the data
#     DataInfo <- DataFrame(
#       dataFiles = lapply(lapply(ModExperiments,dataFiles),path),
#       dataTypes = lapply(ModExperiments,slot,"data"),
#       dataParent = NULL)
#     f <- duplicated(DataInfo)
#     DataInfo[f,"Parent"] <- lapply(which(f),
#                                    function(i){
#                                      which(DataInfo[i,] == DataInfo)[1]
#                                    })
#     Data <- mapply(function(exp,type){
#                      do.call(type,dataFiles(exp))
#                    },
#                    DataInfo[!f,]$dataFiles,
#                    DataInfo[!f,]$dataTypes,
#                    MoreArgs = list(genome = genome,
#                                    annotation = annotation))
#     
#     
#     
#     return("fump")
#   }
# )
