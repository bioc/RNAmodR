#' @include RNAmodR.R
#' @include ModifierSet-class.R
NULL

.from_ModifierSet_to_RSE <- function(from){
  mod <- modifications(from)
  coord <- unname(unique(unlist(mod)))
  coord$score <- NULL
  coord$sd <- NULL
  # rowData
  data <- compareByCoord(from,coord)
  ## remove non sample columns
  m <- which(colnames(data) %in% c("names","positions","mod"))
  data <- data[,-m]
  assays <- S4Vectors::SimpleList(data)
  names(assays) <- modType(from)
  # colData
  colData <- do.call(rbind,lapply(settings(from),DataFrame))
  # remove default settings
  colData <- colData[,!(colnames(colData) %in% c("find.mod"))]
  # rowRanges
  ranges <- split(coord[,colnames(mcols(coord)) != "Parent"],
                 factor(coord$Parent,coord$Parent))
  SummarizedExperiment::SummarizedExperiment(assays = assays,
                                             rowRanges = ranges,
                                             colData = colData)
}

setAs("ModifierSet","RangedSummarizedExperiment",
      .from_ModifierSet_to_RSE)

setAs("ModifierSet","SummarizedExperiment",
      .from_ModifierSet_to_RSE)
