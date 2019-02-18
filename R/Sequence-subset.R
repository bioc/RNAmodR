#' @include RNAmodR.R
#' @include SequenceData-class.R
#' @include SequenceDataList-class.R
#' @include Modifier-subset.R
NULL

.subset_SequenceData_by_GRangesList <- function(x, coord, ...){
}

################################################################################
# This is used for ROC and shares functionality with subsetting
.perform_label <- function(data, coord){
  .check_for_invalid_positions(data,coord)
  # add positions as rownames
  rownames(data@unlistData) <- unlist(lapply(lengths(data),seq_len),
                                      use.names = FALSE)
  # 
  lengths <- lengths(data)
  positions <- start(ranges(coord))
  labels <- IRanges::LogicalList(lapply(lengths,function(l){rep(FALSE,l)}))
  labels <- IRanges::LogicalList(mapply(
    function(l,p){
      l[p] <- TRUE
      l
    },
    labels,
    positions,
    SIMPLIFY = FALSE))
  data@unlistData$labels <- unlist(labels)
  return(data)
}

.label_SequenceData_by_GRangesList <- function(x, coord, ...){
}

################################################################################

.subset_SequenceDataList_by_GRangesList <- function(x, coord, ...){
}

#' @rdname subsetByCoord
#' @export
setMethod("subsetByCoord",
          signature = c(x = "SequenceData", coord = "GRanges"),
          function(x, coord, ...){
            .subset_SequenceData_by_GRangesList(x, coord, ...)
          }
)
#' @rdname subsetByCoord
#' @export
setMethod("subsetByCoord",
          signature = c(x = "SequenceData", coord = "GRangesList"),
          function(x, coord, ...){
            .subset_SequenceData_by_GRangesList(x, coord, ...)
          }
)
#' @rdname subsetByCoord
#' @export
setMethod("subsetByCoord",
          signature = c(x = "SequenceDataList", coord = "GRanges"),
          function(x, coord, ...){
            .subset_SequenceDataList_by_GRangesList(x, coord, ...)
          }
)
#' @rdname subsetByCoord
#' @export
setMethod("subsetByCoord",
          signature = c(x = "SequenceDataList", coord = "GRangesList"),
          function(x, coord, ...){
            .subset_SequenceDataList_by_GRangesList(x, coord, ...)
          }
)
