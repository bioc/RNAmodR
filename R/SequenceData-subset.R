#' @include RNAmodR.R
#' @include SequenceData-class.R
#' @include SequenceDataList-class.R
#' @include Modifier-subset.R
NULL

.subset_SequenceData_by_GRangesList <- function(x, coord, ...){
  args <- .norm_subset_args(list(...), x)
  # converts everything to a GRangesList
  coord <- .norm_coord(coord, args[["type"]])
  if(args[["rawData"]]){
    data <- .norm_sequence_data(as(x,"SplitDataFrameList"))
  } else {
    data <- .norm_aggregate_data(aggregate(x))
  }
  names <- .get_element_names(data, coord, args[["name"]], args[["type"]])
  data <- data[match(names, names(data))]
  coord <- coord[match(names, names(coord))]
  .perform_subset(data, coord, args[["flanking"]], args[["perTranscript"]])
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
  args <- .norm_subset_args(list(...), x)
  # converts everything to a GRangesList
  coord <- .norm_coord(coord, args[["type"]])
  if(args[["rawData"]]){
    data <- as(x,"SplitDataFrameList")
  } else {
    data <- aggregate(x)
  }
  names <- .get_element_names(data, coord, args[["name"]], args[["type"]])
  data <- data[match(names, names(data))]
  coord <- coord[match(names, names(coord))]
  .perform_label(data, coord)
}

.label_SequenceDataList_by_GRangesList <- function(x, coord, ...){
  args <- .norm_subset_args(list(...), x)
  # converts everything to a GRangesList
  coord <- .norm_coord(coord, args[["type"]])
  ans <- lapply(x,
                function(z){
                  if(args[["rawData"]]){
                    data <- as(z,"SplitDataFrameList")
                  } else {
                    data <- aggregate(z)
                  }
                  names <- .get_element_names(data, coord, args[["name"]], args[["type"]])
                  data <- data[match(names, names(data))]
                  coord <- coord[match(names, names(coord))]
                  .perform_label(data, coord)
                })
  labels <- ans[[1]][,"labels",drop=FALSE]
  ans <- lapply(ans,function(a){a@unlistData[,colnames(a@unlistData) != "labels",drop=FALSE]; a})
  ans <- do.call(cbind,c(ans,list(labels)))
  ans
}

################################################################################

.subset_SequenceDataList_by_GRangesList <- function(x, coord, ...){
  args <- .norm_subset_args(list(...),x)
  coord <- .norm_coord(coord,args[["type"]])
  ans <- lapply(x,
                function(z){
                  coord <- .norm_coord(coord, args[["type"]])
                  if(args[["rawData"]]){
                    data <- as(z,"SplitDataFrameList")
                  } else {
                    data <- aggregate(z)
                  }
                  names <- .get_element_names(data, coord, args[["name"]],
                                              args[["type"]])
                  data <- data[match(names, names(data))]
                  coord <- coord[match(names, names(coord))]
                  .perform_subset(data, coord, args[["flanking"]], 
                                  args[["perTranscript"]])
                })
  ans <- do.call(cbind,ans)
  ans
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
