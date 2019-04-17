#' @include RNAmodR.R
#' @include SequenceData-class.R
#' @include SequenceDataSet-class.R
#' @include Modifier-subset.R
NULL

# common utility function for subsetting ---------------------------------------

.norm_subset_args <- function(input,x){
  name <- NA_character_
  if(is(x,"Modifier") || is(x,"ModifierSet")){
    type <- modType(x)
  } else {
    type <- NA_character_
  }
  flanking <- 0L
  perTranscript <- FALSE
  sequenceData <- FALSE
  rawData <- FALSE # only used for subsetting SequenceData
  if(!is.null(input[["name"]])){
    name <- input[["name"]]
    if(!is.character(name) || width(name) == 0L){
      stop("'name' must be a character with a width > 0L.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["type"]])){
    type <- input[["type"]]
    if(!is.na(type)){
      if(!is.character(type) || width(type) == 0L){
        stop("'type' must be a character with a width > 0L.",
             call. = FALSE)
      }
      if(!(type %in% Modstrings::shortName(Modstrings::ModRNAString()))){
        stop("'type' must be one or more elements of shortName(ModRNAString()).",
             call. = FALSE)
      }
    }
  }
  if(!is.null(input[["flanking"]])){
    flanking <- input[["flanking"]]
    if(!is.integer(flanking) || flanking < 0L){
      stop("'flanking' must be a single integer value equal or higher than 0L.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["rawData"]])){ # only used for subsetting SequenceData
    rawData <- input[["rawData"]]
    if(!assertive::is_a_bool(rawData)){
      stop("'rawData' must be a single logical value.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["perTranscript"]])){
    perTranscript <- input[["perTranscript"]]
    if(!assertive::is_a_bool(perTranscript)){
      stop("'perTranscript' must be a single logical value.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["sequenceData"]])){
    sequenceData <- input[["sequenceData"]]
    if(!assertive::is_a_bool(sequenceData)){
      stop("'sequenceData' must be a single logical value.")
    }
  }
  args <- list(name = name, type = type, flanking = flanking,
               rawData = rawData, perTranscript = perTranscript,
               sequenceData = sequenceData)
  args
}

.norm_data <- function(data){
  rn <- rownames(data)
  if(any(lengths(rn) != lengths(data))){
    stop("rownames() of data has to be set.")
  }
  data
}

.norm_coord <- function(coord, type){
  if(is(coord,"GRanges")){
    if(is.null(coord$Parent)){
      stop("Parent column must be present.", call. = FALSE)
    }
  } else if(is(coord,"GRangesList")){
    coord <- unlist(coord, use.names = FALSE)
    coord <- coord[!duplicated(coord)]
    return(.norm_coord(coord,type))
  } else {
    stop("Something went wrong.")
  }
  if(unique(unlist(width(ranges(coord)))) != 1L){
    stop("Elements with a width != 1L are not supported.",
         call. = "FALSE")
  }
  if("mod" %in% colnames(S4Vectors::mcols(coord))){
    if(any(!is.na(type))){
      coord <- coord[!is.na(coord$mod)]
      coord <- coord[coord$mod %in% type]
      if(length(coord) == 0L){
        stop("No modifications of type '",paste(type,collapse = "','"),"' ",
             "found in 'coord'.",
             call. = FALSE)
      }
    }
  }
  coord <- split(coord, factor(coord$Parent, levels = unique(coord$Parent)))
  coord
}

.get_element_names <- function(data, coord, name, type){
  namesData <- names(data)
  namesCoord <- as.character(names(coord))
  if(any(is.na(type))){
    messageType <- "for any modification type"
  } else {
    messageType <- paste0("for modification type '",
                          paste(type, collapse = "','"), "'")
  }
  if(is.na(name)){
    names <- intersect(namesData, namesCoord)
    message <- c("No intersection between names in data of 'x' and Parent in ",
                 "'coord'\n ",messageType)
  } else {
    names <- Reduce(intersect,
                    list(namesData,namesCoord),name)
    message <- c("No intersection between names in data of 'x', Parent in ",
                 "'coord'\n ",messageType," and the selected name.")
  }
  if(length(names) == 0L){
    stop(message,
         call. = FALSE)
  }
  names
}

.check_for_invalid_positions <- function(data, coord){
  lengths <- lengths(data)
  positions <- start(ranges(coord))
  f_names <- names(lengths) %in% names(positions)
  f <- IRanges::LogicalList(mapply(function(i,j){i >= j},
                                   positions,
                                   lengths[f_names],
                                   SIMPLIFY = FALSE))
  if(!any(lengths(BiocGenerics::which(f)) > 0L )){
    return(NULL)
  }
  invalidPositions <- unlist(lapply(coord[f],as.character),
                             use.names = FALSE)
  invalidTypes <- unlist(lapply(coord[f],function(c){c$mod}),
                         use.names = FALSE)
  if(length(invalidPositions) > 10L){
    i <- seq_len(10L)
  } else {
    i <- seq_along(invalidPositions)
  }
  message <- c("'coord' for the following modifications out of range:\n",
               paste0(invalidPositions[i]," for '",invalidTypes[i],"'",
                      collapse = "\n"))
  if(length(invalidPositions) > 10L){
    message <- c(message," and more...")
  }
  stop(message,call. = FALSE)
}

.perform_subset <- function(data, coord, flanking = 0L, perTranscript = FALSE){
  if(!all(names(data) == names(coord))){
    stop("Length and/or order of data and coord do not match.")
  }
  .check_for_invalid_positions(data, coord)
  # construct flanking vector
  flanking <- seq.int(from = -flanking, to = flanking, by = 1L)
  f <- IRanges::CharacterList(
    lapply(start(ranges(coord)),
           function(i){
             unique(unlist(lapply(i,
                                  function(j){
                                    j + flanking
                                  })
             ))
           }))
  if(length(flanking) > 1L){
    f <- f[f > 0L & f <= lengths(data)]
  }
  ff <- IRanges::LogicalList(mapply(function(r,z){r %in% z},rownames(data),f))
  ans <- data[ff]
  if(perTranscript){
    pos <- IRanges::CharacterList(mapply(
      function(d,i){
        BiocGenerics::which(rownames(d) %in% i)
      },
      data,
      f,
      SIMPLIFY = FALSE))
    rownames(ans) <- pos
  }
  return(ans)
}

# subsetting SequenceData ------------------------------------------------------

.subset_SplitDataFrameList_by_GRangesList <- function(data, coord, ...){
  args <- .norm_subset_args(list(...), NULL)
  data <- .norm_data(data)
  coord <- .norm_coord(coord, NA)
  names <- .get_element_names(data, coord, args[["name"]], args[["type"]])
  data <- data[match(names, names(data))]
  coord <- coord[match(names, names(coord))]
  .perform_subset(data, coord, args[["flanking"]], args[["perTranscript"]])
}

.subset_SequenceData_by_GRangesList <- function(x, coord, ...){
  args <- .norm_subset_args(list(...), x)
  # converts everything to a GRangesList
  coord <- .norm_coord(coord, args[["type"]])
  if(args[["rawData"]]){
    data <- .norm_sequence_data(as(x,"SplitDataFrameList"))
  } else {
    data <- .norm_aggregate_data(aggregate(x))
  }
  data <- .norm_data(data)
  names <- .get_element_names(data, coord, args[["name"]], args[["type"]])
  data <- data[match(names, names(data))]
  coord <- coord[match(names, names(coord))]
  .perform_subset(data, coord, args[["flanking"]], args[["perTranscript"]])
}

# subsetting SequenceDataSet ---------------------------------------------------

.subset_SequenceDataSet_by_GRangesList <- function(x, coord, ...){
  args <- .norm_subset_args(list(...),x)
  coord <- .norm_coord(coord,args[["type"]])
  ans <- lapply(x, .subset_SequenceData_by_GRangesList, coord, ...)
  ans <- do.call(cbind,ans)
  ans
}

# subsetting SequenceDataList --------------------------------------------------

.subset_SequenceDataList_by_GRangesList <- function(x, coord, ...){
  args <- .norm_subset_args(list(...),x)
  coord <- .norm_coord(coord,args[["type"]])
  ans <- 
    lapply(x,
           function(z){
             if(is(z,"SequenceData")){
               return(.subset_SequenceData_by_GRangesList(z, coord, ...))
             } else if(is(z,"SequenceDataSet")) {
               return(.subset_SequenceDataSet_by_GRangesList(z, coord, ...))
             } else {
               stop("Something went wrong.")
             }
           })
  ans <- do.call(cbind,ans)
  ans
}

#' @rdname subsetByCoord
#' @export
setMethod("subsetByCoord",
          signature = c(x = "SplitDataFrameList", coord = "GRanges"),
          function(x, coord, ...){
            .subset_SplitDataFrameList_by_GRangesList(x, coord, ...)
          }
)
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
          signature = c(x = "SequenceDataSet", coord = "GRanges"),
          function(x, coord, ...){
            .subset_SequenceDataSet_by_GRangesList(x, coord, ...)
          }
)
#' @rdname subsetByCoord
#' @export
setMethod("subsetByCoord",
          signature = c(x = "SequenceDataSet", coord = "GRangesList"),
          function(x, coord, ...){
            .subset_SequenceDataSet_by_GRangesList(x, coord, ...)
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


# labeling ---------------------------------------------------------------------

.perform_label <- function(data, coord){
  .check_for_invalid_positions(data,coord)
  positions <- start(ranges(coord))
  rn_d <- rownames(data)
  f_rn <- match(names(positions), names(rn_d))
  f_p <- match(names(rn_d)[f_rn], names(positions))
  if(!all(names(rn_d)[f_rn] == names(positions)[f_p])){
    stop("Length and/or order of data and coord do not match.")
  }
  labels <- IRanges::LogicalList(lapply(lengths(data),
                                        function(l){rep(FALSE,l)}))
  labels[f_rn] <- IRanges::LogicalList(mapply("%in%",
                                              rn_d[f_rn],
                                              positions[f_p],
                                              SIMPLIFY = FALSE))
  unlisted_data <- unlist(data, use.names = FALSE)
  unlisted_data$labels <- unlist(labels, use.names = FALSE)
  return(relist(unlisted_data,data))
}

.label_SplitDataFrameList_by_GRangesList <- function(data, coord, ...){
  data <- .norm_data(data)
  coord <- .norm_coord(coord, NA)
  .perform_label(data, coord)
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
  data <- .norm_data(data)
  .perform_label(data, coord)
}

.keep_one_labels_column <- function(data){
  labels <- data[[1]][,"labels",drop=FALSE]
  data <- lapply(data,
                 function(d){
                   d@unlistData <- 
                     d@unlistData[,colnames(d@unlistData) != "labels",
                                  drop=FALSE]
                   d
                 })
  data <- do.call(cbind,c(data,list(labels)))
  data
}

.label_SequenceDataSet_by_GRangesList <- function(x, coord, ...){
  args <- .norm_subset_args(list(...), x)
  # converts everything to a GRangesList
  coord <- .norm_coord(coord, args[["type"]])
  ans <- lapply(x, .label_SequenceData_by_GRangesList, coord, ...)
  .keep_one_labels_column(ans)
}

.label_SequenceDataList_by_GRangesList <- function(x, coord, ...){
  args <- .norm_subset_args(list(...), x)
  # converts everything to a GRangesList
  coord <- .norm_coord(coord, args[["type"]])
  ans <- 
    lapply(x,
           function(z){
             if(is(z,"SequenceData")){
               return(.label_SequenceData_by_GRangesList(z, coord, ...))
             } else if(is(z,"SequenceDataSet")) {
               return(.label_SequenceDataSet_by_GRangesList(z, coord,...))
             } else {
               stop("Something went wrong.")
             }
           })
  .keep_one_labels_column(ans)
}

#' @rdname subsetByCoord
#' @export
setMethod("labelByCoord",
          signature = c(x = "SplitDataFrameList", coord = "GRanges"),
          function(x, coord, ...){
            .label_SplitDataFrameList_by_GRangesList(x, coord, ...)
          }
)
#' @rdname subsetByCoord
#' @export
setMethod("labelByCoord",
          signature = c(x = "SequenceData", coord = "GRanges"),
          function(x, coord, ...){
            .label_SequenceData_by_GRangesList(x, coord, ...)
          }
)
#' @rdname subsetByCoord
#' @export
setMethod("labelByCoord",
          signature = c(x = "SequenceData", coord = "GRangesList"),
          function(x, coord, ...){
            .label_SequenceData_by_GRangesList(x, coord, ...)
          }
)
#' @rdname subsetByCoord
#' @export
setMethod("labelByCoord",
          signature = c(x = "SequenceDataSet", coord = "GRanges"),
          function(x, coord, ...){
            .label_SequenceDataSet_by_GRangesList(x, coord, ...)
          }
)
#' @rdname subsetByCoord
#' @export
setMethod("labelByCoord",
          signature = c(x = "SequenceDataSet", coord = "GRangesList"),
          function(x, coord, ...){
            .label_SequenceDataSet_by_GRangesList(x, coord, ...)
          }
)
#' @rdname subsetByCoord
#' @export
setMethod("labelByCoord",
          signature = c(x = "SequenceDataList", coord = "GRanges"),
          function(x, coord, ...){
            .label_SequenceDataList_by_GRangesList(x, coord, ...)
          }
)
#' @rdname subsetByCoord
#' @export
setMethod("labelByCoord",
          signature = c(x = "SequenceDataList", coord = "GRangesList"),
          function(x, coord, ...){
            .label_SequenceDataList_by_GRangesList(x, coord, ...)
          }
)
