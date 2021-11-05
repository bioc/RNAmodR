#' @include RNAmodR.R
#' @include SequenceData-class.R
#' @include SequenceDataSet-class.R
#' @include Modifier-subset.R
#' @include settings.R
NULL

# common utility function for subsetting ---------------------------------------

.subset_settings <- data.frame(
  variable = c("name",
               "type",
               "merge",
               "flanking",
               "rawData",
               "perTranscript",
               "sequenceData"),
  testFUN = c(".empty_character",
              ".empty_character",
              ".is_a_bool",
              ".not_integer_bigger_equal_than_zero",
              ".is_a_bool",
              ".is_a_bool",
              ".is_a_bool"),
  errorValue = c(TRUE,
                 TRUE,
                 FALSE,
                 TRUE,
                 FALSE,
                 FALSE,
                 FALSE),
  errorMessage = c("'name' must be a character with a width > 0L.",
                   "'type' must be a character with a width > 0L.",
                   "'merge' must be a single logical value.",
                   "'flanking' must be a single integer value equal or higher than 0L.",
                   "'rawData' must be a single logical value.",
                   "'perTranscript' must be a single logical value.",
                   "'sequenceData' must be a single logical value."),
  stringsAsFactors = FALSE)

.norm_subset_args <- function(input,x){
  name <- NA_character_
  if(is(x,"Modifier") || is(x,"ModifierSet")){
    seqtype <- seqtype(x)
    type <- modType(x)
  } else {
    seqtype <- NA_character_
    type <- NA_character_
  }
  merge <- TRUE
  flanking <- 0L
  perTranscript <- FALSE
  sequenceData <- FALSE
  rawData <- FALSE # only used for subsetting SequenceData
  args <- .norm_settings(input, .subset_settings, name, type, merge, flanking,
                         perTranscript, sequenceData, rawData)
  if(all(!is.na(args[["type"]]))){
    if(any(!.is_valid_modType(args[["type"]], seqtype))){
      stop("'type' must be one or more elements of 'shortName(ModRNAString())'",
           " or 'shortName(ModDNAString())'.",
           call. = FALSE)
    }
  }
  args
}

.norm_data <- function(data){
  rn <- rownames(data)
  if(any(lengths(rn) != lengths(data))){
    stop("rownames() of data has to be set.")
  }
  data
}

.norm_coord <- function(coord, type, merge = TRUE){
  if(is(coord,"GRanges")){
    if(is.null(coord$Parent)){
      stop("Parent column must be present.", call. = FALSE)
    }
  } else if(is(coord,"GRangesList")){
    coord <- unlist(coord, use.names = FALSE)
    coord <- coord[!duplicated(coord)]
    return(.norm_coord(coord, type, merge))
  } else {
    stop("")
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
  levels <- unique(coord$Parent)
  # coord <- coord[order(start(coord))]
  if(merge){
    coord <- split(coord, factor(coord$Parent, levels = levels))
  } else {
    coord <- coord[order(factor(coord$Parent, levels = levels))]
    coord <- split(coord, seq_along(coord))
    names(coord) <- mcols(coord, level="within")[,"Parent"]
  }
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
    intersect_names <- intersect(namesData, namesCoord)
    message <- c("No intersection between names in data of 'x' and Parent in ",
                 "'coord'\n ",messageType)
  } else {
    intersect_names <- Reduce(intersect,
                              list(namesData,namesCoord),name)
    message <- c("No intersection between names in data of 'x', Parent in ",
                 "'coord'\n ",messageType," and the selected name.")
  }
  if(length(intersect_names) == 0L){
    stop(message,
         call. = FALSE)
  }
  intersect_names
}

.check_for_invalid_positions <- function(data, coord){
  available_pos <- rownames(data)
  positions <- start(ranges(coord))
  f_names <- match(names(positions),names(available_pos))
  f_names <- f_names[!is.na(f_names)]
  f <- IRanges::LogicalList(Map(function(i,j){!(i %in% j)},
                                positions,
                                available_pos[f_names]))
  if(!any(all(f))){
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
  if(!all(names(coord) %in% names(data))){
    stop("Length and/or order of data and coord do not match.")
  }
  .check_for_invalid_positions(data, coord)
  m <- unlist(match(names(coord),names(data)))
  # construct flanking vector
  flanking <- seq.int(from = -flanking, to = flanking, by = 1L)
  flank <- IRanges::IntegerList(
    lapply(start(ranges(coord)),
           function(i){
             unique(unlist(lapply(i,
                                  function(j){
                                    j + flanking
                                  })
             ))
           }))
  if(length(flanking) > 1L){
    l <- IRanges::IntegerList(as.list(rownames(data)))
    flank <- flank[flank %in% l[m]]
  }
  flank <- as(flank,"CharacterList")
  ff <- IRanges::LogicalList(mapply(function(r,z){r %in% z},
                                    rownames(data)[m],
                                    flank))
  ans <- data[m][ff,,drop=FALSE]
  if(perTranscript){
    pos <- IRanges::CharacterList(Map(
      function(d,i){
        BiocGenerics::which(rownames(d) %in% i)
      },
      data[m],
      flank))
    rownames(ans) <- pos
  }
  return(ans)
}

# subsetting SequenceData ------------------------------------------------------

.subset_SplitDataFrameList_by_GRangesList <- function(data, coord, ...){
  args <- .norm_subset_args(list(...), NULL)
  data <- .norm_data(data)
  coord <- .norm_coord(coord, NA, args[["merge"]])
  element_names <- .get_element_names(data, coord, args[["name"]],
                                      args[["type"]])
  data <- data[match(element_names, names(data))]
  coord <- coord[names(coord) %in% element_names]
  .perform_subset(data, coord, args[["flanking"]], args[["perTranscript"]])
}

.subset_SequenceData_by_GRangesList <- function(x, coord, ...){
  args <- .norm_subset_args(list(...), x)
  # converts everything to a GRangesList
  coord <- .norm_coord(coord, args[["type"]], args[["merge"]])
  if(args[["rawData"]]){
    data <- relist(as(unlist(x, use.names = FALSE),"DFrame"),
                   IRanges::PartitioningByWidth(x))
    data <- .norm_sequence_data(data)
  } else {
    data <- .norm_aggregate_data(aggregate(x))
  }
  data <- .norm_data(data)
  element_names <- .get_element_names(data, coord, args[["name"]],
                                      args[["type"]])
  data <- data[match(element_names, names(data))]
  coord <- coord[names(coord) %in% element_names]
  .perform_subset(data, coord, args[["flanking"]], args[["perTranscript"]])
}

# subsetting SequenceDataSet ---------------------------------------------------

.subset_SequenceDataSet_by_GRangesList <- function(x, coord, ...){
  args <- .norm_subset_args(list(...),x)
  coord <- .norm_coord(coord, args[["type"]], args[["merge"]])
  ans <- lapply(x, .subset_SequenceData_by_GRangesList, coord, ...)
  ans <- do.call(cbind,ans)
  ans
}

# subsetting SequenceDataList --------------------------------------------------

.subset_SequenceDataList_by_GRangesList <- function(x, coord, ...){
  args <- .norm_subset_args(list(...),x)
  coord <- .norm_coord(coord, args[["type"]])
  ans <- 
    lapply(x,
           function(z){
             if(is(z,"SequenceData")){
               return(.subset_SequenceData_by_GRangesList(z, coord, ...))
             } else if(is(z,"SequenceDataSet")) {
               return(.subset_SequenceDataSet_by_GRangesList(z, coord, ...))
             } else {
               stop("")
             }
           })
  ans <- do.call(cbind,ans)
  ans
}

.construct_coord_from_name_from_to <- function(x, name, pos){
  if(!.is_a_string(name)){
    stop("'name' must be a single character value.", call. = FALSE)
  }
  if(!is.integer(pos)){
    stop("'from' and 'to' must be integer and have the same length.", call. = FALSE)
  }
  coord <- GenomicRanges::GRanges(seqnames = seq_len(pos),
                                  ranges = IRanges::IRanges(pos, pos),
                                  Parent = name)
  if(is(x,"Modifier") || is(x,"ModifierSet")){
    S4Vectors::mcols(coord)$mod <- modType(x)[1L]
  }
  coord <- GenomicRanges::GRangesList(coord)
  names(coord) <- name
  coord
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
setMethod("subset",
          signature = c(x = "SequenceData"),
          function(x, name, pos = 1L, ...){
            coord <- .construct_coord_from_name_from_to(x, name, pos)
            .subset_SequenceData_by_GRangesList(x, coord, ...)
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
setMethod("subset",
          signature = c(x = "SequenceDataSet"),
          function(x, name, pos = 1L, ...){
            coord <- .construct_coord_from_name_from_to(x, name, pos)
            .subset_SequenceDataSet_by_GRangesList(x, coord, ...)
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
setMethod("subset",
          signature = c(x = "SequenceDataList"),
          function(x, name, pos = 1L, ...){
            coord <- .construct_coord_from_name_from_to(x, name, pos)
            .subset_SequenceDataList_by_GRangesList(x, coord, ...)
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
  f_rn <- f_rn[!is.na(f_rn)]
  f_p <- f_p[!is.na(f_p)]
  if(!all(names(rn_d)[f_rn] == names(positions)[f_p])){
    stop("Length and/or order of data and coord do not match.")
  }
  labels <- IRanges::LogicalList(lapply(lengths(data),
                                        function(l){rep(FALSE,l)}))
  labels[f_rn] <- IRanges::LogicalList(Map("%in%",
                                           rn_d[f_rn],
                                           positions[f_p]))
  unlisted_data <- unlist(data, use.names = FALSE)
  unlisted_data$labels <- unlist(labels, use.names = FALSE)
  return(relist(unlisted_data, data))
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
    data <- relist(as(unlist(x, use.names = FALSE),"DFrame"),
                   IRanges::PartitioningByWidth(x))
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
               stop("")
             }
           })
  .keep_one_labels_column(ans)
}

#' @rdname RNAmodR-internals
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
