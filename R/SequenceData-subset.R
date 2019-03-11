#' @include RNAmodR.R
#' @include SequenceData-class.R
#' @include SequenceDataList-class.R
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
    if(!is.character(type) || width(type) == 0L){
      stop("'type' must be a character with a width > 0L.",
           call. = FALSE)
    }
    if(!(type %in% Modstrings::shortName(Modstrings::ModRNAString()))){
      stop("'type' must be one or more elements of shortName(ModRNAString()).",
           call. = FALSE)
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

.norm_coord <- function(coord, type){
  if(is(coord,"GRanges")){
    if(is.null(coord$Parent)){
      stop("Parent column must be present.", call. = FALSE)
    }
    coord <- split(coord, coord$Parent)
  } else if(is(coord,"GRangesList")){
    coord <- unname(unlist(coord))
    coord <- coord[!duplicated(coord)]
    return(.norm_coord(coord,type))
  } else {
    stop("Something went wrong.")
  }
  f <- IRanges::LogicalList(lapply(coord,function(c){c$type == "RNAMOD"}))
  coord <- coord[f]
  if(unique(unlist(width(ranges(coord)))) != 1L){
    stop("Elements of type 'RNAMOD' with a width != 1L are not supported.",
         call. = "FALSE")
  }
  if(any(!is.na(type))){
    f <- IRanges::LogicalList(lapply(coord,function(c){c$mod %in% type}))
    coord <- coord[f]
    coord <- coord[vapply(coord,function(c){length(c) > 0L},logical(1))]
    if(length(coord) == 0L){
      stop("No modifications of type '",paste(type,collapse = "','"),"' ",
           "found in 'coord'.",
           call. = FALSE)
    }
  }
  coord
}

.get_element_names <- function(data, coord, name, type){
  namesData <- names(data)
  namesCoord <- as.character(names(coord))
  if(is.na(type)){
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

.check_for_invalid_positions <- function(data,coord){
  lengths <- lengths(data)
  positions <- start(ranges(coord))
  f <- IRanges::LogicalList(mapply(function(i,j){i >= j},
                                   positions,
                                   lengths,
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
  .check_for_invalid_positions(data,coord)
  # construct flanking vector
  flanking <- seq.int(from = -flanking,to = flanking, by = 1L)
  f <- IRanges::IntegerList(lapply(start(ranges(coord)),
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
  ans <- data[f]
  if(perTranscript){
    pos <- IRanges::CharacterList(mapply(
      function(d,i){
        BiocGenerics::which(i == rownames(d))
      },
      data,
      f,
      SIMPLIFY = FALSE))
    rownames(ans) <- pos
  }
  return(ans)
}

# subsetting SequenceData ------------------------------------------------------

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

# subsetting SequenceDataList --------------------------------------------------

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
                  names <- .get_element_names(data, coord, args[["name"]],
                                              args[["type"]])
                  data <- data[match(names, names(data))]
                  coord <- coord[match(names, names(coord))]
                  .perform_label(data, coord)
                })
  labels <- ans[[1]][,"labels",drop=FALSE]
  ans <- lapply(ans,
                function(a){
                  a@unlistData <- 
                    a@unlistData[,colnames(a@unlistData) != "labels",drop=FALSE]
                  a
                })
  ans <- do.call(cbind,c(ans,list(labels)))
  ans
}
