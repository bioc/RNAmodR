#' @include RNAmodR.R
#' @include Modifier-class.R
NULL

#' @name subset
#' @aliases subsetByCoord
#' 
#' @title Subsetting data from a \code{Modifier} or \code{ModifierSet} object.
#' 
#' @description 
#' With \code{subsetByCoord} data from a \code{Modifier} or \code{ModifierSet} 
#' object will be subset to position as defined in \code{coord}.
#' 
#' @param x a \code{Modifier} or \code{ModifierSet} object.
#' @param coord coordinates of position to subset to. Either a \code{GRanges} or
#' a \code{GRangesList} object. For both types the Parent column is expected to
#' match the gene or transcript name.
#' @param ... optional parameters:
#' \itemize{
#' \item{\code{name}}{Limit results to one specific gene or transcript}
#' \item{\code{type}}{the modification type used for subsetting. By default this
#' is derived from the \code{modType(x)}, but it can be overwritten using 
#' \code{type}. It must be a valid shortName for a modification according to
#' \code{shortName(ModRNAString())} and of course present in metadata column 
#' \code{mod} of \code{coord}}
#' }
NULL

.norm_subset_args <- function(input,x){
  name <- NULL
  type <- modType(x)
  flank <- 0L
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
    if(!(type %in% shortName(ModRNAString()))){
      stop("'type' must be one or more elements of shortName(ModRNAString()).",
           call. = FALSE)
    }
  }
  if(!is.null(input[["flank"]])){
    flank <- input[["flank"]]
    if(!is.integer(flank) || flank < 0L){
      stop("'flank' must be a single integer value equal or higher than 0L.",
           call. = FALSE)
    }
  }
  args <- list(name = name, type = type, flank = flank)
  args
}

.norm_coord <- function(coord,
                        type){
  f <- LogicalList(lapply(coord,function(c){c$type == "RNAMOD"}))
  coord <- coord[f]
  if(unique(unlist(width(ranges(coord)))) != 1L){
    stop("Elements of type 'RNAMOD' with a width != 1L are not supported.",
         call. = "FALSE")
  }
  f <- LogicalList(lapply(coord,function(c){c$mod %in% type}))
  coord <- coord[f]
  coord <- coord[vapply(coord,function(c){length(c) > 0L},logical(1))]
  if(length(coord) == 0L){
    stop("No modifications of type '",paste(type,collapse = "','"),"' ",
         "found in 'coord'.",
         call. = FALSE)
  }
  coord
}

.get_element_names <- function(data, coord, name, type){
  namesData <- names(data)
  namesCoord <- as.character(names(coord))
  if(is.null(name)){
    names <- intersect(namesData, namesCoord)
    message <- c("No intersection between names in data of 'x' and Parent in ",
                 "'coord'\n for modification type '",
                 paste(type, collapse = "','"),"'.")
  } else {
    names <- Reduce(intersect,
                    list(namesData,namesCoord),name)
    message <- c("No intersection between names in data of 'x', Parent in ",
                 "'coord' for modification type '",
                 paste(type, collapse = "','"),"' and the selected name.")
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
  f <- LogicalList(mapply(function(i,j){i >= j},
                          positions,
                          lengths,
                          SIMPLIFY = FALSE))
  if(!any(lengths(which(f)) > 0L )){
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
    message <- c(message,"and more...")
  }
  stop(message,call. = FALSE)
}

.perform_subset <- function(data, coord, flank = 0L){
  if(!all(names(data) == names(coord))){
    stop("Length and/or order of data and coord do not match.")
  }
  .check_for_invalid_positions(data,coord)
  # construct flanking vector
  flank <- seq.int(from = -flank,to = flank, by = 1L)
  f <- IntegerList(lapply(start(ranges(coord)),
                          function(i){
                            unique(unlist(lapply(i,function(j){j + flank})))
                          }))
  return(data[f])
}

.subset_Modifier_by_GRangesList <- function(x, coord, ...){
  args <- .norm_subset_args(list(...), x)
  coord <- .norm_coord(coord, args[["type"]])
  data <- aggregateData(x)
  names <- .get_element_names(data, coord, args[["name"]], args[["type"]])
  data <- data[match(names, names(data))]
  coord <- coord[match(names, names(coord))]
  .perform_subset(data, coord, args[["flank"]])
}

################################################################################
# This is used for ROC and shares functionality with subsetting
.perform_label <- function(data,
                           coord){
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

.label_Modifier_by_GRangesList <- function(x, coord, ...){
  args <- .norm_subset_args(list(...), x)
  coord <- .norm_coord(coord, args[["type"]])
  data <- aggregateData(x)
  names <- .get_element_names(data, coord, args[["name"]], args[["type"]])
  data <- data[match(names, names(data))]
  coord <- coord[match(names, names(coord))]
  .perform_label(data, coord)
}
################################################################################

#' @rdname subset
#' @export
setMethod("subsetByCoord",
          signature = c("Modifier","GRanges"),
          function(x,
                   coord,
                   ...){
            coord <- split(coord,
                           coord$Parent)
            .subset_Modifier_by_GRangesList(x,coord,...)
          }
)
#' @rdname subset
#' @export
setMethod("subsetByCoord",
          signature = c("Modifier","GRangesList"),
          function(x,
                   coord,
                   ...){
            .subset_Modifier_by_GRangesList(x,coord,...)
          }
)

.subset_ModifierSet_by_GRangesList <- function(x, coord, ...){
  args <- .norm_subset_args(list(...),x)
  coord <- .norm_coord(coord,args[["type"]])
  lapply(x,
         function(z){
           data <- aggregateData(z)
           names <- .get_element_names(data, coord, args[["name"]],
                                       args[["type"]])
           data <- data[match(names, names(data))]
           coord <- coord[match(names, names(coord))]
           .perform_subset(data, coord, args[["flank"]])
         })
}

#' @rdname subset
#' @export
setMethod("subsetByCoord",
          signature = c("ModifierSet","GRanges"),
          function(x,
                   coord,
                   ...){
            coord <- split(coord,
                           coord$Parent)
            .subset_ModifierSet_by_GRangesList(x,coord,...)
          }
)
#' @rdname subset
#' @export
setMethod("subsetByCoord",
          signature = c("ModifierSet","GRangesList"),
          function(x,
                   coord,
                   ...){
            .subset_ModifierSet_by_GRangesList(x,coord,...)
          }
)
