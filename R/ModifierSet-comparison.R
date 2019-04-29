#' @include RNAmodR.R
#' @include ModifierSet-class.R
NULL

#' @name compareByCoord
#' @aliases compare compareByCoord plotCompare plotCompareByCoord 
#' 
#' @title Comparison of Samples
#' 
#' @description 
#' To compare data of different samples, a
#' \code{\link[=ModifierSet-class]{ModifierSet}} can be used. To select the data
#' alongside the transcripts and their positions a
#' \code{\link[GenomicRanges:GRanges-class]{GRanges}} or a
#' \code{\link[GenomicRanges:GRanges-class]{GRangesList}} needs to be provided.
#' In case of a \code{GRanges} object, the parent column must match the
#' transcript names as defined by the out put of \code{ranges(x)}, whereas in
#' case of a \code{GRangesList} the element names must match the transcript
#' names.
#' 
#' @param x a \code{Modifier} or \code{ModifierSet} object.
#' @param coord coordinates of position to subset to. Either a \code{GRanges} or
#'   a \code{GRangesList} object. For both types the 'Parent' column is expected
#'   to match the transcript name. The \code{GRangesList} object is
#'   unlisted and only non duplicated entries are retained.
#' @param name Only for \code{compare}: the transcript name
#' @param from Only for \code{compare}: start position
#' @param to Only for \code{compare}: end position
#' @param normalize either a single logical or character value. If it is a 
#' character, it must match one of the names in the \code{ModifierSet}.
#' @param ... optional parameters:
#' \itemize{
#' \item{\code{alias}} {a data.frame with two columns, \code{tx_id} and 
#' \code{name}, to convert transcipt ids to another identifier}
#' \item{\code{name}} {Limit results to one specific gene or transcript}
#' \item{\code{sequenceData}} {TRUE or FALSE? Should the aggregate of 
#' sequenceData be used for the comparison instead of the aggregate data if each
#' \code{Modifier} element? (default: \code{sequenceData = FALSE})}
#' \item{\code{compareType}} {a valid score type to use for the comparison. If
#' \code{sequenceData = FALSE} this defaults to \code{mainScore(x)}, whereas
#' if \code{sequenceData = TRUE} all columns will be used by setting 
#' \code{allTypes = TRUE}.}
#' \item{\code{allTypes}} {TRUE or FALSE? Should all available score be 
#' compared? (default: \code{allTypes = sequenceData})}
#' \item{...} {passed on to \code{\link{subsetByCoord}}}
#' }
#' 
#' @return \code{compareByCoord} returns a
#'   \code{\link[S4Vectors:DataFrame-class]{DataFrame}} and
#'   \code{plotCompareByCoord} returns a \code{ggplot} object, which can be
#'   modified further. The \code{DataFrame} contains columns per sample as well
#'   as the columns \code{names}, \code{positions} and \code{mod} incorporated
#'   from the \code{coord} input. If \code{coord} contains a column
#'   \code{Activity} this is included in the results as well.
#' 
#' @examples
#' data(msi,package="RNAmodR")
#' # constructing a GRanges obejct to mark positive positions
#' mod <- modifications(msi)
#' coord <- unique(unlist(mod))
#' coord$score <- NULL
#' coord$sd <- NULL
#' # return a DataFrame
#' compareByCoord(msi,coord)
#' # plot the comparison as a heatmap
#' plotCompareByCoord(msi,coord)
NULL

.norm_alias <- function(input, x){
  alias <- NULL
  if(!is.null(input[["alias"]])){
    alias <- input[["alias"]]
    if(!is.data.frame(alias)){
      stop("'alias' has to be a data.frame with 'tx_id' and 'name' columns.",
           call. = FALSE)
    }
    colnames <- c("tx_id","name")
    if(!all(colnames %in% colnames(alias))){
      stop("'alias' has to be a data.frame with 'tx_id' and 'name' columns.",
           call. = FALSE)
    }
    alias <- alias[,colnames]
    if(any(duplicated(alias$tx_id))){
      stop("Values in 'tx_id' have to be unique.",
           call. = FALSE)
    }
    names <- names(ranges(x))
    if(!all(alias$tx_id %in% names)){
      stop("All values in 'tx_id' have to be valid transcript ids used as ",
           "names for the data.", call. = FALSE)
    }
  }
  list(alias = alias)
}

.norm_compare_args <- function(input, data, x){
  if(is(x,"ModifierSet")){
    compareType <- mainScore(x)
  } else {
    compareType <- NA
  }
  allTypes <- FALSE
  perTranscript <- FALSE
  sequenceData <- FALSE
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
  if(!is.null(input[["compareType"]])){
    compareType <- input[["compareType"]]
    colnames <- unique(unlist(colnames(data[[1]])))
    if(!is.character(compareType) || width(compareType) == 0L ||
       !(compareType %in% colnames)){
      stop("'compareType' must be a character and a valid colname in the ",
           "aggregated data of 'x'.", call. = FALSE)
    }
  }
  if(!is.null(input[["allTypes"]])){
    allTypes <- input[["allTypes"]]
    if(length(allTypes) != 1L ||
       !is.logical(allTypes)){
      stop("'allTypes' must be a single logical value.", call. = FALSE)
    }
  }
  if(allTypes){
    compareType <- names(data[[1]][[1]])
  }
  if(is.na(compareType[1L])){
    stop("'compareType' must be set if 'sequenceData = TRUE' and ",
         "'allTypes = FALSE'", call. = FALSE)
  }
  args <- c(.norm_alias(input, x),
            list(compareType = compareType,
                 perTranscript = perTranscript,
                 sequenceData = sequenceData))
  args
}

.compare_ModifierSet_by_GRangesList <- function(x, coord, normalize, ...){
  coord <- unlist(coord)
  coord <- unname(coord[!duplicated(coord)])
  .compare_ModifierSet_by_GRanges(x, coord, normalize, ...)
}

.assemble_data_per_compare_type <- function(data, coord, sampleNames, alias, 
                                            modType, normalize){
  data <- do.call(cbind,data)
  colnames(data) <- sampleNames
  if(!is(data,"CompressedSplitDataFrameList")){
    data <- IRanges::SplitDataFrameList(data)
  }
  coord <- coord[match(names(data), names(coord))]
  # keep rownames/names and unlist data
  positions <- rownames(data)
  names <- as.character(S4Vectors::Rle(names(data), lengths(data)))
  data <- unlist(data)
  # add names and positions column as factors
  data$names <- factor(names)
  data$positions <- factor(as.integer(unlist(positions)))
  rownames(data) <- NULL
  # add activity information if present
  coord <- unlist(coord)
  if(any(!is.na(modType))){
    coord <- coord[coord$mod %in% modType,]
  }
  if(!is.null(coord$Activity) || !is.null(coord$mod)){
    f <- unlist(lapply(unique(as.character(data$names)),
                       function(n){
                         d <- data[data$names == n,]
                         d$positions %in% start(coord)[coord$Parent == n]
                       }))
    if(!is.null(coord$Activity)){
      data$Activity <- ""
      data$Activity[f] <- unlist(lapply(coord$Activity, paste, collapse = "/"))
    }
    if(!is.null(coord$mod)){
      data$mod <- ""
      data$mod[f] <- unlist(coord$mod)
    }
  }
  # convert ids to names for labeling if present
  if(!is.null(alias)){
    m <- match(as.character(alias$tx_id),names(data))
    names(data)[m[!is.na(m)]] <- as.character(alias$name)[!is.na(m)]
  }
  #
  data <- .normlize_data_against_one_sample(data, normalize)
  data
}

.compare_ModifierSet_by_GRanges <- function(x, coord, normalize, ...){
  coord <- .norm_coord(coord, modType(x))
  data <- subsetByCoord(x, coord, ...)
  args <- .norm_compare_args(list(...), data, x)
  # restructure to different compare types
  sampleNames <- names(data)
  compareTypes <- args[["compareType"]]
  if(args[["sequenceData"]]){
    modType <- NA
  } else {
    modType <- modType(x)
  }
  data <- lapply(compareTypes,
                 function(ct){
                   lapply(data,
                          function(d){
                            d[,ct,drop = FALSE]
                          })
                 })
  names(data) <- compareTypes
  data <- lapply(data, .assemble_data_per_compare_type, coord, sampleNames,
                 args[["alias"]], modType, normalize)
  if(length(data) == 1L){
    return(data[[1L]])
  }
  data
}


#' @rdname compareByCoord
#' @export
setMethod("compareByCoord",
          signature = c("ModifierSet","GRanges"),
          function(x, coord, normalize, ...){
            .compare_ModifierSet_by_GRanges(x, coord, normalize, ...)
          }
)

#' @rdname compareByCoord
#' @export
setMethod("compareByCoord",
          signature = c("ModifierSet","GRangesList"),
          function(x, coord, normalize, ...){
            .compare_ModifierSet_by_GRangesList(x, coord, normalize, ...)
          }
)

.normlize_data_against_one_sample <- function(data, normalize){
  
  if(!missing(normalize)){
    colnames <- colnames(data)
    colnames <- colnames[!(colnames %in% c("positions","names","mod","Activity"))]
    if(is.character(normalize)){
      assertive::assert_is_a_string(normalize)
      if(!(normalize %in% colnames)){
        stop("Data column '",normalize,"' not found in data. Available columns",
             " are '",paste(colnames, collapse = "','"),"'.",
             call. = FALSE)
      }
      data[,colnames] <- as.data.frame(data[,colnames,drop = FALSE]) - 
        data[,normalize]
    } else if(is.logical(normalize)){
      assertive::assert_is_a_bool(normalize)
      if(normalize == TRUE){
        data[,colnames] <- as.data.frame(data[,colnames,drop = FALSE]) - 
          apply(data[,colnames],1,max)
      }
    } else {
      stop("'normalize' must be a single character or a logical value.",
           call. = FALSE)
    }
  }
  data
}

.norm_compare_plot_args <- function(input){
  limits <- NA
  if(!is.null(input[["limits"]])){
    limits <- input[["limits"]]
    if(!is.numeric(limits) | length(limits) != 2L){
      stop("'limits' must be numeric vector with the length == 2.",
           call. = FALSE)
    }
  }
  args <- list(limits = limits)
  args
}

.create_position_labels <- function(positions, mod, activity){
  if(is.factor(positions)){
    positions <- as.numeric(as.character(positions))
  }
  list <- list(as.character(positions),
               mod,
               activity)
  spacer <- lapply(list,
                   function(el){
                     if(is.null(el)) return(NULL)
                     length <- nchar(el)
                     missingLength <- max(length) - length
                     unlist(lapply(missingLength,
                                   function(n){
                                     paste0(rep(" ",n),collapse = "")
                                   }))
                   })
  sep <- lapply(seq_along(list),
                   function(i){
                     if(i > 1L){
                       rep(" - ",length(list[[i]]))
                     } else {
                       rep("",length(list[[i]]))
                     }
                   })
  labels <- mapply(paste0, spacer, list, sep, SIMPLIFY = FALSE)
  labels <- Reduce(paste0, rev(labels))
  f <- factor(labels, levels = unique(labels))
  stats::reorder(f,positions)
}

.create_sample_labels <- function(labels){
  labels <- as.character(labels)
  labels <- gsub("\\.", " ",labels)
  factor(labels, levels = unique(labels))
}

.plot_compare_ModifierSet_by_GRangesList <- function(x, coord, normalize, ...){
  coord <- unlist(coord)
  coord <- unname(coord[!duplicated(coord)])
  .plot_compare_ModifierSet_by_GRanges(x, coord, normalize, ...)
}

#' @importFrom ggplot2 ggplot geom_raster
#' @importFrom reshape2 melt
.plot_compare_ModifierSet_by_GRanges <- function(x, coord, normalize,  ...){
  args <- .norm_compare_plot_args(list(...))
  data <- .compare_ModifierSet_by_GRanges(x, coord, normalize, ...)
  data$labels <- .create_position_labels(data$positions, data$mod,
                                         data$Activity)
  # melt data an plot
  data$labels <- factor(data$labels, levels = rev(levels(data$labels)))
  data$positions <- NULL
  data$mod <- NULL
  data$Activity <- NULL
  data <- reshape2::melt(as.data.frame(data), id.vars = c("names","labels"))
  data$variable <- .create_sample_labels(data$variable)
  # adjust limits
  if(!missing(normalize) && normalize != FALSE){
    max <- max(max(data$value),abs(min(data$value)))
    max <- max(max,0.5)
    max <- round(max,1) + 0.1
    limits <- c(-max,max)
  } else {
    limits <- c(0,ceiling(max(data$value)))
  }
  if(!is.na(args[["limits"]])){
    limits <- args[["limits"]]
  }
  # plot
  ggplot2::ggplot(data) + 
    ggplot2::geom_raster(mapping = ggplot2::aes_(x = ~variable,
                                                 y = ~labels,
                                                 fill = ~value)) +
    ggplot2::facet_grid(names ~ ., scales = "free", space = "free") +
    ggplot2::scale_fill_gradientn(name = "Score",
                                  colours = rev(colorRamps::matlab.like(100)),
                                  limits = limits) +
    ggplot2::scale_y_discrete(name = "Positions",
                              expand = c(0,0)) +
    ggplot2::scale_x_discrete(name = "Samples",
                              position = "top",
                              expand = c(0,0)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0),
                   axis.text.x.top = ggplot2::element_text(angle = 30,vjust = 0.5))
}

#' @rdname compareByCoord
#' @export
setMethod("plotCompareByCoord",
          signature = c("ModifierSet","GRanges"),
          function(x, coord, normalize, ...){
            .plot_compare_ModifierSet_by_GRanges(x, coord, normalize, ...)
          }
)

#' @rdname compareByCoord
#' @export
setMethod("plotCompareByCoord",
          signature = c("ModifierSet","GRangesList"),
          function(x, coord, normalize, ...){
            .plot_compare_ModifierSet_by_GRangesList(x, coord, normalize, ...)
          }
)
