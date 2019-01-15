#' @include RNAmodR.R
#' @include ModifierSet-class.R
NULL

#' @name compare
#' @aliases compareByCoord
#' 
#' @title compare
#' 
#' @description 
#' title
#' 
#' @param x a \code{Modifier} or \code{ModifierSet} object.
#' @param coord coordinates of position to subset to. Either a \code{GRanges} or
#' a \code{GRangesList} object. For both types the Parent column is expected to
#' match the gene or transcript name.
#' @param normalize either a single logical or character value. If it is a 
#' character, it must match one of the names in the \code{ModifierSet}.
#' @param ... optional parameters:
#' \itemize{
#' \item{\code{name}}{Limit results to one specific gene or transcript}
#' \item{...}{passed on to \code{\link{subsetByCoord}}}
#' }
#' 
#' @return for \code{compareByCoord} as \code{\link{DataFrameList}} and for
#' \code{plotCompareByCoord} a \code{ggplot} object, which can be modified 
#' further.
NULL

.norm_compare_args <- function(input,data,x){
  compareType <- mainScore(x)
  if(!is.null(input[["compareType"]])){
    compareType <- input[["compareType"]]
    colnames <- unique(unlist(colnames(data[[1]])))
    if(!is.character(compareType) || width(compareType) == 0L ||
       !(compareType %in% colnames)){
      stop("'compareType' must be a character and a valid colname in the 
           aggregated data of 'x'.",
           call. = FALSE)
    }
  }
  args <- list(compareType = compareType)
  args
}

.compare_ModifierSet_by_GRangesList <- function(x,coord,...){
  data <- subsetByCoord(x,coord,...)
  args <- .norm_compare_args(list(...),
                             data,
                             x)
  # subset to compare type
  sampleNames <- names(data)
  data <- lapply(data,
                 function(d){
                   d[,args[["compareType"]],drop = FALSE]
                 })
  data <- do.call(cbind,data)
  if(!is(data,"CompressedSplitDataFrameList")){
    data <- SplitDataFrameList(data)
  }
  colnames(data) <- sampleNames
  # convert ids to names for labeling if present
  ranges <- unlist(ranges(x[[1]])[names(data)])
  if(!is.null(ranges$Name) && 
     any(!is.na(ranges$Name))){
    names <- ranges[match(names(data),ranges$ID)]$Name
    f <- !is.na(names)
    names(data)[f] <- names[f]
  }
  # keep rownames/names and unlist data
  positions <- rownames(data)
  names <- as.character(Rle(names(data),lengths(data)))
  data <- unlist(data)
  # add names and positions column as factors
  data$names <- factor(names)
  data$positions <- factor(as.integer(unlist(positions)))
  rownames(data) <- NULL
  # add activity information if present
  coord <- unlist(coord)
  coord <- coord[coord$mod %in% modType(x),]
  if(!is.null(coord$Activity)){
    data$Activity <- unlist(lapply(coord$Activity,
                                   paste,
                                   collapse = "/"))
  }
  #
  data
}


#' @rdname compare
#' @export
setMethod("compareByCoord",
          signature = c("ModifierSet","GRanges"),
          function(x,
                   coord,
                   ...){
            coord <- split(coord,
                           coord$Parent)
            .compare_ModifierSet_by_GRangesList(x,coord,...)
          }
)

#' @rdname compare
#' @export
setMethod("compareByCoord",
          signature = c("ModifierSet","GRangesList"),
          function(x,
                   coord,
                   ...){
            .compare_ModifierSet_by_GRangesList(x,coord,...)
          }
)

.create_position_labels <- function(list){
  spacer <- lapply(list,
                   function(el){
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
  factor(labels, levels = labels)
}

.create_sample_labels <- function(labels){
  labels <- as.character(labels)
  labels <- gsub("\\.", " ",labels)
  factor(labels)
}

#' @importFrom ggplot2 ggplot geom_raster
#' @importFrom reshape2 melt
.plot_compare_ModifierSet_by_GRangesList <- function(x,
                                                     coord,
                                                     normalize,
                                                     ...){
  data <- .compare_ModifierSet_by_GRangesList(x,coord,...)
  if(!missing(normalize)){
    colnames <- colnames(data)
    colnames <- colnames[!(colnames %in% c("positions","names","Activity"))]
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
  data$labels <- .create_position_labels(list(as.character(data$positions),
                                              data$Activity))
  # melt data an plot
  data$labels <- factor(data$labels, levels = rev(data$labels))
  data$positions <- NULL
  data$Activity <- NULL
  data <- reshape2::melt(as.data.frame(data), id.vars = c("names","labels"))
  data$variable <- .create_sample_labels(data$variable)
  # adjust limits
  if(!missing(normalize) && normalize != FALSE){
    max <- max(max(data$value),abs(min(data$value)))
    max <- max(max,0.5)
    limits <- c(round(-max,1),round(max,1))
  } else {
    limits <- c(NA,round(max(data$value)))
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
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0))
}

#' @rdname compare
#' @export
setMethod("plotCompareByCoord",
          signature = c("ModifierSet","GRanges"),
          function(x,
                   coord,
                   normalize,
                   ...){
            coord <- split(coord,
                           coord$Parent)
            .plot_compare_ModifierSet_by_GRangesList(x,
                                                     coord,
                                                     normalize,
                                                     ...)
          }
)

#' @rdname compare
#' @export
setMethod("plotCompareByCoord",
          signature = c("ModifierSet","GRangesList"),
          function(x,
                   coord,
                   normalize,
                   ...){
            .plot_compare_ModifierSet_by_GRangesList(x,
                                                     coord,
                                                     normalize,
                                                     ...)
          }
)