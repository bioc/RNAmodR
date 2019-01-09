#' @include RNAmodR.R
#' @include ModifierSet-class.R
NULL

#' @name compareByCoord
#' 
#' @title compareByCoord
#' 
#' @description 
#' title
#' 
#' @param ... optional parameters:
#' \itemize{
#' \item{\code{name}}{Limit results to one specific gene or transcript}
#' \item{...}{passed on to \code{\link{subsetByCoord}}}
#' }
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
  browser()
  data <- subsetByCoord(x,coord,...)
  args <- .norm_compare_args(list(...),data,x)
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
  ranges <- unlist(ranges(x[[1]])[names(data)])
  names(data) <- ranges[match(names(data),ranges$ID)]$Name
  positions <- rownames(data)
  names <- as.character(Rle(names(data),lengths(data)))
  data <- unlist(data)
  data$positions <- unlist(positions)
  data$names <- names
  data$positions <- factor(as.integer(data$positions))
  data$names <- factor(data$names)
  rownames(data) <- NULL
  data
}


#' @rdname compareByCoord
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

#' @rdname compareByCoord
#' @export
setMethod("compareByCoord",
          signature = c("ModifierSet","GRangesList"),
          function(x,
                   coord,
                   ...){
            .compare_ModifierSet_by_GRangesList(x,coord,...)
          }
)


#' @importFrom ggplot2 ggplot geom_raster
#' @importFrom reshape2 melt
.plot_compare_ModifierSet_by_GRangesList <- function(x,coord,...){
  data <- .compare_ModifierSet_by_GRangesList(x,coord,...)
  data
  browser()
  # melt data an plot
  data <- reshape2::melt(as.data.frame(data), id.vars = c("names","positions"))
  ggplot2::ggplot(data) + 
    ggplot2::geom_raster(mapping = ggplot2::aes_(x = ~variable,
                                                 y = ~positions,
                                                 fill = ~value)) +
    ggplot2::facet_grid(names ~ ., scales = "free", space = "free")
}

#' @rdname compareByCoord
#' @export
setMethod("plotCompareByCoord",
          signature = c("ModifierSet","GRanges"),
          function(x,
                   coord,
                   ...){
            coord <- split(coord,
                           coord$Parent)
            .plot_compare_ModifierSet_by_GRangesList(x,coord,...)
          }
)

#' @rdname compareByCoord
#' @export
setMethod("plotCompareByCoord",
          signature = c("ModifierSet","GRangesList"),
          function(x,
                   coord,
                   ...){
            .plot_compare_ModifierSet_by_GRangesList(x,coord,...)
          }
)