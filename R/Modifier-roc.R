#' @include RNAmodR.R
#' @include Modifier-class.R
#' @include ModifierSet-class.R
NULL

#' @name plotROC
#' 
#' @title ROCR functions for \code{Modifier} and \code{ModifierSet} objects
#' 
#' @description 
#' \code{plotROC} streamlines labeling, prediction, performance and plotting
#' to test the peformance of a \code{Modifier} object and the data analyzed via
#' the functionallity from the \code{ROCR} package.
#' 
#' The data from \code{x} will be labeled as positive using the \code{coord}
#' arguments. The other arguments will be passed on to the specific \code{ROCR}
#' functions.
#' 
#' By default the \code{prediction.args} include three values:
#' \itemize{
#' \item{\code{measure = "tpr"}}
#' \item{\code{x.measure = "fpr"}}
#' \item{\code{score = mainScore(x)}}
#' }
#' The remaining arguments are not predefined.
#' 
#' @param x a \code{Modifier} or a \code{ModifierSet} object
#' @param coord coordinates of position to label as positive. Either a 
#' \code{GRanges} or a \code{GRangesList} object. For both types the Parent 
#' column is expected to match the gene or transcript name.
#' @param prediction.args arguments which will be used for calling 
#' \code{\link[ROCR:prediction]{prediction}} form the \code{ROCR} package
#' @param performance.args arguments which will be used for calling 
#' \code{\link[ROCR:performance]{performance}} form the \code{ROCR} package
#' @param plot.args arguments which will be used for calling \code{plot} on the
#' performance object of the \code{ROCR} package
#' @param ... additional arguments
NULL

.norm_prediction_args <- function(input){
  if(missing(input)){
    input <- list()
  }
  args <- input
  args
}

.norm_performance_args <- function(input, x){
  if(missing(input)){
    input <- list()
  }
  measure <- "tpr"
  x.measure <- "fpr"
  score <- mainScore(x)
  if(!is.null(input[["measure"]])){
    measure <- input[["measure"]]
    if(!assertive::is_a_string(measure)){
      stop("'measure' must a single character compatible with ",
           "?ROCR::performance.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["x.measure"]])){
    x.measure <- input[["x.measure"]]
    if(!assertive::is_a_string(x.measure)){
      stop("'x.measure' must a single character compatible with ",
           "?ROCR::performance.",
           call. = FALSE)
    }
  }
  if(!is.null(input[["score"]])){
    score <- input[["score"]]
    if(!assertive::is_a_string(score) ||
       !(score %in% colnames(aggregateData(x)[[1]]))){
      stop("'score' must a single character and a valid column name in ",
           "aggregateData().",
           call. = FALSE)
    }
  }
  args <- list(measure = measure,
               x.measure = x.measure,
               score = score)
  args <- c(args,
            input[!(names(input) %in% names(args))])
  args
}

.norm_plot_args <- function(input){
  if(missing(input)){
    input <- list()
  }
  colorize.palette <- NULL
  if(!is.null(input[["colorize.palette"]])){
    colorize.palette <- input[["colorize.palette"]]
    if(!assertive::is_a_string(colorize.palette)){
      stop("'colorize.palette' must a single character compatible with ",
           "?ROCR::plot.performance.",
           call. = FALSE)
    }
  }
  args <- list(colorize.palette = colorize.palette)
  args
}

.get_prediction_data_Modifier <- function(x, coord){
  data <- .label_Modifier_by_GRangesList(x,coord)
  colnames <- colnames(data@unlistData)
  colnames <- colnames[stringr::str_detect(colnames,"score")]
  data <- lapply(seq_along(colnames),
                 function(i){
                   c <- colnames[i]
                   c <- c("labels",c)
                   d <- data[,c]
                   colnames(d) <- c("labels","predictions")
                   d <- unlist(d)
                   rownames(d) <- NULL
                   d
                 })
  names(data) <- colnames
  data
}

.get_prediction_data_ModifierSet <- function(x, coord){
  browser()
  
}

#' @importFrom graphics par abline title legend plot.new
#' @importFrom colorRamps matlab.like
#' @importFrom ROCR prediction performance
.plot_ROCR <- function(data, prediction.args, performance.args, plot.args){
  # add argument logical vector
  n <- length(data)
  # save mfrow setting
  mfrow_bak <- graphics::par("mfrow")
  mfrow_col <- ceiling(sqrt(n))
  mfrow_row <- ceiling(n / mfrow_col)
  graphics::par(mfrow = c(mfrow_row,mfrow_col))
  n_remaining <- (mfrow_row * mfrow_col) - n
  #
  if(is.null(plot.args[["colorize.palette"]])){
    plot.args[["colorize.palette"]] <- NULL
  }
  #
  mapply(
    function(d,
             name,
             colour,
             prediction.args,
             performance.args,
             plot.args){
      pred <- do.call(ROCR::prediction,
                      c(list(d$predictions,
                             d$labels),
                        prediction.args))
      perf <- do.call(ROCR::performance,
                      c(list(pred),
                        performance.args))
      do.call(graphics::plot,
              c(list(perf,
                     colorize = TRUE,
                     lwd = 3),
                plot.args))
      graphics::abline(a = 0, b = 1)
      graphics::title(main = name)
      auc <- unlist(slot(performance(pred,"auc"),"y.values"))
      auc <- paste(c("AUC = "),round(auc,2L),sep="")
      graphics::legend(0.55, 0.25, auc, border = "white", cex = 1,
                       box.col = "white")
    },
    data,
    names(data),
    MoreArgs = list(prediction.args = prediction.args,
                    performance.args = performance.args,
                    plot.args = plot.args),
    SIMPLIFY = FALSE)
  for(i in seq_len(n_remaining)){
    graphics::plot.new()
  }
  graphics::par(mfrow = mfrow_bak)
  invisible(NULL)
  
}

#' @rdname plotROC
#' @export
setMethod(
  f = "plotROC", 
  signature = signature(x = "Modifier"),
  definition = function(x, coord, prediction.args, performance.args, plot.args){
    data <- .get_prediction_data_Modifier(x, coord)
    .plot_ROCR(data,
               .norm_prediction_args(prediction.args),
               .norm_performance_args(performance.args, x),
               .norm_plot_args(plot.args))
  }
)

#' @rdname plotROC
#' @export
setMethod(
  f = "plotROC", 
  signature = signature(x = "ModifierSet"),
  definition = function(x, coord, prediction.args, performance.args, plot.args){
    browser()
    data <- .get_prediction_data_ModifierSet(x, coord)
    .plot_ROCR(data,
               .norm_prediction_args(prediction.args),
               .norm_performance_args(performance.args),
               .norm_plot_args(plot.args))
  }
)
