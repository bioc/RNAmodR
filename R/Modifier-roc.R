#' @include RNAmodR.R
#' @include Modifier-class.R
#' @include ModifierSet-class.R
NULL


#' @name roc
#' 
#' @title Performance
#' 
#' @description 
#' title
NULL

.norm_prediction_args <- function(input){
  if(missing(input)){
    input <- list()
  }
  args <- list()
  args
}

.norm_performance_args <- function(input){
  if(missing(input)){
    input <- list()
  }
  measure <- "tpr"
  x.measure <- "fpr"
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
               x.measure = x.measure)
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

.get_prediction_data <- function(x,
                                 coord){
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

#' @importFrom colorRamps matlab.like
#' @importFrom ROCR prediction performance
.plot_ROCR <- function(data,
                       prediction.args,
                       performance.args,
                       plot.args){
  # add argument logical vector
  n <- length(data)
  # save mfrow setting
  mfrow_bak <- par("mfrow")
  mfrow_col <- ceiling(sqrt(n))
  mfrow_row <- ceiling(n / mfrow_col)
  par(mfrow = c(mfrow_row,mfrow_col))
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
      do.call(plot,
              c(list(perf,
                     colorize = TRUE,
                     lwd = 3),
                plot.args))
      abline(a=0, b= 1)
      title(main = name)
      auc <- unlist(slot(performance(pred,"auc"),"y.values"))
      auc <- paste(c("AUC = "),round(auc,2L),sep="")
      legend(0.55,
             0.25,
             auc,
             border="white",
             cex=1,
             box.col = "white")
    },
    data,
    names(data),
    MoreArgs = list(prediction.args = prediction.args,
                    performance.args = performance.args,
                    plot.args = plot.args))
  for(i in seq_len(n_remaining)){
    plot.new()
  }
  par(mfrow = mfrow_bak)
  invisible(NULL)
  
}

#' @rdname roc
#' 
#' @export
setMethod(
  f = "plotROC", 
  signature = signature(x = "Modifier"),
  definition = function(x,
                        coord,
                        redo,
                        redo.args,
                        prediction.args,
                        performance.args,
                        plot.args) {
    if(missing(redo) || 
       (assertive::is_a_bool(redo) && redo == FALSE)){
      data <- .get_prediction_data(x,coord)
    } else {
      browser()
    }
    .plot_ROCR(data,
               .norm_prediction_args(prediction.args),
               .norm_performance_args(performance.args),
               .norm_plot_args(plot.args))
  }
)

#' @rdname roc
#' @export
setMethod(
  f = "plotROC", 
  signature = signature(x = "ModifierSet"),
  definition = function(x,
                        coord,
                        redo,
                        redo.args,
                        prediction.args,
                        performance.args,
                        plot.args) {
    browser()
    data <- list()
    .plot_ROCR(data,
               .norm_prediction_args(prediction.args),
               .norm_perfomance_args(performance.args),
               .norm_plot_args(plot.args))
  }
)