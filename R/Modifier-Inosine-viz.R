#' @include RNAmodR.R
#' @include Modifier-Inosine-class.R
NULL

RNAMODR_I_PLOT_DATA_COLOURS <- c("score" = "#ABABAB") 
RNAMODR_I_PLOT_DATA_NAMES <- c(score = "Score Inosine")

#' @rdname ModInosine-functions
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "ModInosine",
                        coord = "GRanges"),
  definition = function(x, coord, type = "score", window.size = 15L, ...) {
    if(missing(type)){
      type <- "score"
    }
    type <- match.arg(type, "score")
    callNextMethod(x = x, coord = coord, type = type, window.size = window.size,
                   ...)
  }
)
#' @rdname ModInosine-functions
#' @export
setMethod(
  f = "visualizeData",
  signature = signature(x = "ModInosine"),
  definition = function(x, name, from, to, type = "score", ...) {
    callNextMethod(x = x, name, from, to, type = type, ...)
  }
)

.norm_viz_mod_inosine_args <- function(...){
  input <- list(...)
  colour <- input[["colour"]]
  if(is.na(colour) || length(colour) != length(RNAMODR_I_PLOT_DATA_COLOURS)){
    colour <- RNAMODR_I_PLOT_DATA_COLOURS
  }
  input <- list(colour = colour)
  input
}

.clean_mcols_mod_inosine <- function(x, data, name){
  d <- mcols(data@unlistData)
  seq <- sequences(x)
  letters <- unlist(unname(strsplit(as.character(seq),"")))
  f <- d$score < 0
  d$sd[f] <- 0
  d$score[f] <- 0
  d <- d[colnames(d) %in% c("score","sd")]
  d$score[letters != "A"] <- 0
  mcols(data@unlistData) <- d
  data
}

#' @name ModInosine-functions
setMethod(
  f = "getDataTrack",
  signature = signature(x = "ModInosine"),
  definition = function(x, ...) {
    args <- .norm_viz_mod_inosine_args(...)
    name <- .norm_viz_name(args[["name"]])
    data <- .get_data_for_visualization(x, name)
    data <- .clean_mcols_mod_inosine(x, data, name)
    data <- unlist(data)
    lim <- c(min(mcols(data)$score), max(mcols(data)$score))
    dtscore <- Gviz::DataTrack(range = data[,"score"],
                               groups = "score",
                               name = RNAMODR_I_PLOT_DATA_NAMES["score"],
                               col = args[["colour"]],
                               type = "histogram",
                               ylim = lim)
    Gviz::displayPars(dtscore)$background.title <- "#FFFFFF"
    Gviz::displayPars(dtscore)$fontcolor.title <- "#000000"
    Gviz::displayPars(dtscore)$col.axis <- "#000000"
    Gviz::displayPars(dtscore) <- args
    # dtsd <- Gviz::DataTrack(range = data[,"sd"],
    #                            groups = "sd",
    #                            col = args[["colour"]],
    #                            type = "histogram",
    #                            ylim = lim)
    # Gviz::displayPars(dtsd)$background.title <- "#FFFFFF"
    # Gviz::displayPars(dtsd)$fontcolor.title <- "#000000"
    # Gviz::displayPars(dtsd)$col.axis <- "#000000"
    # Gviz::displayPars(dtsd) <- args
    # ot <- OverlayTrack(list(dtscore,dtsd))
    dtscore
  }
)

#' @rdname ModInosine-functions
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "ModSetInosine",
                        coord = "GRanges"),
  definition = function(x, coord, type = "score", window.size = 15L, ...) {
    if(missing(type)){
      type <- "score"
    }
    type <- match.arg(type, "score")
    callNextMethod(x = x, coord = coord, type = type, window.size = window.size,
                   ...)
  }
)
#' @rdname ModInosine-functions
#' @export
setMethod(
  f = "visualizeData",
  signature = signature(x = "ModSetInosine"),
  definition = function(x, name, from, to, type = "score", ...) {
    if(missing(type)){
      type <- "score"
    }
    type <- match.arg(type, "score")
    callNextMethod(x = x, name, from, to, type = type, ...)
  }
)

