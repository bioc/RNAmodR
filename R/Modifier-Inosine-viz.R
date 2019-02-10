#' @include RNAmodR.R
#' @include Modifier-Inosine-class.R
NULL

RNAMODR_I_PLOT_DATA <- c("score")
RNAMODR_I_PLOT_DATA_DEFAULT <- c("score")

RNAMODR_I_PLOT_DATA_COLOURS <- c("score" = "#ABABAB") 
RNAMODR_I_PLOT_DATA_NAMES <- c(score = "Score Inosine")

.norm_viz_mod_inosine_args <- function(input, type){
  if(!all(type %in% RNAMODR_I_PLOT_DATA)){
    stop("Type '",type,"' is not valid. Valid types are: '",
         paste0(RNAMODR_I_PLOT_DATA, collapse = "','"),"'.",
         call. = FALSE)
  }
  colour <- input[["colour"]]
  if(!is.null(input[["colour"]])){
    colour <- .norm_viz_colour(input[["colour"]], type)
  } else {
    colour <- RNAMODR_I_PLOT_DATA_COLOURS[type]
  }
  input <- list(type = type,
                colour = colour)
  input
}

.clean_mcols_mod_inosine <- function(x, data){
  d <- mcols(data@unlistData)
  seq <- sequences(x)[names(data)]
  letters <- unlist(unname(strsplit(as.character(seq),"")))
  f <- d$score < 0
  d$sd[f] <- 0
  d$score[f] <- 0
  d <- d[colnames(d) %in% c("score","sd")]
  d$score[letters != "A"] <- 0
  mcols(data@unlistData) <- d
  data
}

#' @rdname ModInosine-functions
#' @export
setMethod(
  f = "getDataTrack",
  signature = signature(x = "ModInosine"),
  definition = function(x, name, type, ...) {
    args <- .norm_viz_mod_inosine_args(list(...), type)
    data <- .get_data_for_visualization(x, name)
    data <- .clean_mcols_mod_inosine(x, data)
    data <- unlist(data)
    lim <- c(min(mcols(data)$score), max(mcols(data)$score))
    dtscore <- Gviz::DataTrack(range = data[,"score"],
                               groups = "score",
                               name = RNAMODR_I_PLOT_DATA_NAMES["score"],
                               col = args[["colour"]]["score"],
                               type = "histogram",
                               ylim = lim)
    Gviz::displayPars(dtscore)$background.title <- "#FFFFFF"
    Gviz::displayPars(dtscore)$fontcolor.title <- "#000000"
    Gviz::displayPars(dtscore)$col.axis <- "#000000"
    Gviz::displayPars(dtscore) <- args[names(args) != "type"]
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
    list("score" = dtscore)
  }
)

#' @rdname ModInosine-functions
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "ModInosine",
                        coord = "GRanges"),
  definition = function(x, coord, type = "score", window.size = 15L, ...) {
    if(missing(type)){
      type <- RNAMODR_I_PLOT_DATA_DEFAULT
    }
    type <- match.arg(type, RNAMODR_I_PLOT_DATA)
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
    if(missing(type)){
      type <- RNAMODR_I_PLOT_DATA_DEFAULT
    }
    type <- match.arg(type, RNAMODR_I_PLOT_DATA)
    callNextMethod(x = x, name, from, to, type = type, ...)
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
      type <- RNAMODR_I_PLOT_DATA_DEFAULT
    }
    type <- match.arg(type, RNAMODR_I_PLOT_DATA)
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
      type <- RNAMODR_I_PLOT_DATA_DEFAULT
    }
    type <- match.arg(type, RNAMODR_I_PLOT_DATA)
    callNextMethod(x = x, name, from, to, type = type, ...)
  }
)

