#' @include RNAmodR.R
#' @include Modifier-Inosine-class.R
NULL

RNAMODR_I_PLOT_BASES_COLOURS <- 
  c("G" = biovizBase::getBioColor("RNA_BASES_N")[["G"]],
    "A" = biovizBase::getBioColor("RNA_BASES_N")[["A"]],
    "U" = biovizBase::getBioColor("RNA_BASES_N")[["U"]],
    "C" = biovizBase::getBioColor("RNA_BASES_N")[["C"]])
RNAMODR_I_PLOT_DATA_COLOURS <- c("score" = "#ABABAB") 
RNAMODR_I_PLOT_DATA_NAMES <- c(score = "Score Inosine",
                               bases = "Bases called [%]")

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

#' @name ModInosine-internals
setMethod(
  f = ".dataTracks",
  signature = signature(x = "ModInosine", data = "GRanges", seqdata = "GRanges",
                        sequence = "XString"),
  definition = function(x, data, seqdata, sequence, args) {
    requireNamespace("Gviz")
    n <- ncol(mcols(data))
    colour.bases <- args[["colour.bases"]]
    if(is.na(colour.bases) || 
       length(colour.bases) != length(RNAMODR_I_PLOT_BASES_COLOURS)){
      colour.bases <- RNAMODR_I_PLOT_BASES_COLOURS
    } else {
      if(any(!(names(colour.bases) %in% names(RNAMODR_I_PLOT_BASES_COLOURS)))){
        stop("Unrecognized names for additional argument 'colour.bases'. ",
             "The names must be ",
             paste(names(RNAMODR_I_PLOT_BASES_COLOURS),collapse = ","),".",
             call. = FALSE)
      }
      colour.bases <- colour.bases[match(names(colour.bases),
                                         names(RNAMODR_I_PLOT_BASES_COLOURS))]
    }
    colour <- args[["colour"]]
    if(is.na(colour) || length(colour) != n){
      colour <- RNAMODR_I_PLOT_DATA_COLOURS
    }
    # DataTrack for sequence data
    mcols(seqdata) <- 
      mcols(seqdata)[,stringr::str_detect(colnames(mcols(seqdata)),"means")]
    colnames(mcols(seqdata)) <- gsub("means.treated.","",colnames(mcols(seqdata)))
    mcols(seqdata) <- 
      DataFrame(as.data.frame(mcols(seqdata)[,colnames(mcols(seqdata)) != "."]) * 100)
    colnames(mcols(seqdata))[colnames(mcols(seqdata)) == "T"] <- "U"
    mcols(seqdata) <- mcols(seqdata)[,match(colnames(mcols(seqdata)),
                                            names(colour.bases))]
    dtbases <- DataTrack(seqdata,
                         groups = colnames(mcols(seqdata)),
                         name = RNAMODR_I_PLOT_DATA_NAMES["bases"],
                         col = colour.bases[order(names(colour.bases))],
                         type = "histogram",
                         stackedBars = TRUE)
    displayPars(dtbases)$background.title <- "#FFFFFF"
    displayPars(dtbases)$fontcolor.title <- "#000000"
    displayPars(dtbases)$col.axis <- "#000000"
    # DataTrack for score
    lim <- c(0,max(mcols(data)$score))
    mcols(data)$score[mcols(data)$score < 0] <- 0
    mcols(data)$score[strsplit(as.character(sequence),"")[[1]] == "G"] <- 0
    dtscore <- DataTrack(data,
                         data = "score",
                         name = RNAMODR_I_PLOT_DATA_NAMES["score"],
                         fill = colour,
                         type = "histogram",
                         ylim = lim)
    displayPars(dtscore)$background.title <- "#FFFFFF"
    displayPars(dtscore)$fontcolor.title <- "#000000"
    displayPars(dtscore)$col.axis <- "#000000"
    # Return as a list
    list(score = dtscore,
         seqdata = dtbases)
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

