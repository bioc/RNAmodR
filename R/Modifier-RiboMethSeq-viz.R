#' @include RNAmodR.R
#' @include Modifier-RiboMethSeq-class.R
NULL

RNAMODR_RMS_PLOT_DATA <- c("ends",
                           "scoreA",
                           "scoreB",
                           "scoreRMS")
RNAMODR_RMS_PLOT_DATA_NAMES <- c(ends = "5'- & 3'-ends",
                                 scoreA = "Score A",
                                 scoreB = "Score B",
                                 scoreRMS = "Score RiboMethSeq")
RNAMODR_RMS_PLOT_DATA_COLOURS <- c(ends = "#FBB4AE",
                                   scoreA = "#B3CDE3",
                                   scoreB = "#CCEBC5",
                                   scoreRMS = "#DECBE4")

#' @rdname ModRiboMethSeq
#' @export
setMethod(
  f = "visualizeDataByCoord",
  signature = signature(x = "ModRiboMethSeq",
                        coord = "GRanges"),
  definition = function(x,
                        coord,
                        type = c("ends","scoreA","scoreB","scoreRMS"),
                        window.size = 15L,
                        ...) {
    if(missing(type)){
      type <- RNAMODR_RMS_PLOT_DATA
    }
    callNextMethod(x = x,
                   coord = coord,
                   type = type,
                   window.size = window.size,
                   ...)
  }
)

setMethod(
  f = ".dataTracksByCoord",
  signature = signature(x = "ModRiboMethSeq",
                        data = "GRanges"),
  definition = function(x,
                        data,
                        args) {
    requireNamespace("Gviz")
    n <- ncol(mcols(data))
    colour <- args[["colour"]]
    if(is.na(colour) || length(colour) != n){
      colour <- RNAMODR_RMS_PLOT_DATA_COLOURS
    }
    dts <- lapply(seq_len(n),
                  function(i){
                    column <- colnames(mcols(data)[i])
                    colour <- colour[column]
                    name <- RNAMODR_RMS_PLOT_DATA_NAMES[column]
                    DataTrack(data,
                              data = column,
                              name = name,
                              fill = colour,
                              type = "histogram")
                  })
    dts
  }
)
