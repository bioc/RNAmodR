#' @include RNAmodR.R
#' @include Modifier-RiboMethSeq-class.R
NULL

RNAMODR_RMS_PLOT_DATA <- c("ends",
                           "scoreA",
                           "scoreB",
                           "scoreRMS")
RNAMODR_RMS_PLOT_DATA_NAMES <- c("5'- & 3'-ends",
                                 "Score A",
                                 "Score B",
                                 "Score RiboMethSeq")

#' @rdname ModRiboMethSeq
#' @export
setMethod(
  f = "visualizeData", 
  signature = signature(x = "ModRiboMethSeq"),
  definition = function(x,
                        i,
                        type = c("ends","scoreA","scoreB","scoreRMS")) {
    requireNamespace("Gviz")
    type <- type[type %in% RNAMODR_RMS_PLOT_DATA]
    if(missing(i)){
      i <- 1L
    }
    browser()
    data <- aggregateData(x)[[i]]
    r <- .get_parent_annotations(ranges(x))
    data <- .norm_positions(data,start,end)
    
    dtl <- lapply(seq_along(type),
                  function(z){
                    DataTrack(start = seq_len(nrow(data)), 
                              width = 0.5,
                              genome = genome(r)[i],
                              chromosome = "chrNA",
                              name = RNAMODR_RMS_PLOT_DATA_NAMES[z], 
                              data = data[,z]) 
                  })
    plotTracks(dtl, type = "h")
  }
)
