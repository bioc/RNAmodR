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
                        type = c("ends","scoreA","scoreB","scoreRMS"),
                        start,
                        end) {
    requireNamespace("Gviz")
    type <- type[type %in% RNAMODR_RMS_PLOT_DATA]
    if(missing(i)){
      i <- 1L
    }
    # get plotting data
    data <- aggregateData(x)[[i]]
    # get coordinates
    coord <- .norm_viz_coord(data,start,end)
    # get plotting sequence
    seq <- sequences(x)[[i]]
    #
    r <- .get_parent_annotations(ranges(x))[i]
    chromosome <- .norm_viz_chromosome(r)
    genome <- .norm_viz_genome(r)
    #
    dtl <- lapply(seq_along(type),
                  function(z){
                    DataTrack(start = seq_len(nrow(data)), 
                              end = seq_len(nrow(data)), 
                              genome = genome,
                              chromosome = rep(chromosome,nrow(data)),
                              name = RNAMODR_RMS_PLOT_DATA_NAMES[z], 
                              data = data[,z],
                              type = "h") 
                  })
    seqTrack <- SequenceTrack(DNAStringSet(c("chrNA" = seq)),
                              chromosome = chromosome,
                              noLetters = TRUE)
    plotTracks(c(dtl,
                 list(seqTrack)),
               from = coord$start,
               to = coord$end)
  }
)
