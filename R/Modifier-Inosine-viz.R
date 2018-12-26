#' @include RNAmodR.R
#' @include Modifier-Inosine-class.R
NULL

RNAMODR_NUCLEOTIDE_COLOUR <- 
  c("." = biovizBase::getBioColor("RNA_BASES_N")[["N"]],
    "G" = biovizBase::getBioColor("RNA_BASES_N")[["G"]],
    "A" = biovizBase::getBioColor("RNA_BASES_N")[["A"]],
    "U" = biovizBase::getBioColor("RNA_BASES_N")[["U"]],
    "C" = biovizBase::getBioColor("RNA_BASES_N")[["C"]])

#' @rdname ModInosine
#' @export
setMethod(
  f = "visualizeData", 
  signature = signature(x = "ModInosine"),
  definition = function(x,
                        i,
                        start,
                        end){
    requireNamespace("Gviz")
    if(missing(i)){
      i <- 1L
    }
    data <- aggregateData(x)[[i]]
    # get coordinates
    coord <- .norm_viz_coord(data,start,end)
    # get plotting data
    data <- aggregateData(x)[[i]]
    data <- data[,grepl("means",
                        colnames(data))]
    colnames(data) <- c(".","G","A","U","C")
    # get plotting sequence
    seq <- sequences(x)[[i]]
    #
    r <- .get_parent_annotations(ranges(x))[i]
    chromosome <- .norm_viz_chromosome(r)
    genome <- .norm_viz_genome(r)
    # get plotting setting
    groups <- colnames(data)
    colour <- RNAMODR_NUCLEOTIDE_COLOUR[match(names(RNAMODR_NUCLEOTIDE_COLOUR),
                                              groups)]
    # create plot
    dataTrack <- DataTrack(start = seq_len(nrow(data)), 
                           width = 0.5,
                           genome = genome,
                           chromosome = rep(chromosome,nrow(data)),
                           name = "rel. nucleotide occurance",
                           data = t(as.matrix(data)),
                           type = "hist",
                           groups = factor(groups,levels = groups),
                           col =  colour)
    seqTrack <- SequenceTrack(DNAStringSet(c("chrNA" = seq)),
                              chromosome = chromosome,
                              noLetters = TRUE)
    plotTracks(list(dataTrack,seqTrack),
               from = coord$start,
               to = coord$end)
  }
)
