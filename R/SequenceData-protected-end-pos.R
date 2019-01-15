#' @include RNAmodR.R
#' @include SequenceData-class.R
#' @include SequenceData-viz.R
NULL

RNAMODR_PROT_SEQDATA_PLOT_DATA <- c("means","sds")
RNAMODR_PROT_SEQDATA_PLOT_DATA_NAMES <- c(means = "5'- & 3'-ends",
                                          sds = "s.d.")
RNAMODR_PROT_SEQDATA_PLOT_DATA_COLOURS <- c(means = "#FBB4AE",
                                            sds = "#808080")

#' @name ProtectedEndSequenceData
#' 
#' @title ProtectedEndSequenceData
#' 
#' @description
#' title
NULL

#' @rdname ProtectedEndSequenceData
#' @export
setClass(Class = "ProtectedEndSequenceData",
         contains = "SequenceData",
         prototype = list(minQuality = 5L))


# ProtectedEndSequenceData -----------------------------------------------------

#' @rdname ProtectedEndSequenceData
#' @export
ProtectedEndSequenceData <- function(bamfiles,
                                     fasta,
                                     gff,
                                     ...){
  browser()
  args <- .get_mod_data_args(...)
  ans <- new("ProtectedEndSequenceData",
             bamfiles,
             fasta,
             gff,
             args)
  ranges <- .load_annotation(ans@gff)
  sequences <- .load_transcript_sequences(ans@fasta,
                                          ranges)
  param <- .assemble_scanBamParam(ranges,
                                  ans@minQuality,
                                  ans@chromosomes)
  message("Loading protected end data from BAM files...")
  data <- lapply(ans@bamfiles,
                 FUN = .get_position_data_of_transcript_ends,
                 ranges = ranges,
                 param = param,
                 type = "protected_ends",
                 args = args)
  names(data) <- paste0("protectedend.",
                        names(ans@bamfiles),
                        ".",
                        seq_along(ans@bamfiles))
  .postprocess_read_data(ans,
                         data,
                         ranges,
                         sequences)
}

#' @name EndSequenceData
#' @export
setMethod("aggregate",
          signature = c(x = "ProtectedEndSequenceData"),
          function(x,
                   condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_end_data_mean_sd(x,condition)
          }
)

.fix_na_mcols <- function(seqdata){
  mcols(seqdata) <- DataFrame(lapply(mcols(seqdata),
                                     function(l){
                                       l[is.na(l)] <- 0
                                       l
                                     }))
  seqdata
}

setMethod(
  f = ".dataTracks",
  signature = signature(x = "ProtectedEndSequenceData",
                        data = "missing",
                        seqdata = "GRanges",
                        sequence = "XString"),
  definition = function(x,
                        seqdata,
                        sequence,
                        args) {
    requireNamespace("Gviz")
    seqdata <- .fix_na_mcols(seqdata)
    n <- ncol(mcols(seqdata)) / 2L
    names <- RNAMODR_PROT_SEQDATA_PLOT_DATA_NAMES
    colours <- RNAMODR_PROT_SEQDATA_PLOT_DATA_COLOURS
    dts <- lapply(seq_len(n),
                  function(i){
                    i <- i * 2L - 1
                    d <- seqdata[,seq.int(from = i, to = i + 1L)]
                    colnames(mcols(d)) <- RNAMODR_PROT_SEQDATA_PLOT_DATA
                    dtmeans <- DataTrack(d,
                                         data = "means",
                                         name = names["means"],
                                         fill = colours["means"],
                                         type = "histogram")
                    displayPars(dtmeans)$background.title <- "#FFFFFF"
                    displayPars(dtmeans)$fontcolor.title <- "#000000"
                    displayPars(dtmeans)$col.axis <- "#000000"
                    # mcols(d)$sdm <- mcols(d)$means - mcols(d)$sds
                    # mcols(d)$sdp <- mcols(d)$means + mcols(d)$sds
                    # dtsds <- DataTrack(d,
                    #                    data = c("sdm","sdp"),
                    #                    name = names["sds"],
                    #                    fill = colours["sds"],
                    #                    type = "confint")
                    # displayPars(dtsds)$background.title <- "#FFFFFF"
                    # displayPars(dtsds)$fontcolor.title <- "#000000"
                    # displayPars(dtsds)$col.axis <- "#000000"
                    # ot <- OverlayTrack(trackList = list(dtmeans,dtsds))
                    # displayPars(ot)$background.title <- "#FFFFFF"
                    # displayPars(ot)$fontcolor.title <- "#000000"
                    # displayPars(ot)$col.axis <- "#000000"
                    # ot
                    dtmeans 
                  })
    dts
  }
)
