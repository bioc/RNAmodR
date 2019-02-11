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
#' \code{ProtectedEndSequenceData} implements
#' \code{\link[=SequenceData-class]{SequenceData}} to contain and aggregate the
#' start and ends of reads per position along the transcripts. 
#' \code{ProtectedEndSequenceData} behaves specially and offsets the start 
#' position by -1 to align the information on the 5'-3'-phosphate bonds to one
#' position. The \code{ProtectedEndSequenceData} class is specifically used
#' for the RiboMethSeq method.
#' 
#' \code{aggregate} calculates the mean and sd for samples in the \code{control}
#' and \code{treated} condition separatly.
#' title
#' 
#' @param bamfiles,annotation,sequences,seqinfo,... See 
#' \code{\link[=SequenceData-class]{SequenceData}}
#' @param x a \code{ProtectedEndSequenceData}
#' @param name For \code{getDataTrack}: a valid transcript name. Must be a name
#' of \code{ranges(x)}
#' @param condition For \code{\link{aggregate}}: condition for which the data 
#' should be aggregated.
#' 
#' @return a \code{ProtectedEndSequenceData} object
#' 
#' @examples 
#' # Construct a ProtectedEndSequenceData object
#' RNAmodR.files <- RNAmodR.Data::RNAmodR.files
#' annotation <- rtracklayer::GFF3File(RNAmodR.files[["example_RMS.gff3"]])
#' sequences <- Rsamtools::FaFile(RNAmodR.files[["example_RMS.fasta"]])
#' files <- c(control = RNAmodR.files[["example_RMS1.bam"]],
#'            treated = RNAmodR.files[["example_RMS2.bam"]])
#' pesd <- ProtectedEndSequenceData(files, annotation = annotation,
#'                                  sequences = sequences)
#' # aggregate data
#' aggregate(pesd)
NULL

#' @rdname ProtectedEndSequenceData
#' @export
setClass(Class = "ProtectedEndSequenceData",
         contains = "SequenceData",
         prototype = list(minQuality = 5L))

#' @rdname ProtectedEndSequenceData
#' @export
ProtectedEndSequenceData <- function(bamfiles, annotation, sequences, seqinfo, 
                                     ...){
  SequenceData("ProtectedEnd", bamfiles = bamfiles, annotation = annotation,
               sequences = sequences, seqinfo = seqinfo, ...)
}

# ProtectedEndSequenceData -----------------------------------------------------

#' @rdname RNAmodR-internals
setMethod(".getData",
          signature = c(x = "ProtectedEndSequenceData",
                        grl = "GRangesList",
                        sequences = "XStringSet",
                        param = "ScanBamParam"),
          definition = function(x, grl, sequences, param, args){
            message("Loading protected end data from BAM files ... ",
                    appendLF = FALSE)
            data <- lapply(bamfiles(x),
                           FUN = .get_position_data_of_transcript_ends,
                           grl = grl,
                           param = param,
                           type = "protected_ends",
                           args = args)
            names(data) <- paste0("protectedend.",x@condition,".",x@replicate)
            data
          }
)


# aggregation ------------------------------------------------------------------

#' @rdname ProtectedEndSequenceData
#' @export
setMethod("aggregate",
          signature = c(x = "ProtectedEndSequenceData"),
          function(x, condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_list_data_mean_sd(x,condition)
          }
)

# data visualization -----------------------------------------------------------

RNAMODR_PLOT_SEQ_PROTEND_NAMES <- c("protend" = "mean")

#' @rdname ProtectedEndSequenceData
#' @export
setMethod(
  f = "getDataTrack",
  signature = signature(x = "ProtectedEndSequenceData"),
  definition = function(x, name, ...) {
    args <- list(...)
    # DataTrack for sequence data
    seqdata <- .get_data_for_visualization(x, name)
    # clean meta data columns
    seqdata <- .clean_mcols_end(seqdata)
    seqdata <- unlist(seqdata)
    conditions <- unique(x@condition)
    if("control" %in% conditions){
      d <- seqdata[,stringr::str_detect(colnames(mcols(seqdata)),"control")]
      colnames(mcols(d)) <- gsub(".control","",colnames(mcols(d)))
      dt.control <- Gviz::DataTrack(
        range = d,
        group = "means",
        name = paste0(RNAMODR_PLOT_SEQ_PROTEND_NAMES["protend"],
                      "\ncontrol"),
        type = "histogram")
      Gviz::displayPars(dt.control)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.control)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.control)$col.axis <- "#000000"
      Gviz::displayPars(dt.control) <- args
      track <- list("ProtectedEnd" = dt.control)
    }
    if("treated" %in% conditions){
      d <- seqdata[,stringr::str_detect(colnames(mcols(seqdata)),"treated")]
      colnames(mcols(d)) <- gsub(".treated","",colnames(mcols(d)))
      dt.treated <- Gviz::DataTrack(
        range = d,
        group = "means",
        name = paste0(RNAMODR_PLOT_SEQ_PROTEND_NAMES["protend"],
                      "\ntreated"),
        type = "histogram")
      Gviz::displayPars(dt.treated)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.treated)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.treated)$col.axis <- "#000000"
      Gviz::displayPars(dt.treated) <- args
      track <- list("ProtectedEnd" = dt.treated)
    }
    if(length(conditions) == 2L){
      track <- list("ProtectedEnd" = dt.control,
                    "ProtectedEnd" = dt.treated)
    }
    track
  }
)
