#' @include RNAmodR.R
#' @include SequenceData-class.R
#' @include SequenceData-viz.R
NULL

RNAMODR_PROT_SEQDATA_PLOT_DATA <- c("means","sds")
RNAMODR_PROT_SEQDATA_PLOT_DATA_NAMES <- c(means = "5'- & 3'-ends",
                                          sds = "s.d.")
RNAMODR_PROT_SEQDATA_PLOT_DATA_COLOURS <- c(means = "#FBB4AE",
                                            sds = "#808080")

#' @name ProtectedEndSequenceData-class
#' @aliases ProtectedEndSequenceData
#' 
#' @title ProtectedEndSequenceData
#' 
#' @description
#' \code{ProtectedEndSequenceData} implements
#' \code{\link[=SequenceData-class]{SequenceData}} to contain and aggregate the
#' start and ends of reads per position along a transcript.
#' \code{ProtectedEndSequenceData} offsets the start position by -1 to align the
#' information on the 5'-3'-phosphate bonds to one position. The
#' \code{ProtectedEndSequenceData} class is implemented specifically as required
#' for the \code{RiboMethSeq} method.
#' 
#' The objects of type \code{ProtectedEndSequenceData} contain three columns per
#' data file named using the following naming convention
#' \code{protectedend.condition.replicate}.
#' 
#' \code{aggregate} calculates the mean and sd for samples in the \code{control}
#' and \code{treated} condition separatly.
#' 
#' @param bamfiles,annotation,seqinfo,grl,sequences,param,args,... See 
#' \code{\link[=SequenceData-class]{SequenceData}} and
#' \code{\link[=SequenceData-functions]{SequenceData-functions}}
#' @param x a \code{ProtectedEndSequenceData}
#' @param name For \code{\link[=visualizeDataByCoord]{getDataTrack}}: a valid 
#' transcript name. Must be a name of \code{ranges(x)}
#' @param condition For \code{\link{aggregate}}: condition for which the data 
#' should be aggregated.
#' 
#' @return a \code{ProtectedEndSequenceData} object
#' 
#' @examples 
#' # Construction of a ProtectedEndSequenceData object
#' library(RNAmodR.Data)
#' library(rtracklayer)
#' annotation <- GFF3File(RNAmodR.Data.example.man.gff3())
#' sequences <- RNAmodR.Data.example.man.fasta()
#' files <- c(treated = RNAmodR.Data.example.bam.1())
#' pesd <- ProtectedEndSequenceData(files, annotation = annotation,
#'                                  sequences = sequences)
NULL

#' @rdname ProtectedEndSequenceData-class
#' @export
setClass(Class = "ProtectedEndSequenceData",
         contains = "SequenceData",
         prototype = list(minQuality = 5L,
                          dataDescription = "protected end data"))

#' @rdname ProtectedEndSequenceData-class
#' @export
ProtectedEndSequenceData <- function(bamfiles, annotation, sequences, seqinfo, 
                                     ...){
  SequenceData("ProtectedEnd", bamfiles = bamfiles, annotation = annotation,
               sequences = sequences, seqinfo = seqinfo, ...)
}

# ProtectedEndSequenceData -----------------------------------------------------

#' @rdname ProtectedEndSequenceData-class
#' @export
setMethod("getData",
          signature = c(x = "ProtectedEndSequenceData",
                        bamfiles = "BamFileList",
                        grl = "GRangesList",
                        sequences = "XStringSet",
                        param = "ScanBamParam"),
          definition = function(x, bamfiles, grl, sequences, param, args){
            data <- lapply(bamfiles,
                           FUN = .get_position_data_of_transcript_ends,
                           grl = grl,
                           param = param,
                           type = "protected_ends",
                           args = args)
            names(data) <- rep("protectedend",length(data))
            data
          }
)


# aggregation ------------------------------------------------------------------

#' @rdname ProtectedEndSequenceData-class
#' @export
setMethod("aggregateData",
          signature = c(x = "ProtectedEndSequenceData"),
          function(x, condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_list_data_mean_sd(x,condition)
          }
)

# data visualization -----------------------------------------------------------

RNAMODR_PLOT_SEQ_PROTEND_NAMES <- c("protend" = "mean")

#' @rdname ProtectedEndSequenceData-class
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
