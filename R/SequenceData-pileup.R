#' @include RNAmodR.R
#' @include SequenceData-class.R
NULL

#' @name PileupSequenceData-class
#' @aliases PileupSequenceData
#' 
#' @title PileupSequenceData
#' 
#' @description
#' The \code{PileupSequenceData} aggregates the pileup of called bases per 
#' position.
#' 
#' \code{PileupSequenceData} contains five columns per data file named using the
#' following naming convention \code{pileup.condition.replicate}. The five
#' columns are distinguished by additional identifiers \code{-}, \code{G},
#' \code{A}, \code{T} and \code{C}.
#' 
#' \code{aggregate} calculates the mean and sd for each nucleotide in the
#' \code{control} and \code{treated} condition separatly. The results are then
#' normalized to a row sum of 1.
#' 
#' @param bamfiles,annotation,seqinfo,grl,sequences,param,args,... See 
#' \code{\link[=SequenceData-class]{SequenceData}} and
#' \code{\link[=SequenceData-functions]{SequenceData-functions}}
#' @param x a \code{PileupSequenceData}
#' @param name For \code{\link[=plotDataByCoord]{getDataTrack}}: a valid 
#' transcript name. Must be a name of \code{ranges(x)}
#' @param condition For \code{\link{aggregate}}: condition for which the data 
#' should be aggregated.
#' @param df,ranges,sequence,replicate inputs for creating a 
#' \code{SequenceDataFrame}. See 
#' \code{\link[=SequenceDataFrame-class]{SequenceDataFrame}}.
#' 
#' @return a \code{PileupSequenceData} object
#' 
#' @examples
#' # Construction of a PileupSequenceData object
#' library(RNAmodR.Data)
#' library(rtracklayer)
#' annotation <- GFF3File(RNAmodR.Data.example.man.gff3())
#' sequences <- RNAmodR.Data.example.man.fasta()
#' files <- c(treated = RNAmodR.Data.example.wt.1())
#' psd <- PileupSequenceData(files, annotation = annotation,
#'                           sequences = sequences)
NULL

#' @rdname PileupSequenceData-class
#' @export
setClass(Class = "PileupSequenceDataFrame",
         contains = "SequenceDFrame")
#' @rdname PileupSequenceData-class
#' @export
PileupSequenceDataFrame <- function(df, ranges, sequence, replicate, condition,
                                    bamfiles, seqinfo){
  .SequenceDataFrame("Pileup", df, ranges, sequence, replicate, condition,
                     bamfiles, seqinfo)
}
#' @rdname PileupSequenceData-class
#' @export
setClass(Class = "PileupSequenceData",
         contains = "SequenceData",
         slots = c(unlistData = "PileupSequenceDataFrame"),
         prototype = list(unlistData = PileupSequenceDataFrame(),
                          unlistType = "PileupSequenceDataFrame",
                          minQuality = 5L,
                          dataDescription = "Pileup data"))

#' @rdname PileupSequenceData-class
#' @export
PileupSequenceData <- function(bamfiles, annotation, sequences, seqinfo, ...){
  SequenceData("Pileup", bamfiles = bamfiles, annotation = annotation,
               sequences = sequences, seqinfo = seqinfo, ...)
}

setSequenceDataCoercions("Pileup")

# PileupSequenceData ----------------------------------------------------------------

.fill_up_pileup_data <- function(pileup,grl,irl){
  which_label <- .get_which_label(irl)
  unlisted_which_label <- unlist(which_label, use.names = FALSE)
  which_label <- IRanges::CharacterList(
    split(unlisted_which_label,
          factor(names(unlisted_which_label),
                 levels = names(grl))))
  if(any(names(which_label) != names(grl))){
    stop("")
  }
  which_label <- IRanges::CharacterList(
    Map(function(w,l){unname(unlist(Map(rep,w,l)))},
        which_label,
        width(grl)))
  strand <- Map(rep, .get_strand_u_GRangesList(grl), sum(width(grl)))
  pileup_add <- data.frame(pos = unlist(.seqs_rl(grl)),
                           strand = unlist(strand),
                           nucleotide = "-",
                           count = 0L,
                           which_label = unlist(which_label),
                           row.names = NULL)
  rbind(pileup_add,pileup)
}

#' @importFrom Rsamtools pileup PileupParam
#' @importFrom reshape2 dcast
.get_position_data_of_transcript_pileup <- function(bamFile, grl, sequences,
                                                    param, args = list()){
  # get user set argumenst
  pileupArgs <- args[c("max_depth", 
                       "min_base_quality", 
                       "min_mapq", 
                       "min_nucleotide_depth", 
                       "min_minor_allele_depth", 
                       "distinguish_strands", 
                       "distinguish_nucleotides", 
                       "ignore_query_Ns",
                       "left_bins", 
                       "query_bins", 
                       "cycle_bins")]
  pileupArgs <- pileupArgs[!vapply(pileupArgs,is.null,logical(1))]
  # get data per chromosome
  pileupParam <- do.call("PileupParam",pileupArgs)
  pileup <- Rsamtools::pileup(bamFile,
                              scanBamParam = param,
                              pileupParam = pileupParam)
  pileup <- pileup[,c("pos","strand","nucleotide","count","which_label")]
  pileup <- .fill_up_pileup_data(pileup,grl,Rsamtools::bamWhich(param))
  pileup <- reshape2::dcast(pileup, which_label + pos + strand ~ nucleotide,
                            sum, value.var = "count")
  cols <- c("which_label","pos","strand","-","A","C","G","T")
  pileup <- pileup[,cols]
  pileup <- S4Vectors::DataFrame(pileup)
  colnames(pileup) <- cols
  # split into data per transcript which is defined by the which_label column
  # format: chromosome:start-end
  # merge results from different exons by creating a custom PartitioningByEnd
  # object and using it to relist
  pileup <- .splitPileupAsList_transcript(pileup, grl)
  if(length(pileup) != length(grl)){
    stop("")
  }
  # keep only data for correct strand
  strands_u <- .get_strand_u_GRangesList(grl)
  pileup <- pileup[pileup[,"strand"] == strands_u]
  # sort rev on minus strand
  f <- strands_u == "-"
  if(any(f)){
    pileup[f] <- pileup[f][order(pileup[f,"pos"],decreasing = FALSE)]
  }
  pileup <- pileup[,c("-","G","A","T","C")]
  pileup
}

#' @rdname PileupSequenceData-class
#' @export
setMethod("getData",
          signature = c(x = "PileupSequenceData",
                        bamfiles = "BamFileList",
                        grl = "GRangesList",
                        sequences = "XStringSet",
                        param = "ScanBamParam"),
          definition = function(x, bamfiles, grl, sequences, param, args){
            data <- lapply(bamfiles,
                           FUN = .get_position_data_of_transcript_pileup,
                           grl = grl,
                           sequences = sequences,
                           param = param,
                           args = args)
            names(data) <- rep("pileup",length(data))
            data
          }
)

# summary ----------------------------------------------------------------------

setMethod("summary",
          signature = "PileupSequenceData",
          function(object){
            .get_summary_MultiColSequenceData(object)
          })

# aggregation ------------------------------------------------------------------

# aggregate
# - calculate percentage
# - calculate mean per observation
# - calculate sd per observation
#' @importFrom matrixStats rowSds
.aggregate_data_frame_percentage_mean_sd <- function(x,condition){
  conditions <- conditions(x)
  f <- .subset_to_condition(conditions, condition)
  df <- as(unlist(x,use.names=FALSE),"DFrame")[,f,drop=FALSE]
  conditions_u <- unique(conditions[f])
  replicates <- replicates(x)[f]
  # set up some base values
  sample_width <- length(replicates[conditions[f] == conditions_u[1] & 
                                      replicates == unique(replicates)[1]])
  colNames <- strsplit(colnames(df)[seq_len(sample_width)],"\\.")
  colNames <- IRanges::CharacterList(colNames)[as.list(lengths(colNames))]
  # get percentage per replicate
  for(con in conditions_u){
    ff <- conditions[f] == con
    for(i in unique(replicates[ff])){
      df[,ff][,replicates[ff] == i] <- 
        as.data.frame(df[,ff,drop = FALSE][,replicates[ff] == i,drop = FALSE]) / 
        rowSums(as.data.frame(df[,ff,drop = FALSE][,replicates[ff] == i,drop = FALSE]))
    }
  }
  # get means
  means <- do.call(
    c,
    lapply(conditions_u,
           function(con){
             ff <- conditions[f] == con
             ncol <- ncol(df[,ff,drop = FALSE]
                          [,replicates[ff] == unique(replicates[ff])[1],
                            drop = FALSE])
             seqAdd <- seq.int(from = 0, 
                               to = ncol(df[,ff,drop=FALSE]) - 1, 
                               by = sample_width)
             means <- IRanges::NumericList(
               lapply(seq_len(ncol),
                      function(i){
                        unname(rowMeans(as.data.frame(df[,ff,drop=FALSE][,i + seqAdd,drop=FALSE]),
                                        na.rm = TRUE))
                      }))
             names(means) <- paste0("means.", con, ".", colNames)
             means
           }))
  # get sds
  sds <- do.call(
    c,
    lapply(conditions_u,
           function(con){
             ff <- conditions[f] == con
             ncol <- ncol(df[,ff,drop = FALSE]
                          [,replicates[ff] == unique(replicates[ff])[1],
                            drop = FALSE])
             seqAdd <- seq.int(from = 0, 
                               to = ncol(df[,ff,drop=FALSE]) - 1, 
                               by = sample_width)
             means <- IRanges::NumericList(
               lapply(seq_len(ncol),
                      function(i){
                        unname(matrixStats::rowSds(as.matrix(df[,ff,drop=FALSE][,i + seqAdd,drop=FALSE]),
                                                   na.rm = TRUE))
                      }))
             names(means) <- paste0("sds.", con, ".", colNames)
             means
           }))
  # merge data
  ans <- cbind(do.call(DataFrame, means),
               do.call(DataFrame, sds))
  ans <- relist(ans, IRanges::PartitioningByEnd(x))
  positions <- .seqs_rl_strand(ranges(x))
  rownames(ans) <- IRanges::CharacterList(positions)
  ans
}

#' @rdname PileupSequenceData-class
#' @export
setMethod("aggregateData",
          signature = c(x = "PileupSequenceData"),
          function(x, condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_data_frame_percentage_mean_sd(x,condition)
          }
)


# data visualization -----------------------------------------------------------

RNAMODR_PLOT_BASES_COLOURS <- c("G" = "#D7191C",
                                "A" = "#ABD9E9",
                                "U" = "#2C7BB6",
                                "C" = "#FDAE61")
RNAMODR_PLOT_SEQ_PILEUP_NAMES <- c("bases" = "Bases called [%]")

.norm_viz_pileup_args <- function(input){
  colour.bases <- input[["colour.bases"]]
  if(is.null(colour.bases) || 
     length(colour.bases) != length(RNAMODR_PLOT_BASES_COLOURS)){
    colour.bases <- RNAMODR_PLOT_BASES_COLOURS
  } else {
    if(any(!(names(colour.bases) %in% names(RNAMODR_PLOT_BASES_COLOURS)))){
      stop("Unrecognized names for additional argument 'colour.bases'. ",
           "The names must be ",
           paste(names(RNAMODR_PLOT_BASES_COLOURS),collapse = ","),".",
           call. = FALSE)
    }
    colour.bases <- colour.bases[match(names(colour.bases),
                                       names(RNAMODR_PLOT_BASES_COLOURS))]
  }
  input <- list(colour.bases = colour.bases)
  input
}

.clean_mcols_pileup <- function(seqdata, colour.bases){
  d <- mcols(seqdata@unlistData)
  d <- d[,!grepl("\\.\\.",colnames(d)),drop=FALSE]
  d <- d[,!grepl("sds.",colnames(d)),drop=FALSE]
  conditions <- c("control","treated")
  for(con in conditions){
    f <- grepl(con,colnames(d))
    if(any(f)){
      f <- f & grepl("means",colnames(d))
      colnames(d)[f] <- gsub("means.","",colnames(d)[f])
      d[,f] <- DataFrame(as.data.frame(d[,f]) / rowSums(as.data.frame(d[,f])) * 100)
      colnames(d)[f][colnames(d[,f]) == paste0(con,".T")] <- paste0(con,".U")
      d[,f] <- d[,f][,match(colnames(d)[f], paste0(con,".",names(colour.bases)))]
    }
  }
  mcols(seqdata@unlistData) <- d
  seqdata
}

#' @rdname PileupSequenceData-class
#' @importFrom Gviz DataTrack
#' @export
setMethod(
  f = "getDataTrack",
  signature = signature(x = "PileupSequenceData"),
  definition = function(x, name, ...) {
    args <- .norm_viz_pileup_args(list(...))
    colour.bases <- args[["colour.bases"]]
    # DataTrack for sequence data
    seqdata <- .get_data_for_visualization(x, name)
    # clean meta data columns
    seqdata <- .clean_mcols_pileup(seqdata, colour.bases)
    seqdata <- unlist(seqdata)
    conditions_u <- unique(conditions(x))
    if("control" %in% conditions_u){
      d <- seqdata[,grepl("control",colnames(mcols(seqdata)))]
      colnames(mcols(d)) <- gsub("control.","",colnames(mcols(d)))
      dt.control <- Gviz::DataTrack(range = d,
                                    groups = factor(colnames(mcols(d))),
                                    name = paste0(RNAMODR_PLOT_SEQ_PILEUP_NAMES["bases"],
                                                  "\ncontrol"),
                                    col = colour.bases[order(names(colour.bases))],
                                    type = "histogram",
                                    stackedBars = TRUE)
      Gviz::displayPars(dt.control)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.control)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.control)$col.axis <- "#000000"
      Gviz::displayPars(dt.control) <- args
      track <- list("Pileup" = dt.control)
    }
    if("treated" %in% conditions_u){
      d <- seqdata[,grepl("treated",colnames(mcols(seqdata)))]
      colnames(mcols(d)) <- gsub("treated.","",colnames(mcols(d)))
      dt.treated <- Gviz::DataTrack(range = d,
                                    groups = factor(colnames(mcols(d))),
                                    name = paste0(RNAMODR_PLOT_SEQ_PILEUP_NAMES["bases"],
                                                  "\ntreated"),
                                    col = colour.bases[order(names(colour.bases))],
                                    type = "histogram",
                                    stackedBars = TRUE)
      Gviz::displayPars(dt.treated)$background.title <- "#FFFFFF"
      Gviz::displayPars(dt.treated)$fontcolor.title <- "#000000"
      Gviz::displayPars(dt.treated)$col.axis <- "#000000"
      Gviz::displayPars(dt.treated) <- args
      track <- list("Pileup" = dt.treated)
    }
    if(length(conditions_u) == 2L){
      track <- list("Pileup" = dt.control,
                    "Pileup" = dt.treated)
    }
    track
  }
)

# special funtions -------------------------------------------------------------

#' @rdname PileupSequenceData-class
#' @export
setGeneric(name = "pileupToCoverage",
           signature = "x",
           def = function(x) standardGeneric("pileupToCoverage"))

.aggregate_pile_up_to_coverage <- function(data){
  unlisted_data <- unlist(data,use.names=FALSE)
  replicates <- unique(as.integer(interaction(conditions(data),
                                              replicates(data))))
  ans  <- IRanges::IntegerList(
    lapply(seq_along(replicates),
           function(i){
             rowSums(as.data.frame(unlisted_data[,replicates == i]))
           }))
  names(ans) <- paste0("replicate.",replicates)
  ans <- do.call(S4Vectors::DataFrame,ans)
  ans <- relist(ans, data)
  rownames(ans) <- rownames(data)
  ans
}

#' @rdname PileupSequenceData-class
#' @export
setMethod("pileupToCoverage",
          signature = "PileupSequenceData",
          definition = function(x){
            .aggregate_pile_up_to_coverage(x)
          })
