#' @include RNAmodR.R
#' @include SequenceData-class.R
NULL

#' @name PileupSequenceData
#' 
#' @title PileupSequenceData
#' 
#' @description
#' title
#' 
#' @param bamfiles,annotation,sequences,seqinfo,... See 
#' \code{\link[=SequenceData-class]{SequenceData}}
#' @param x a \code{CoverageSequenceData}
#' @param condition For \code{\link{aggregate}}: condition for which the data 
#' should be aggregated.
#' 
NULL

#' @rdname PileupSequenceData
#' @export
setClass(Class = "PileupSequenceData",
         contains = "SequenceData",
         prototype = list(minQuality = 5L))

#' @rdname PileupSequenceData
#' @export
PileupSequenceData <- function(bamfiles, annotation, sequences, seqinfo, ...){
  SequenceData("Pileup", bamfiles = bamfiles, annotation = annotation,
               sequences = sequences, seqinfo = seqinfo, ...)
}

# PileupSequenceData ----------------------------------------------------------------

.pileup_colnames <- c("pos","-","G","A","T","C")
.pileup_measure_colnames <- c("-","G","A","T","C")

#' @importFrom reshape2 dcast melt
.fill_up_pileup_data <- function(d,seq){
  pos <- seq
  if(is.null(d)){
    colnames <- .pileup_colnames
    d <- data.frame(pos = pos,
                    "-" = 0,
                    "G" = 0,
                    "A" = 0,
                    "T" = 0,
                    "C" = 0)
  } else {
    missingPos <- pos[!(pos %in% d$pos)]
    missingCols <- .pileup_colnames[!(.pileup_colnames %in% colnames(d))]
    d <- reshape2::melt(d,
                        id.vars = "pos",
                        measure.vars = .pileup_measure_colnames[
                          .pileup_measure_colnames %in% colnames(d)])
    if(length(missingPos) > 0){
      d <- rbind(d,
                 data.frame(pos = missingPos,
                            variable = "-",
                            value = 0))
    }
    mpos <- max(d$pos)
    if(length(missingCols) > 0){
      d <- do.call(
        rbind,
        c(list(d),
          lapply(
            missingCols,
            function(c){
              data.frame(pos = mpos,
                         variable = c,
                         value = ifelse(is.null(d[d$pos == mpos,"value"]),
                                        0,
                                        d[d$pos == mpos,"value"]))
            })))
    }
    d <- reshape2::dcast(d, pos ~ variable, fun.aggregate = sum, fill = 0)
    colnames <- colnames(d)
  }
  df <- S4Vectors::DataFrame(d)
  colnames(df) <- colnames
  df <- df[,.pileup_colnames]
  df
}

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
                       "include_deletions", 
                       "include_insertions", 
                       "left_bins", 
                       "query_bins", 
                       "cycle_bins")]
  pileupArgs <- pileupArgs[!vapply(pileupArgs,is.null,logical(1))]
  # get data per chromosome
  pileupParam <- do.call("PileupParam",pileupArgs)
  pileup <- Rsamtools::pileup(bamFile,
                              scanBamParam = param,
                              pileupParam = pileupParam)
  pileup <- S4Vectors::DataFrame(pileup)
  # split into data per transcript which is defined by the which_label column
  # format: chromosome:start-end
  pileup <- split(pileup, pileup$which_label)
  if(length(pileup) != length(grl)){
    stop("Something went wrong.")
  }
  # sanitize pilup data
  # - keep only data for correct strand
  # - fillup empty positions with zero
  strands_u <- .get_strand_u_GRangesList(grl)
  seqs <- .seqs_rl(grl)
  pileup <- IRanges::SplitDataFrameList(
    mapply(
      function(d,seq,strand){
        ans <- NULL
        d <- d[d$strand == strand,]
        if(nrow(d) > 0) {
          ans <- reshape2::dcast(as.data.frame(d), pos ~ nucleotide, sum,
                                 value.var = "count")
        }
        ans <- .fill_up_pileup_data(ans,seq)
        # remove pos column since we don't need this anymore. seq_along == pos
        ans$pos <- NULL
        ans
      },
      pileup,
      seqs,
      strands_u,
      SIMPLIFY = FALSE))
  names(pileup) <- names(grl)
  pileup
}

#' @rdname RNAmodR-internals
setMethod(".getData",
          signature = c(x = "PileupSequenceData",
                        grl = "GRangesList",
                        sequences = "XStringSet",
                        param = "ScanBamParam"),
          definition = function(x, grl, sequences, param, args){
            message("Loading Pileup data from BAM files ... ",
                    appendLF = FALSE)
            files <- bamfiles(x)
            data <- lapply(files,
                           FUN = .get_position_data_of_transcript_pileup,
                           grl = grl,
                           sequences = sequences,
                           param = param,
                           args = args)
            names(data) <- paste0("pileup.",
                                  names(files),
                                  ".",
                                  seq_along(files))
            data
          }
)


# aggregation ------------------------------------------------------------------

# aggregate
# - calculate percentage
# - calculate mean per observation
# - calculate sd per observation
#' @importFrom matrixStats rowSds
.aggregate_data_frame_percentage_mean_sd <- function(x,condition){
  f <- .subset_to_condition(x@condition, condition)
  df <- x@unlistData[f]
  conditions <- unique(x@condition[f])
  replicates <- x@replicate[f]
  # set up some base values
  sample_width <- length(replicates[x@condition[f] == conditions[1] & 
                                      replicates == unique(replicates)[1]])
  colNames <- strsplit(colnames(df)[seq_len(sample_width)],"\\.")
  colNames <- IRanges::CharacterList(colNames)[as.list(lengths(colNames))]
  # get percentage per replicate
  for(con in conditions){
    ff <- x@condition[f] == con
    for(i in unique(replicates[ff])){
      df[,ff][,replicates[ff] == i] <- 
        as.data.frame(df[,ff,drop = FALSE][,replicates[ff] == i,drop = FALSE]) / 
        rowSums(as.data.frame(df[,ff,drop = FALSE][,replicates[ff] == i,drop = FALSE]))
    }
  }
  # get means
  means <- do.call(
    c,
    lapply(conditions,
           function(con){
             ff <- x@condition[f] == con
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
    lapply(conditions,
           function(con){
             ff <- x@condition[f] == con
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
  ans <- IRanges::SplitDataFrameList(ans)
  ans@partitioning <- x@partitioning
  ans
}

#' @rdname PileupSequenceData
#' @export
setMethod("aggregate",
          signature = c(x = "PileupSequenceData"),
          function(x, condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_data_frame_percentage_mean_sd(x,condition)
          }
)


# data visualization -----------------------------------------------------------

RNAMODR_PLOT_BASES_COLOURS <- 
  c("G" = biovizBase::getBioColor("RNA_BASES_N")[["G"]],
    "A" = biovizBase::getBioColor("RNA_BASES_N")[["A"]],
    "U" = biovizBase::getBioColor("RNA_BASES_N")[["U"]],
    "C" = biovizBase::getBioColor("RNA_BASES_N")[["C"]])
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
  d <- d[,!stringr::str_detect(colnames(d),"\\.\\."),drop=FALSE]
  d <- d[,!stringr::str_detect(colnames(d),"sds."),drop=FALSE]
  conditions <- c("control","treated")
  for(con in conditions){
    f <- stringr::str_detect(colnames(d),con)
    if(any(f)){
      f <- f & stringr::str_detect(colnames(d),"means")
      colnames(d)[f] <- gsub("means.","",colnames(d)[f])
      d[,f] <- DataFrame(as.data.frame(d[,f]) / rowSums(as.data.frame(d[,f])) * 100)
      colnames(d)[f][colnames(d[,f]) == paste0(con,".T")] <- paste0(con,".U")
      d[,f] <- d[,f][,match(colnames(d)[f], paste0(con,".",names(colour.bases)))]
    }
  }
  mcols(seqdata@unlistData) <- d
  seqdata
}

#' @rdname PileupSequenceData
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
    conditions <- unique(x@condition)
    if("control" %in% conditions){
      d <- seqdata[,stringr::str_detect(colnames(mcols(seqdata)),"control")]
      colnames(mcols(d)) <- gsub("control.","",colnames(mcols(d)))
      dt.control <- Gviz::DataTrack(range = d,
                                    groups = colnames(mcols(d)),
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
    if("treated" %in% conditions){
      d <- seqdata[,stringr::str_detect(colnames(mcols(seqdata)),"treated")]
      colnames(mcols(d)) <- gsub("treated.","",colnames(mcols(d)))
      dt.treated <- Gviz::DataTrack(range = d,
                                    groups = colnames(mcols(d)),
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
    if(length(conditions) == 2L){
      track <- list("Pileup" = dt.control,
                    "Pileup" = dt.treated)
    }
    track
  }
)
