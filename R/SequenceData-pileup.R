#' @include RNAmodR.R
#' @include SequenceData-class.R
NULL

#' @name PileupSequenceData
#' 
#' @title PileupSequenceData
#' 
#' @description
#' title
NULL

#' @rdname PileupSequenceData
#' @export
setClass(Class = "PileupSequenceData",
         contains = "SequenceData",
         prototype = list(minQuality = 5L))

# PileupSequenceData ----------------------------------------------------------------

.pileup_colnames <- c("pos","-","G","A","T","C")
.pileup_measure_colnames <- c("-","G","A","T","C")

#' @importFrom reshape2 dcast melt
.fill_up_pileup_data <- function(d,r){
  pos <- BiocGenerics::start(r):BiocGenerics::end(r)
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
  parentRanges <- RNAmodR:::.get_parent_annotations(ranges)
  pileupArgs <- pileupArgs[!vapply(pileupArgs,is.null,logical(1))]
  # get data per chromosome
  pileupParam <- do.call("PileupParam",pileupArgs)
  pileup <- Rsamtools::pileup(bamFile,
                              scanBamParam = param,
                              pileupParam = pileupParam)
  pileup <- S4Vectors::DataFrame(pileup)
  # reformat data
  pileup$seqnames <- as.character(pileup$seqnames)
  pileup$nucleotide <- as.character(pileup$nucleotide)
  # split into data per transcript which is defined by the which_label column
  # format: chromosome:start-end
  pileup <- split(pileup,
                  pileup$which_label)
  # sanitize pilup data
  # - keep only data for correct strand
  # - fillup empty positions with zero
  strands <- as.character(BiocGenerics::strand(parentRanges))
  rl <- split(parentRanges,seq_along(parentRanges))
  pileup <- IRanges::SplitDataFrameList(
    mapply(
      function(d,r,strand){
        ans <- NULL
        d <- d[d$strand == strand,]
        if(nrow(d) > 0) {
          ans <- reshape2::dcast(
            as.data.frame(d),
            pos ~ nucleotide,
            sum,
            value.var = "count")
        }
        ans <- .fill_up_pileup_data(ans,r)
        # remove pos column since we don't need this anymore. seq_along == pos
        ans$pos <- NULL
        ans
      },
      pileup,
      rl,
      strands,
      SIMPLIFY = FALSE))
  names(pileup) <- parentRanges$ID
  pileup
}

setMethod(".get_Data",
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

#' @name PileupSequenceData
#' @export
PileupSequenceData <- function(bamfiles, annotation, sequences, seqinfo, ...){
  SequenceData("Pileup", files = bamfiles, annotation = annotation,
               sequences = sequences, seqinfo = seqinfo, ...)
}

# aggregation ------------------------------------------------------------------

# aggregate
# - calculate percentage
# - calculate mean per observation
# - calculate sd per observation
#' @importFrom matrixStats rowSds
.aggregate_data_frame_percentage_mean_sd <- function(x,condition){
  df <- .subset_to_condition(x@unlistData,
                             x@conditions,
                             condition)
  # set up some base values
  replicates <- unique(x@replicate)
  ncol <- ncol(df[,x@replicate == 1L,drop = FALSE])
  seqAdd <- seq.int(from = 0, to = ncol(df) - 1, by = ncol)
  colNames <- strsplit(colnames(df)[seq_len(ncol)],"\\.")
  colNames <- IRanges::CharacterList(colNames)[as.list(lengths(colNames))]
  # get percentage per replicate
  for(i in seq_along(replicates)){
    df[,x@replicate == i] <- 
      as.data.frame(df[,x@replicate == i]) / 
      rowSums(as.data.frame(df[,x@replicate == i]))
  }
  # get means
  means <- IRanges::NumericList(lapply(seq_len(ncol),
                                       function(i){
                                         rowMeans(as.data.frame(df[,i + seqAdd]),
                                                  na.rm = TRUE)
                                       }))
  names(means) <- paste0("means.",colNames)
  # get sds
  sds <- IRanges::NumericList(lapply(seq_len(ncol),
                                     function(i){
                                       matrixStats::rowSds(as.matrix(df[,i + seqAdd]),
                                                           na.rm = TRUE)
                                     }))
  names(sds) <- paste0("sds.",colNames)
  # merge data
  ans <- cbind(do.call(DataFrame, means),
               do.call(DataFrame, sds))
  ans <- IRanges::SplitDataFrameList(ans)
  ans@partitioning <- x@partitioning
  ans
}

#' @name PileupSequenceData
#' @export
setMethod("aggregate",
          signature = c(x = "PileupSequenceData"),
          function(x, condition = c("Both","Treated","Control")){
            condition <- tolower(match.arg(condition))
            .aggregate_data_frame_percentage_mean_sd(x,condition)
          }
)


# data visualization -----------------------------------------------------------
setMethod(
  f = ".dataTracks",
  signature = signature(x = "PileupSequenceData",
                        data = "missing",
                        seqdata = "GRanges",
                        sequence = "XString"),
  definition = function(x, seqdata, sequence,  args) {
    requireNamespace("Gviz")
    browser()
  }
)
