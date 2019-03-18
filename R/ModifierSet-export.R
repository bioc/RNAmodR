#' @include RNAmodR.R
#' @include Modifier-class.R
#' @include ModifierSet-class.R
NULL

#' @name ModifierSet-export
#' 
#' @title Exporting data to other classes and files
#' 
#' @description 
#' \code{Modifier}, \code{SequenceData}, \code{SequenceDataSet} and 
#' \code{SequenceDataList} can be exported into several other formats.
#' 
#' \code{export.bw} exports to a bigWig file as implemented by the
#' \code{rtracklayer} package.
#' 
#' \code{SummarizedExperiment} converts the results in a \code{ModifierSet} into
#' a \code{\link[SummarizedExperiment:RangedSummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#' 
#' @param assays a \code{\link[=ModifierSet-class]{ModifierSet}} object
#' @param object a \code{\link[=Modifier-class]{Modifier}}, a 
#' \code{\link[=SequenceData-class]{SequenceData}}, 
#' \code{\link[=SequenceDataSet-class]{SequenceDataSet}} or a
#' \code{\link[=SequenceDataList-class]{SequenceDataList}} object
#' @param con See \code{\link[rtracklayer:BigWigFile]{export.bw}}
#' @param type Which column(s) of data should be exported?
#' (A wig files only supports one type of data)
#' @param ... See \code{\link[rtracklayer:BigWigFile]{export.bw}} or 
#' optional arguments for \code{SummarizedExperiment}:
#' \itemize{
#' \item{modificationsOnly} {\code{TRUE} or \code{FALSE}: Should only data on
#' positions with a modification found be included in the 
#' \code{RangedSummarizedExperiment}? (default = 
#' \code{modificationsOnly = TRUE})}
#' \item{flanking} {If \code{modificationsOnly = FALSE} how many neighboring
#' positions should be included in the \code{RangedSummarizedExperiment}?
#' (default = \code{flanking = 10L})}
#' \item{allPositions} {\code{TRUE} or \code{FALSE}: Should data on
#' *all* positions of the transcripts be included in the 
#' \code{RangedSummarizedExperiment}? This overrides the previous two options.
#' It is potentially very resource hungry. Be careful. (default = 
#' \code{allPositions = FALSE})}
#' \item{sequenceData} {\code{TRUE} or \code{FALSE}: Should data from the 
#' \code{SequenceData} objects (utilized in the \code{ModifierSet}) instead of 
#' aggregated data from the \code{ModifierSet} object be included in the 
#' \code{RangedSummarizedExperiment}? (default = 
#' \code{sequenceData = FALSE})}
#' }
NULL

.add_seqlengths_for_bigWig_export <- function(gr, x){
  seqlengths <- GenomeInfoDb::seqlengths(seqinfo(x))
  f <- vapply(seqlengths, is.na, logical(1))
  if(any(f)){
    seqlengths[f] <- lengths(gr)[f]
  }
  GenomeInfoDb::seqlengths(gr) <- seqlengths
  gr
}

.expand_ranges_per_position <- function(ranges, positions, x){
  if(!is(positions,"IntegerList")){
    positions <- as(positions,"IntegerList")
  }
  seqnames <- mapply(rep, seqnames(ranges), lengths(positions))
  gr <- GRangesList(mapply(
    function(sn, start){
      GenomicRanges::GRanges(seqnames = sn,
                             ranges = IRanges::IRanges(start = start,
                                                       width = 1L))
    },
    seqnames,
    positions))
  gr <- .add_seqlengths_for_bigWig_export(gr, x)
  gr <- unlist(gr, use.names = FALSE)
  gr
}

# export to bigWig file --------------------------------------------------------

.get_GRanges_for_bigWig_export <- function(x, type){
  ranges <- ranges(x)
  if(is(x,"Modifier")){
    data <- aggregateData(x)
  } else {
    data <- aggregate(x)
  }
  if(is(x,"SequenceDataSet")){
    data <- do.call(cbind, unname(data))
  }
  type <- .norm_score_type(type, colnames(data@unlistData), multiple = FALSE)
  ranges <- .expand_ranges_per_position(ranges, rownames(data), x)
  mcols(ranges)$score <- unlist(data)[,type]
  ranges
}

#' @rdname ModifierSet-export
#' @importFrom rtracklayer export.bw
#' @export
setMethod("export.bw","Modifier",
          function(object, con, type, ...){
            object <- .get_GRanges_for_bigWig_export(object, type)
            export.bw(object, con, ...)
          }
)

#' @rdname ModifierSet-export
#' @export
setMethod("export.bw","SequenceData",
          function(object, con, type, ...){
            object <- .get_GRanges_for_bigWig_export(object, type)
            export.bw(object, con, ...)
          }
)

#' @rdname ModifierSet-export
#' @export
setMethod("export.bw","SequenceDataSet",
          function(object, con, type, ...){
            object <- .get_GRanges_for_bigWig_export(object, type)
            export.bw(object, con, ...)
          }
)

# export to wig file -----------------------------------------------------------

.get_GRangesList_for_Wig_export <- function(x, type){
  ranges <- ranges(x)
  if(is(x,"Modifier")){
    data <- aggregateData(x)
  } else {
    data <- aggregate(x)
  }
  if(is(x,"SequenceDataSet")){
    data <- do.call(cbind, unname(data))
  }
  type <- .norm_score_type(type, colnames(data@unlistData), multiple = TRUE)
  ranges <- .expand_ranges_per_position(ranges, rownames(data), x)
  grl <- GRangesList(lapply(type,
                            function(t){
                              gr <- ranges
                              mcols(gr)$score <- unlist(data)[,t]
                              gr
                            }))
  names(grl) <- type
  grl
}

#' @rdname ModifierSet-export
#' @importFrom rtracklayer export.wig
#' @export
setMethod("export.wig","Modifier",
          function(object, con, type, ...){
            object <- .get_GRangesList_for_Wig_export(object, type)
            export.wig(object, con, ...)
          }
)

#' @rdname ModifierSet-export
#' @export
setMethod("export.wig","SequenceData",
          function(object, con, type, ...){
            object <- .get_GRangesList_for_Wig_export(object, type)
            export.wig(object, con, ...)
          }
)

#' @rdname ModifierSet-export
#' @export
setMethod("export.wig","SequenceDataSet",
          function(object, con, type, ...){
            object <- .get_GRangesList_for_Wig_export(object, type)
            export.wig(object, con, ...)
          }
)

# export to RangedSummarizedexperiment -----------------------------------------

.get_rse_args <- function(input){
  modificationsOnly <- TRUE
  flanking <- 10L
  allPositions <- FALSE
  sequenceData <- FALSE
  if(!is.null(input[["modificationsOnly"]])){
    modificationsOnly <- input[["modificationsOnly"]]
    if(!assertive::is_a_bool(modificationsOnly)){
      stop("'modificationsOnly' must be a single logical value.")
    }
  }
  if(!is.null(input[["flanking"]])){
    flanking <- input[["flanking"]]
    if(!is.integer(flanking) || 
       flanking < 0L ||
       length(flanking) != 1){
      stop("'flanking' must be a single positive integer value.")
    }
  }
  if(!is.null(input[["allPositions"]])){
    allPositions <- input[["allPositions"]]
    if(!assertive::is_a_bool(allPositions)){
      stop("'allPositions' must be a single logical value.")
    }
  }
  if(!is.null(input[["sequenceData"]])){
    sequenceData <- input[["sequenceData"]]
    if(!assertive::is_a_bool(sequenceData)){
      stop("'sequenceData' must be a single logical value.")
    }
  }
  args <- list(modificationsOnly = modificationsOnly, flanking = flanking, 
               allPositions = allPositions, sequenceData = sequenceData)
  args
}

.get_unified_modifications <- function(assays){
  mod <- modifications(assays)
  mod@unlistData$score <- NULL
  mod <- unlist(mod)
  mod <- unname(mod[!duplicated(mod)])
  mod
}

#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
.from_ModifierSet_to_RangedSummarizedExperiment <- function(assays, ...){
  input <- list(...)
  browser()
  args <- .get_rse_args(input)
  if(!args[["allPositions"]]){
    mod <- .get_unified_modifications(assays)
    if(args[["modificationsOnly"]]){
      data <- compareByCoord(assays, mod, sequenceData = args[["sequenceData"]],
                             allTypes = TRUE)
    } else {
      data <- compareByCoord(assays, mod, sequenceData = args[["sequenceData"]],
                             flanking = args[["flanking"]], allTypes = TRUE)
    }
    
    rowData <- data[,(colnames(data) %in% c("names","positions","mod"))]
    colnames(rowData) <- c("Parent","positions","mod")
    data <- data[,!(colnames(data) %in% c("names","positions","mod"))]
    f <- match(rowData$Parent,as.character(mod$Parent))
    seqnames <- seqnames(mod)[f]
    strand <- strand(mod)[f]
  } else {
    
  }
  ir <- IRanges::IRanges(start = as.integer(rowData$positions),
                         width = 1L)
  gr <- GRanges(seqnames,
                ranges = ir,
                strand = strand,
                rowData[,colnames(rowData) != "positions"])
  ans <- SummarizedExperiment(assayList, rowRanges = gr)
  return(ans)
}


#' @rdname ModifierSet-export
#' @export
setMethod("SummarizedExperiment", "ModifierSet",
          function(assays, ...) 
            .from_ModifierSet_to_RangedSummarizedExperiment(assays, ...))
