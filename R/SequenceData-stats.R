#' @include SequenceData-class.R
#' @include Modifier-class.R
#' @include ModifierSet-class.R
NULL

#' @name stats
#' 
#' @title Retrieving information about used reads in RNAmodR
#' 
#' @description 
#' \code{stats} returns information about reads used in the RNAmodR analysis.
#' Three modes are available depending on which type of object is provided. If a
#' \code{\link[=SequenceData-class]{SequenceData}} object is provided, a
#' \code{\link[Rsamtools:BamFile-class]{BamFile}} or
#' \code{\link[Rsamtools:BamFile-class]{BamFileList}} must be provided as well. If a
#' \code{\link[=Modifier-class]{Modifier}} object is used, the bam files
#' returned from the \code{bamfiles} function are used. This is also the case,
#' if a \code{\link[=ModifierSet-class]{ModifierSet}} object is used.
#' 
#' @param x a \code{\link[=SequenceData-class]{SequenceData}}, 
#'   \code{\link[=Modifier-class]{Modifier}} or 
#'   \code{\link[=ModifierSet-class]{ModifierSet}} object
#' @param file a \code{\link[Rsamtools:BamFile-class]{BamFile}} or 
#'   \code{\link[Rsamtools:BamFile-class]{BamFileList}}, if \code{x} is a 
#'   \code{\link[=SequenceData-class]{SequenceData}} object.
#' @param ... optional parameters used as stated 
#'   \code{\link[=SequenceData-class]{here}} (except \code{minQuality}),
#'   if \code{x} is a \code{\link[=SequenceData-class]{SequenceData}} object.
#'   
#' @return a \code{DataFrame}, \code{DataFrameList} or \code{SimpleList} with 
#'   the results in aggregated form
#'   
#' @examples 
#' library(RNAmodR.Data)
#' library(rtracklayer)
#' sequences <- RNAmodR.Data.example.AAS.fasta()
#' annotation <- GFF3File(RNAmodR.Data.example.AAS.gff3())
#' files <- list("SampleSet1" = c(treated = RNAmodR.Data.example.wt.1(),
#'                                treated = RNAmodR.Data.example.wt.2(),
#'                                treated = RNAmodR.Data.example.wt.3()),
#'               "SampleSet2" = c(treated = RNAmodR.Data.example.bud23.1(),
#'                                treated = RNAmodR.Data.example.bud23.2()),
#'               "SampleSet3" = c(treated = RNAmodR.Data.example.trm8.1(),
#'                                treated = RNAmodR.Data.example.trm8.2()))
#' msi <- ModSetInosine(files, annotation = annotation, sequences = sequences)
#' # smallest chunk of information
#' stats(sequenceData(msi[[1L]]),bamfiles(msi[[1L]])[[1L]])
#' # partial information
#' stats(sequenceData(msi[[1L]]),bamfiles(msi[[1L]]))
#' # the whole stats
#' stats(msi)
NULL

#' @importFrom Rsamtools idxstatsBam
.get_read_stats <- function(bamFile, data, hits, grl){
  # get basic stats per seqname
  idxstats <- S4Vectors::DataFrame(Rsamtools::idxstatsBam(bamFile))
  init <- IRanges::IntegerList(as.list(rep(NA_integer_,nrow(idxstats))))
  names(init) <- idxstats$seqnames 
  idxstats$used <- init
  idxstats$used_distro <- as(init,"SimpleList")
  # add number of reads per range and per seqname
  used <- lengths(split(data[queryHits(hits)],subjectHits(hits)))
  used <- split(used, seqnames(unlist(grl, use.names = FALSE)))
  idxstats$used[names(used)] <- used 
  # add length distribution of reads per range and per seqname
  distro <- GenomicAlignments::qwidth(data)
  distro <- split(distro[queryHits(hits)],subjectHits(hits))
  distro <- IRanges::IntegerList(lapply(distro,table))
  distro <- split(distro, seqnames(unlist(grl, use.names = FALSE)))
  idxstats$used_distro[names(distro)] <- distro 
  # convert distro to SimpleList
  idxstats$used_distro <- as(idxstats$used_distro,"SimpleList")
  idxstats
}

#' @rdname stats
#' @export
setMethod(f = "stats", 
          signature = signature(x = "SequenceData", file = "BamFile"),
          definition = function(x, file, ...){
            args <- .get_SequenceData_args(list(...))
            param <- .assemble_scanBamParam(ranges(x), x@minQuality, seqinfo(x))
            # get data
            data <- .load_bam_alignment_data(file, param, args)
            # get hits
            hits <- GenomicAlignments::findOverlaps(data, ranges(x))
            # get stats
            stats <- .get_read_stats(file, data, hits, ranges(x))
            stats
          })

#' @rdname stats
#' @export
setMethod(f = "stats", 
          signature = signature(x = "SequenceData", file = "BamFileList"),
          definition = function(x, file, ...){
            FUN <- function(file, x){
              stats(x, file, ...)
            }
            IRanges::DataFrameList(lapply(file,FUN,x = x))
          })

#' @rdname stats
#' @export
setMethod(f = "stats", 
          signature = signature(x = "Modifier", file = "missing"),
          definition = function(x){
            sequenceData <- sequenceData(x)
            if(is(sequenceData,"SequenceDataSet")){
              sequenceData <- sequenceData[[1L]]
            }
            if(is(sequenceData,"SequenceDataList")){
              sequenceData <- sequenceData[[1L]]
            }
            do.call("stats", c(list(sequenceData,bamfiles(x)),
                               settings(x)))
          })

#' @rdname stats
#' @export
setMethod(f = "stats", 
          signature = signature(x = "ModifierSet", file = "missing"),
          definition = function(x){
            S4Vectors::SimpleList(lapply(x,stats))
          })
