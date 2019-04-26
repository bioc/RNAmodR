#' @include RNAmodR.R
#' @include Modifier-class.R
#' @include ModifierSet-class.R
NULL

# documentation ----------------------------------------------------------------

#' @name ModInosine
#' @aliases Inosine ModifierInosine ICE-Seq
#' 
#' @author Felix G.M. Ernst [aut]
#' 
#' @title ModInosine
#' 
#' @description 
#' Inosine can be detected in RNA-Seq data by the conversion of A positions to 
#' G. This conversion is detected by \code{ModInosine} and used to search for 
#' Inosine positions. \code{dataType} is \code{"PileupSequenceData"}.
#' 
#' Only samples named \code{treated} are used for this analysis, since the
#' A to G conversion is common feature among the reverse transcriptases usually
#' emploied. Let us know, if that is not the case, and the class needs to be
#' modified.
#' 
#' Further information on \code{\link[=ModInosine-functions]{Functions}} of 
#' \code{ModInosine}.
#' 
#' @details
#' \code{ModInosine} score: the scores for reported Inosine positions are 
#' between 0 and 1. They are calculated as the relative amount of called G bases 
#' (\code{(G / N)}) and only kept for A positions.
#' 
#' @param x the input which can be of the different types depending on whether
#' a \code{ModRiboMethSeq} or a \code{ModSetRiboMethSeq} object is to be 
#' constructed. For more information have a look at the documentation of
#' the \code{\link[RNAmodR:Modifier-class]{Modifier}} and 
#' \code{\link[RNAmodR:ModifierSet-class]{ModifierSet}} classes.
#' @param annotation annotation data, which must match the information contained
#' in the BAM files. This is parameter is only required if \code{x} if not a 
#' \code{Modifier} object.
#' @param sequences sequences matching the target sequences the reads were 
#' mapped onto. This must match the information contained in the BAM files. This
#' is parameter is only required if \code{x} if not a \code{Modifier} object.
#' @param seqinfo An optional \code{\link[GenomeInfoDb:Seqinfo-class]{Seqinfo}} 
#' argument or character vector, which can be coerced to one, to subset the 
#' sequences to be analyzed on a per chromosome basis.
#' @param ... Optional arguments overwriting default values, which are
#' \itemize{
#' \item{minCoverage:} {The minimal coverage at the position as integer value 
#' (default: \code{minCoverage = 10L}).}
#' \item{minReplicate:} {minimum number of replicates needed for the analysis 
#' (default: \code{minReplicate = 1L}).}
#' \item{minScore:} {minimum score to identify Inosine positions de novo 
#' (default: \code{minScore = 0.4}).}
#' }
#' 
#' @return a \code{ModInosine} or \code{ModSetInosine} object
#' 
#' @examples
#' # construction of ModInosine object
#' library(rtracklayer)
#' sequences <- RNAmodR.Data.example.AAS.fasta()
#' annotation <- GFF3File(RNAmodR.Data.example.AAS.gff3())
#' files <- c(treated = RNAmodR.Data.example.wt.1(),
#'            treated = RNAmodR.Data.example.wt.2(),
#'            treated = RNAmodR.Data.example.wt.3())
#' mi <- ModInosine(files,annotation = annotation ,sequences = sequences)
#' # construction of ModSetInosine object
#' files <- list("SampleSet1" = c(treated = RNAmodR.Data.example.wt.1(),
#'                                treated = RNAmodR.Data.example.wt.2(),
#'                                treated = RNAmodR.Data.example.wt.3()),
#'               "SampleSet2" = c(treated = RNAmodR.Data.example.bud23.1(),
#'                                treated = RNAmodR.Data.example.bud23.2()),
#'               "SampleSet3" = c(treated = RNAmodR.Data.example.trm8.1(),
#'                                treated = RNAmodR.Data.example.trm8.2()))
#' msi <- ModSetInosine(files, annotation = annotation, sequences = sequences)
NULL

#' @name ModInosine-internals
#' @aliases .dataTracks,ModInosine,GRanges,GRanges,XString-method
#' 
#' @title ModInosine internal functions
#' 
#' @description
#' These functions are not intended for general use, but are used for 
#' additional package development.
#' 
#' @param x,data,seqdata,sequence,args internally used arguments
NULL

#' @name ModInosine-functions
#' 
#' @title Functions for ModInosine
#' 
#' @description
#' All of the functions of \code{\link[RNAmodR:Modifier-class]{Modifier}} and
#' the \code{\link[RNAmodR:ModifierSet-class]{ModifierSet}} classes are 
#' inherited by the \code{ModInosine} and \code{ModSetInosine} classes.
#' 
#' Check below for the specifically implemented functions.
#' 
#' @param x a \code{\link[RNAmodR:Modifier-class]{Modifier}} or a 
#' \code{\link[RNAmodR:ModifierSet-class]{ModifierSet}} object. For more details
#' see also the man pages for the functions mentioned below.
#' @param value See \code{\link[RNAmodR:Modifier-functions]{settings}}
#' @param coord,name,from,to,type,window.size,... See 
#' \code{\link{visualizeData}}
#' 
#' @details 
#' \code{ModInosine} specific arguments for \link{visualizeData}:
#' \itemize{
#' \item{\code{colour.bases} - }{a named character vector of \code{length = 4} 
#' for the colours of the individual bases. The names are expected to be 
#' \code{c("G","A","U","C")}}
#' }
#' 
#' @return 
#' \itemize{
#' \item{\code{settings}} {See \code{\link[=Modifier-functions]{settings}}.}
#' \item{\code{aggregate}} {See \code{\link{aggregate}}.}
#' \item{\code{modify}} {See \code{\link{modify}}.}
#' \item{\code{getDataTrack}} {a list of 
#' \code{\link[Gviz:DataTrack-class]{DataTrack}} object.}
#' \item{\code{visualizeData}} {See \code{\link{visualizeDataByCoord}}.}
#' \item{\code{visualizeDataByCoord}} {See \code{\link{visualizeDataByCoord}}.}
#' }
#' 
#' @examples 
#' data(msi,package="RNAmodR")
#' mi <- msi[[1]]
#' settings(mi)
#' aggregate(mi)
#' modify(mi)
#' getDataTrack(mi, "1", mainScore(mi))
NULL

# class ------------------------------------------------------------------------

#' @rdname ModInosine
#' @export
setClass("ModInosine",
         contains = c("Modifier"),
         prototype = list(mod = "I",
                          score = "score",
                          dataType = "PileupSequenceData"))

# constructor ------------------------------------------------------------------

# Create Modifier class from file character, fasta and gff file
#' @rdname ModInosine
#' @export
ModInosine <- function(x, annotation, sequences, seqinfo, ...){
  Modifier("ModInosine", x = x, annotation = annotation, sequences = sequences,
           seqinfo = seqinfo, ...)
}

# settings ---------------------------------------------------------------------

.norm_inosine_args <- function(input){
  minScore <- 0.4
  if(!is.null(input[["minScore"]])){
    minScore <- input[["minScore"]]
    if(!is.numeric(minScore) | minScore < 0 | minScore > 100){
      stop("'minScore' must be numeric with a value between 0 and 100.",
           call. = FALSE)
    }
  }
  args <- .norm_args(input)
  args <- c(args,
            list(minScore = minScore))
  args
}

#' @rdname ModInosine-functions
#' @export
setReplaceMethod(f = "settings", 
                 signature = signature(x = "ModInosine"),
                 definition = function(x, value){
                   x <- callNextMethod()
                   value <- .norm_inosine_args(value)
                   x@arguments[names(value)] <- unname(value)
                   x
                 })

# functions --------------------------------------------------------------------

.calculate_inosine_score <- function(x, letters){
  data <- unlist(x)
  letters <- unlist(letters)
  scores <- 
    data$means.treated.G / (data$means.treated.A + data$means.treated.T + 
                              data$means.treated.C +data$means.treated.G)
  scores[is.infinite(scores) | is.na(scores)] <- 0
  scores[letters != "A"] <- 0
  ans <- S4Vectors::DataFrame(value = scores)
  relist(ans, x)
}

.aggregate_inosine <- function(x){
  mod <- aggregate(sequenceData(x))
  letters <- IRanges::CharacterList(strsplit(as.character(sequences(x)),""))
  score <- .calculate_inosine_score(mod, letters)
  ans <- cbind(S4Vectors::DataFrame(score = unlist(score)$value,
                                    row.names = NULL))
  ans <- relist(ans, mod)
  rownames(ans) <- rownames(mod)
  ans
}

#' @rdname ModInosine-functions
#' @export
setMethod(f = "aggregateData", 
          signature = signature(x = "ModInosine"),
          definition = 
            function(x){
              .aggregate_inosine(x)
            }
)

.get_inosine_score <- function(data){
  list(score = data$score)
}

.find_inosine <- function(x){
  if(!hasAggregateData(x)){
    stop("Something went wrong.")
  }
  letters <- IRanges::CharacterList(strsplit(as.character(sequences(x)),""))
  # get the aggregate data
  mod <- getAggregateData(x)
  coverage <- pileupToCoverage(sequenceData(x))
  # get arguments
  minCoverage <- settings(x,"minCoverage")
  minReplicate <- settings(x,"minReplicate")
  minScore <- settings(x,"minScore")
  # construct logical vector for passing the coverage threshold
  coverage <- rowSums(as.data.frame(unlist(coverage))) >= minCoverage
  coverage <- relist(unname(coverage),mod@partitioning)
  # find inosine positions by looking for A to G conversion at position with 
  # enough coverage
  grl <- ranges(x)
  modifications <- mapply(
    function(m,c,l,r){
      m <- m[l == "A" &
               m$score >= minScore & 
               c,,drop=FALSE] # coverage check
      if(nrow(m) == 0L) return(NULL)
      ans <- constructModRanges(r, m, modType = "I",
                                RNAmodR:::.get_inosine_score, "RNAmodR",
                                "RNAMOD")
      ans
    },
    mod,
    coverage,
    letters,
    grl,
    SIMPLIFY = FALSE)
  f <- !vapply(modifications,
               is.null,
               logical(1))
  modifications <- mapply(
    function(m,name){
      m$Parent <- name
      m
    },
    modifications[f],
    names(grl)[f],
    SIMPLIFY = FALSE)
  modifications <- GenomicRanges::GRangesList(modifications)
  unname(unlist(modifications))
}

#' @rdname ModInosine-functions
#' @export
setMethod("findMod",
          signature = c(x = "ModInosine"),
          function(x){
            .find_inosine(x)
          }
)


# ModSetInosine ----------------------------------------------------------------

#' @rdname ModInosine
#' @export
setClass("ModSetInosine",
         contains = "ModifierSet",
         prototype = list(elementType = "ModInosine"))

#' @rdname ModInosine
#' @export
ModSetInosine <- function(x, annotation = NA, sequences = NA, seqinfo = NA, 
                          ...){
  RNAmodR::ModifierSet("ModInosine", x = x, annotation = annotation,
                       sequences = sequences, seqinfo = seqinfo, ...)
}
