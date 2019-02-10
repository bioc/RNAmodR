#' @include RNAmodR.R
#' @include Modifier-class.R
#' @include ModifierSet-class.R
NULL

# documentation /---------------------------------------------------------------

#' @name ModInosine
#' @aliases Inosine ModifierInosine ICE-Seq
#' 
#' @author Felix G.M. Ernst [aut]
#' 
#' @title ModInosine
#' 
#' @description 
#' Inosine can be detected in sequence by the conversion of A positions to G.
#' \code{ModInosine} uses to search for Inosine positions.
#' 
#' Only samples named \code{treated} are used for this analysis, since the
#' A to G conversion is common feature among the reverse transcriptases usually
#' emploied. Let us know, if that is not the case, and the class needs to be
#' modified.
#' 
#' \code{\link{ModInosine-functions}}
#' 
#' @details
#' \code{ModInosine} score: the score for reported Inosine positions are between
#' 1 and 100. It is calculated as \code{min(log2(percent(G) /
#' (percent(A)+percent(T)+percent(C))) * 10,100)}.
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
#' @param force See \code{\link{aggregate}}
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

.aggregate_pile_up_to_coverage <- function(data){
  df <- data@unlistData
  replicates <- unique(data@replicate)
  ans  <- IRanges::IntegerList(lapply(seq_along(replicates),
                                      function(i){
                                        rowSums(as.data.frame(df[,data@replicate == i]))
                                      }))
  names(ans) <- paste0("replicate.",replicates)
  ans <- do.call(S4Vectors::DataFrame,ans)
  ans <- IRanges::SplitDataFrameList(ans)
  ans@partitioning <- data@partitioning
  ans
}

.calculate_inosine_score <- function(x, letters){
  data <- x@unlistData
  letters <- letters@unlistData
  scores <- 
    data$means.treated.G / (data$means.treated.A + data$means.treated.T + 
                              data$means.treated.C +data$means.treated.G)
  scores[is.infinite(scores) | is.na(scores)] <- 0
  scores[letters != "A"] <- 0
  ans <- S4Vectors::DataFrame(value = scores)
  ans <- IRanges::SplitDataFrameList(ans)
  ans@partitioning <- x@partitioning
  ans
}

.aggregate_inosine <- function(x){
  message("Aggregating data and calculating scores...")
  data <- seqData(x)
  mod <- aggregate(data)
  letters <- IRanges::CharacterList(strsplit(as.character(sequences(x)),""))
  score <- .calculate_inosine_score(mod, letters)
  ans <- cbind(S4Vectors::DataFrame(score = unlist(score)$value,
                                    row.names = NULL))
  ans <- IRanges::SplitDataFrameList(ans)
  ans@partitioning <- mod@partitioning
  ans
}

#' @rdname ModInosine-functions
#' @export
setMethod(f = "aggregate", 
          signature = signature(x = "ModInosine"),
          definition = 
            function(x, force = FALSE){
              if(missing(force)){
                force <- FALSE
              }
              if(!hasAggregateData(x) || force){
                x@aggregate <- .aggregate_inosine(x)
                x@aggregateValidForCurrentArguments <- TRUE
              }
              x
            }
)

.get_inosine_score <- function(data){
  list(score = data$score)
}

.find_inosine <- function(x){
  letters <- IRanges::CharacterList(strsplit(as.character(sequences(x)),""))
  # get the aggregate data
  mod <- aggregateData(x)
  coverage <- .aggregate_pile_up_to_coverage(seqData(x))
  # get arguments
  minCoverage <- settings(x,"minCoverage")
  minReplicate <- settings(x,"minReplicate")
  minScore <- settings(x,"minScore")
  # construct logical vector for passing the coverage threshold
  coverage <- apply(as.matrix(unlist(coverage)),1,
                    function(row){
                      sum(row > minCoverage) >= minReplicate
                    })
  coverage <- split(unname(coverage),mod@partitioning)
  # find inosine positions by looking for A to G conversion at position with 
  # enough coverage
  grl <- ranges(x)
  modifications <- mapply(
    function(m,c,l,r){
      m <- m[l == "A" &
               m$score >= minScore & 
               c,,drop=FALSE] # coverage check
      if(nrow(m) == 0L) return(NULL)
      ans <- .constructModRanges(r, m, modType = "I",
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
setMethod("modify",
          signature = c(x = "ModInosine"),
          function(x, force = FALSE){
            # get the aggregate data
            x <- aggregate(x, force)
            x@modifications <- .find_inosine(x)
            x@modificationsValidForCurrentArguments <- TRUE
            message("done.")
            x
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
ModSetInosine <- function(x, annotation, sequences, seqinfo, ...){
  ModifierSet("ModInosine", x = x, annotation = annotation,
              sequences = sequences, seqinfo = seqinfo, ...)
}
