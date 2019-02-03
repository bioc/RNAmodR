#' @include RNAmodR.R
#' @include Modifier-class.R
#' @include ModifierSet-class.R
NULL

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
#' \code{\link{ModInosine-functions}}
#' 
#' @details
#' \code{ModInosine} score: the score for reported Inosine positions are between
#' 1 and 100. It is calculated as \code{min(log2(percent(G) /
#' (percent(A)+percent(T)+percent(C))) * 10,100)}.
#' 
NULL

#' @rdname ModInosine
#' @export
setClass("ModInosine",
         contains = c("Modifier"),
         prototype = list(mod = "I",
                          score = "score",
                          dataType = "PileupSequenceData"))

# settings ---------------------------------------------------------------------

#' @name ModInosine-functions
#' 
#' @title Functions for ModInosine
#' 
#' @description
#' All of the functions of \code{\link[RNAmodR:Modifier]{Modifier}} and the
#' \code{\link[RNAmodR:Modifier]{ModifierSet}} classes are inherited by the 
#' \code{ModInosine} and \code{ModSetInosine} classes.
#' 
#' Check below for the specifically implemented functions.
#' 
#' @param x a \code{\link[RNAmodR:Modifier]{Modifier}} or a 
#' \code{\link[RNAmodR:ModifierSet]{ModifierSet}} object. For more details see 
#' also the man pages for the functions mentioned below.
#' @param value See \code{\link{settings}}
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

.norm_inosine_args <- function(input){
  minScore <- 10
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

# constructor ------------------------------------------------------------------

# Create Modifier class from file character, fasta and gff file
#' @rdname ModInosine
#' @export
ModInosine <- function(x, annotation, sequences, seqinfo, ...){
  Modifier("ModInosine", x = x, annotation = annotation, sequences = sequences,
           seqinfo = seqinfo, ...)
}

# functions --------------------------------------------------------------------

.calculate_inosine_score <- function(x){
  data <- x@unlistData
  df <- data.frame(x = data$means.G,
                   y = data$means.A + data$means.T + data$means.C,
                   dx = data$sds.G,
                   dy = data$sds.A + data$sds.T + data$sds.C)
  f <- z ~ log2(x/y) * 10
  scores <- .mutate_with_error(df,f)
  scores$z <- vapply(scores$z,min,numeric(1),100)
  scores$z[is.infinite(scores$z) | is.na(scores$z)] <- 0
  scores$dz[is.infinite(scores$dz) | is.na(scores$dz)] <- 0
  ans <- S4Vectors::DataFrame(value = scores$z,
                              sd = scores$dz)
  ans <- IRanges::SplitDataFrameList(ans)
  ans@partitioning <- x@partitioning
  ans
}

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

.aggregate_inosine <- function(x){
  mod <- aggregate(seqData(x))
  data <- seqData(x)
  coverage <- .aggregate_pile_up_to_coverage(data)
  score <- .calculate_inosine_score(mod)
  ans <- cbind(S4Vectors::DataFrame(score = unlist(score)$value,
                                    sd = unlist(score)$sd,
                                    row.names = NULL),
               unlist(coverage))
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
  list(score = data$score,
       sd = data$sd)
}

.find_inosine <- function(x){
  message("Searching for Inosine ... ", appendLF = FALSE)
  letters <- IRanges::CharacterList(strsplit(as.character(sequences(x)),""))
  # get the aggregate data
  mod <- aggregateData(x)
  # get arguments
  minCoverage <- settings(x,"minCoverage")
  minReplicate <- settings(x,"minReplicate")
  minScore <- settings(x,"minScore")
  # construct logical vector for passing the coverage threshold
  coverage <- mod[,seq_len(unique(IRanges::ncol(mod)) - 2) + 2,drop = FALSE]
  coverage <- apply(as.matrix(unlist(coverage)),1,
                    function(row){
                      sum(row > minCoverage) >= minReplicate
                    })
  coverage <- split(unname(coverage),mod@partitioning)
  # find inosine positions by looking for A to G conversion at position with 
  # enough coverage
  modifications <- mapply(
    function(m,c,l,r){
      rownames(m) <- seq_len(width(r))
      m <- m[l == "A" &
               (m$score - m$sd) >= minScore & 
               c,] # coverage check
      if(nrow(m) == 0L) return(NULL)
      ans <- .constructModRanges(r,
                                 m,
                                 modType = "I",
                                 RNAmodR:::.get_inosine_score,
                                 "RNAmodR",
                                 "RNAMOD")
      ans
    },
    mod[,seq_len(2)],
    coverage,
    letters,
    ranges(x),
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
          function(x, force){
            # get the aggregate data
            x <- aggregate(x, force)
            x@modifications <- .find_inosine(x)
            x@modificationsValidForCurrentArguments <- TRUE
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
