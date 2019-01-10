#' @include RNAmodR.R
#' @include Modifier-class.R
#' @include ModifierSet-class.R
NULL

#' @name ModInosine
#' @aliases Inosine ModifierInosine
#' 
#' @title ModInosine
#' @description 
#' title
#' 
#' \code{score}: the score for reported Inosine positions are between 1 and 100.
#' It is calculated as \code{min(log2(percent(G)/percent(A)) * 10,100)}.
#' 

NULL

#' @rdname ModInosine
#' @export
setClass("ModInosine",
         contains = c("Modifier"),
         prototype = list(mod = "I",
                          score = "score",
                          dataClass = "PileupSequenceData"))

setMethod(
  f = "initialize", 
  signature = signature(.Object = "ModInosine"),
  definition = function(.Object,
                        bamfiles,
                        fasta,
                        gff) {
    .Object <- callNextMethod(.Object,
                              bamfiles,
                              fasta,
                              gff)
    return(.Object)
  }
)

# constructors -----------------------------------------------------------------

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

setGeneric( 
  name = "ModInosine",
  def = function(x,
                 ...) standardGeneric("ModInosine")
)
# Create Modifier class from file character, fasta and gff file
#' @rdname ModInosine
#' @export
setMethod("ModInosine",
          signature = c(x = "character"),
          function(x,
                   fasta,
                   gff,
                   modifications = NULL,
                   ...){
            args <- .norm_inosine_args(list(...))
            ans <- .ModFromCharacter("ModInosine",
                                     x,
                                     fasta,
                                     gff,
                                     args)
            ans <- .norm_modifications(ans,
                                       args)
            ans
          }
)

# Create Modifier class from bamfiles, fasta and gff file
#' @rdname ModInosine
#' @export
setMethod("ModInosine",
          signature = c(x = "BamFileList"),
          function(x,
                   fasta,
                   gff,
                   modifications = NULL,
                   ...){
            args <- .norm_inosine_args(list(...))
            ans <- .ModFromCharacter("ModInosine",
                                     x,
                                     fasta,
                                     gff,
                                     args)
            ans <- .norm_modifications(ans,
                                       args)
            ans
          }
)

# Create Modifier class from existing SequenceData
#' @rdname ModInosine
#' @export
setMethod("ModInosine",
          signature = c(x = "SequenceData"),
          function(x,
                   modifications = NULL,
                   ...){
            args <- .norm_inosine_args(list(...))
            ans <- .ModFromSequenceData("ModInosine",
                                        x,
                                        args)
            ans <- .norm_modifications(ans,
                                       args)
            ans
          }
)

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
  ans <- DataFrame(value = scores$z,
                   sd = scores$dz)
  ans <- SplitDataFrameList(ans)
  ans@partitioning <- x@partitioning
  ans
}

.aggregate_pile_up_to_coverage <- function(data){
  df <- data@unlistData
  replicates <- unique(data@replicate)
  ans  <- IntegerList(lapply(seq_along(replicates),
                             function(i){
                               rowSums(as.data.frame(df[,data@replicate == i]))
                             }))
  names(ans) <- paste0("replicate.",replicates)
  ans <- do.call("DataFrame",ans)
  ans <- SplitDataFrameList(ans)
  ans@partitioning <- data@partitioning
  ans
}

.aggregate_inosine <- function(x,...){
  mod <- aggregate(seqData(x))
  data <- seqData(x)
  coverage <- .aggregate_pile_up_to_coverage(data)
  score <- .calculate_inosine_score(mod)
  ans <- cbind(DataFrame(score = unlist(score)$value,
                         sd = unlist(score)$sd,
                         row.names = NULL),
               unlist(coverage))
  ans <- SplitDataFrameList(ans)
  ans@partitioning <- mod@partitioning
  ans
}

#' @name Modifier
#' @export
setMethod(f = "aggregate", 
          signature = signature(x = "ModInosine"),
          definition = 
            function(x,
                     ...){
              if(!hasAggregateData(x)){
                x@aggregate <- .aggregate_inosine(x,...)
              }
              x
            }
)

.get_inosine_score <- function(data){
  list(score = data$score,
       sd = data$sd)
}

.find_inosine <- function(x,
                          args){
  message("Searching for Inosine...")
  letters <- CharacterList(strsplit(as.character(sequences(x)),""))
  ranges <- ranges(x)
  ranges <- split(.get_parent_annotations(ranges),
                  seq_along(ranges))
  # get the aggregate data
  mod <- aggregateData(x)
  # get arguments
  minCoverage <- args[["minCoverage"]]
  minReplicate <- args[["minReplicate"]]
  minScore <- args[["minScore"]]
  # construct logical vector for passing the coverage threshold
  coverage <- mod[,seq_len(unique(ncol(mod)) - 2) + 2,drop = FALSE]
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
      ans <- .construct_mod_ranges(r,
                                   m,
                                   modType = "I",
                                   .get_inosine_score,
                                   "RNAmodR",
                                   "RNAMOD")
      ans
    },
    mod[,seq_len(2)],
    coverage,
    letters,
    ranges)
  modifications <- GRangesList(modifications[!vapply(modifications,
                                                     is.null,
                                                     logical(1))])
  unname(unlist(modifications))
}

#' @rdname ModInosine
#' @export
setMethod("modify",
          signature = c(x = "ModInosine"),
          function(x,
                   ...){
            # get the aggregate data
            x <- aggregate(x)
            x@modifications <- .find_inosine(x,
                                             .norm_inosine_args(list(...)))
            message("done.")
            x
          }
)

# ModSetInosine ------------------------------------------------------------

#' @rdname ModInosine
#' @export
setClass("ModSetInosine",
         contains = "ModifierSet",
         prototype = list(elementType = "ModInosine"))

#' @rdname ModRiboMethSeq
#' @export
ModSetInosine <- function(x, fasta = NA, gff = NA){
  ModifierSet("ModInosine", x, fasta = fasta, gff = gff)
}