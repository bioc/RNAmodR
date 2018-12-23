#' @include RNAmodR.R
#' @include Modifier-class.R
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
                          dataType = "PileupSequenceData"))
POS_DATA_TYPE_INOSINE <- "PileupSequenceData"


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
  args <- .norm_args(input)
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

.get_inosine_score <- function(data){
  df <- data.frame(x = data$means.G,
                   y = data$means.A,
                   dx = data$sds.G,
                   dy = data$sds.A)
  f <- z ~ log2(x/y) * 10
  scores <- .mutate_with_error(df,f)
  list(score = scores$z,
       sd = scores$dz)
}

.find_inosine <- function(x,
                          args){
  data <- seqdata(x)
  mod <- aggregate(data)
  message("Searching for Inosine...")
  letters <- CharacterList(strsplit(as.character(sequences(data)),""))
  ranges <- split(.get_parent_annotations(ranges(data)),
                  seq_along(ranges(data)))
  coverage <- .aggregate_pile_up_to_coverage(data)
  minCoverage <- args[["minCoverage"]]
  minReplicate <- args[["minReplicate"]]
  coverage <- apply(as.matrix(unlist(coverage)),1,
                    function(row){
                      sum(row > minCoverage) >= minReplicate
                    })
  coverage <- split(unname(coverage),data@partitioning)
  modifications <- mapply(
    function(m,c,l,r){
      rownames(m) <- seq_len(width(r))
      # m <- m[l == "A",]
      m <- m[l == "A" &
               m$means.G > m$means.A &
               !is.na(m$means.G) &
               !is.na(m$means.A) & c,]
      if(nrow(m) == 0L) return(NULL)
      .construct_mod_ranges(r,m,modType = "I",.get_inosine_score,
                            "RNAmodR","RNAMOD")
    },
    mod,
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
            x@modifications <- .find_inosine(x,
                                             .norm_inosine_args(list(...)))
            message("done.")
            x
          }
)
