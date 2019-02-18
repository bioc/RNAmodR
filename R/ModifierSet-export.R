#' @include RNAmodR.R
#' @include Modifier-class.R
#' @include ModifierSet-class.R
NULL

#' @name ModifierSet-export
#' 
#' @title Exporting data to other classes and files
#' 
#' @description 
#' title
#' 
#' @param assays a \code{\link[=ModifierSet-class]{ModifierSet}} object
#' @param ... Optional arguments:
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
#' }
NULL

.get_rse_args <- function(input){
  modificationsOnly <- TRUE
  flanking <- 10L
  allPositions <- FALSE
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
  args <- list(modificationsOnly = modificationsOnly,
               flanking = flanking,
               allPositions = allPositions)
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
      data <- compareByCoord(assays, mod)
    } else {
      data <- compareByCoord(assays, mod, flanking = args[["flanking"]])
    }
    rowData <- data[,(colnames(data) %in% c("names","positions","mod"))]
    colnames(rowData) <- c("Parent","positions","mod")
    data <- data[,!(colnames(data) %in% c("names","positions","mod"))]
    seqnames <- seqnames(mod)[match(rowData$names,as.character(mod$Parent))]
    strand <- strand(mod)[match(rowData$names,as.character(mod$Parent))]
  } else {
    
  }
  ans <- SummarizedExperiment(list(aggregate = data), 
                              rowRanges = GRanges(seqnames,
                                                  ranges = IRanges::IRanges(start = as.integer(rowData$positions), width = 1L),
                                                  strand = strand,
                                                  rowData[,colnames(rowData) != "positions"]))
  return(ans)
}


#' @rdname ModifierSet-export
#' @export
setMethod("SummarizedExperiment", "ModifierSet",
          function(assays, ...) 
            .from_ModifierSet_to_RangedSummarizedExperiment(assays, ...))
