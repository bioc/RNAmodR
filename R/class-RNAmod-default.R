#' @include class-RNAmod-type.R
NULL

#' @rdname mod
#'
#' @description 
#' \code{mod_m7G}
#'
#' @return
#' @export
#'
#' @examples
setClass("mod_default",
         contains = "mod")


#' @rdname convertReadsToPositions
#'
#' @description
#' \code{mod_m7G}: calls the default method
#'
#' @return
#' @export
#'
#' @examples
setMethod(
  f = "convertReadsToPositions",
  signature = signature(object = "mod_default",
                        counts = "numeric",
                        gff = "GRanges",
                        data = "DataFrame"),
  definition = function(object,
                        counts,
                        gff,
                        data) {
    browser()
    strand <- as.character(strand(gff))
    if(strand == "-"){
      pos <- data$pos + data$qwidth - 1
    } else {
      pos <- data$pos
    }
    # Normalize counts per positions against million of reads in BamFile
    posData <- table(pos)/(counts/10^6)
    
    return(posData)
  }
)


#' @rdname parseMod
#' 
#' @description 
#' \code{mod_m7G}
#' 
#' @return
#' @export
#' 
#' @importFrom stringr str_locate_all
#'
#' @examples
setMethod(
  f = "parseMod",
  signature = signature(object = "mod_default",
                        gff = "GRanges",
                        seq = "FaFile",
                        data = "list"),
  definition = function(object,
                        gff,
                        seq,
                        data) {
    return(NA)
  }
)

#' @rdname mergePositionsOfReplicates
#'
#' @description
#' \code{mod_m7G}
#'
#' @return
#' @export
#'
#' @examples
setMethod(
  f = "mergePositionsOfReplicates",
  signature = signature(object = "mod_default",
                        gff = "GRanges",
                        seq = "FaFile",
                        data = "list"),
  definition = function(object,
                        gff,
                        seq,
                        data) {
    # Process only genes found in all datasets
    IDs <- lapply(data,names)
    IDs <- Reduce(intersect, IDs)
    res <- lapply(IDs,
                  FUN = .merge_positions,
                  data)
    names(res) <- IDs
    res <- res[!is.null(res)]
    
    # If not results are present return NA instead of NULL
    if(is.null(res)){
      return(NA)
    }
    return(res)
  }
)

# returns the position data for m7G analysis.
# each entry in list a result for a gene of all replicates
.get_default_data <- function(ID,data){
  res <- lapply(data,function(x){
    return(x[[ID]][["default"]])
  })
  return(res)
}


# merge positions in one transcript
.merge_positions <- function(ID,data){
  data <- .get_default_data(ID,data)
  positions <- unique(unlist(lapply(data,names)))
  res <- lapply(positions,
                FUN = .merge_position,
                data)
  df <- data.frame(pos = unlist(lapply(res,"[[","pos")),
                   mean = unlist(lapply(res,"[[","mean")),
                   sd = unlist(lapply(res,"[[","sd")))
  return(df)
}

# merge position data for one position
.merge_position <- function(pos,data){
  data <- unlist(lapply(data, "[[",pos))
  data[is.na(data)] <- 0
  return(list(pos = pos,
              mean = mean(data),
              sd = stats::sd(data)))
}

