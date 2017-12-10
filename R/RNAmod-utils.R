#' @include RNAmod.R
NULL



#' @rdname setupWorkEnvir
#'
#' @param experimentName a name for the experiment, e.g. "xyz"
#'
#' @return TRUE
#' @export
#'
#' @examples
#' experimentName <- "test"
#' setupWorkEnvir(experimentName)
setMethod(
  f = "setupWorkEnvir", 
  signature = signature(experimentName = "character"),
  definition = function(experimentName = "experiment") {
    assertive::assert_is_scalar(experimentName)
    assertive::assert_is_a_non_missing_nor_empty_string(experimentName)
    
    createFolder = function(folder){
      if (!dir.exists(folder)) {
        dir.create(folder, recursive = TRUE)
      }
    }
    
    dataFolder <- paste0(experimentName, "/data/")
    resultFolder <- paste0(experimentName, "/results/")
    
    createFolder(dataFolder)
    createFolder(resultFolder)
    
    file.copy(from = system.file("extdata", "example_experiment_layout.csv",
                                 package = "RNAmod"),
              to = paste0(dataFolder,"experiment_layout.csv") )
    
    message(paste0("work enviroment '",experimentName,"' created."))
    
    return(invisible(TRUE))
  }
)

# reading BAM input ------------------------------------------------------------

#' converts results from scanBam to DataFrame
#' returns a simplified DataFrame for read data read out from a BAM
#' file
#' 
#' @importFrom stats setNames
#' @importClassesFrom S4Vectors DataFrame
.convert_bam_to_DataFrame <- function(dataList, 
                                      param){
  .unlist <- function (x){
    ## do.call(c, ...) coerces factor to integer, which is undesired
    x1 <- x[[1L]]
    if (is.factor(x1)) {
      structure(unlist(x), class = "factor", levels = levels(x1))
    } else {
      do.call(c, x)
    }
  } 
  dataList <- unname(dataList) # names not useful in unlisted result
  elts <- stats::setNames(bamWhat(param), bamWhat(param))
  
  lst <- lapply(elts, function(elt) .unlist(lapply(dataList, "[[", elt)))
  
  df <- do.call(S4Vectors::DataFrame, lst)
  return(df)
}

#' scans BAM files and converts them to DataFrame each 
#' 
#' @import Biostrings
.get_bam_data <- function(files, 
                          param){
  
  res <- c()
  files <- files[grepl(".bam",files)]
  if( length(files) > 0){
    
    FUN <-  function(x, param){
      requireNamespace("Rsamtools", quietly = TRUE)
      
      y <- scanBam(x, param=param)
      y <- .convert_bam_to_DataFrame(y, param)
      return(y)
    }
    
    res <- bplapply(files,
                    FUN,
                    param = param)
    names(res) <- files
  }
  return(res)
  
}

# Returns the number of reads contained in a bam file for parameters set in
# param
.get_bam_read_count <- function(files, quality){
  
  res <- c()
  files <- files[grepl(".bam",files)]
  if( length(files) > 0){
    for(i in 1:length(files)) {
      res[[i]] <- countBam(files[[i]], 
                           param = ScanBamParam(mapqFilter = quality))$records
    }
  }
  return(res)
}


#' return parameters to be used by scanBam
#' gRangeInput a GRanges object containing the ranges for search in the BAM file
#' quality quality argument used for scanBamParam
#' 
#' @importFrom GenomeInfoDb seqnames
#' @importFrom Rsamtools ScanBamParam
.assemble_scanBamParam <- function(gRangeInput,
                                   quality){
  gRangeList <- GRangesList()
  for(i in 1:length(unique(GenomeInfoDb::seqnames(gRangeInput)))){
    gRangeList <- append(gRangeList, GRangesList(gRangeInput[GenomeInfoDb::seqnames(gRangeInput) == unique(GenomeInfoDb::seqnames(gRangeInput))[i],]) )
  }
  names(gRangeList) <- unique(GenomeInfoDb::seqnames(gRangeInput))
  which <- gRangeList
  what <- c("rname", "strand", "pos", "qwidth", "seq", "mapq")
  param <- Rsamtools::ScanBamParam(which=which, 
                                   what=what, 
                                   mapqFilter = quality)
  return(param)
}
