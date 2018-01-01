#' @include class-RNAmod-type.R
NULL

RNAMOD_M7G_P_THRESHOLD <- 0.05
RNAMOD_M7G_SIGMA_THRESHOLD <- 5


#' @rdname mod
#'
#' @description 
#' \code{mod_m7G}
#'
#' @return
#' @export
#'
#' @examples
setClass("mod_m7G",
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
  signature = signature(object = "mod_m7G",
                        counts = "numeric",
                        gff = "GRanges",
                        data = "DataFrame"),
  definition = function(object,
                        counts,
                        gff,
                        data) {
    return(NA)
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
  signature = signature(object = "mod_m7G",
                        gff = "GRanges",
                        seq = "FaFile",
                        data = "list"),
  definition = function(object,
                        gff,
                        seq,
                        data) {
    # browser()
    # Process only genes found in all datasets
    IDs <- lapply(data,names)
    IDs <- Reduce(intersect, IDs)
    
    # detect modification per transcript
    res <- lapply(IDs,
    # res <- bplapply(IDs,
                  FUN = .analyze_transcript,
                  data,
                  gff,
                  seq)
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
.get_m7G_data <- function(ID,data){
  res <- lapply(data,function(x){
    return(x[[ID]][["default"]])
  })
  return(res)
}


# merge positions in one transcript
.analyze_transcript <- function(ID,data,gff,fafile){
  data <- .get_m7G_data(ID,data)
  browser()
  
  # get sequence of transcript and subset gff for single transcript data
  gff <- gff[(is.na(S4Vectors::mcols(gff)$ID) & 
                S4Vectors::mcols(gff)$Name == ID) |
               (!is.na(S4Vectors::mcols(gff)$ID) & 
                  S4Vectors::mcols(gff)$ID == ID),]
  strand <- unique(as.character(strand(gff)))
  seq <- getSeq(fafile,gff)[[1]]
  
  # detect all G positions
  loc <- stringr::str_locate_all(as.character(seq), "G")
  loc <- loc[[1]][,"start"]
  
  # Convert local G position to global positions
  locations <- .convert_local_to_global_locations(gff, loc)
  
  # if a transcript is encountered having ID present in the global variable
  # RNAmod_sample_transcript trigger sample plotting to access quality on
  # know locations
  # ID %in% options("RNAmod_sample_transcripts")
  
  name <- NULL
  if(ID %in% options("RNAmod_sample_transcripts")){
    name <- paste0(ID,
                   "_m7G_")
  }
  # Retrieve m7G positions
  modifications <- lapply(locations,
                          .check_for_m7G,
                          data,
                          strand,
                          name)
  
  # name the locations based on sequence position
  names(modifications) <- paste0("G_",loc)
  modifications <-  modifications[!vapply(modifications,is.null,logical(1))]
  return(modifications)
}

.check_for_m7G <- function(location, data, strand, name = NULL){
  # if( location == 456158){
  #   browser()
  # }
  
  # m7G expects a high number of reads at the N+1 position
  if(as.character(strand) == "-"){
    testLocation <- location-1
  } else {
    testLocation <- location+1
  }
  # number of replicates
  nbSamples <- length(data)
  # merge data for positions
  
  testData <- lapply(data,function(dataPerReplicate){
    return(dataPerReplicate[names(dataPerReplicate) == testLocation])
  })
  testData <- unlist(testData)
  baseData <- lapply(data,function(dataPerReplicate){
    # This might also be an option
    # return(dataPerReplicate[names(dataPerReplicate) < (location+50) & 
    #               names(dataPerReplicate) > (location-50) &
    #               names(dataPerReplicate) != testLocation])
    
    return(dataPerReplicate[names(dataPerReplicate) != testLocation])
  })
  baseData <- unlist(baseData)
  # Fill missing positions with 0?
  
  # if no read is at the expected location
  if(length(testData) == 0){
    return(NULL)
  }
  # if no reads can be used comparison
  if(length(baseData) < (3 * nbSamples)){
    return(NULL)
  }
  
  # overall mean and sd
  sd <-  stats::sd(unlist(baseData))
  mean <-  mean(unlist(baseData))
  # Use the sigma level as value for signal strength
  sigStrength <- (as.numeric(as.character(testData)) - mean) %/% sd
  sigStrength.mean <- mean(sigStrength)
  sigStrength.sd <- stats::sd(sigStrength)
  
  # Since normality of distribution can not be assumed use the MWU
  # generate p.value for single position
  p.value <- wilcox.test(baseData, testData)$p.value
  
  if( sigStrength.mean > RNAMOD_M7G_SIGMA_THRESHOLD &&
        p.value < RNAMOD_M7G_P_THRESHOLD  ){
  # if( sigStrength.mean > RNAMOD_M7G_SIGMA_THRESHOLD  ){
    
    # if location is among sample location (name is not null)
    # plot the data
    if(!is.null(name)){
      .plot_sample_data(baseData, 
                        testData, 
                        paste0(name,location))
    }
    
    # Return data
    return(list(location = location,
                signal = sigStrength.mean,
                signal.sd = sigStrength.sd,
                p.value = p.value,
                nbsamples = nbSamples))
  }
  return(NULL)
}

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
  signature = signature(object = "mod_m7G",
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
