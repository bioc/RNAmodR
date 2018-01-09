#' @include class-RNAmod-type.R
NULL

RNAMOD_D_P_THRESHOLD <- 0.05
RNAMOD_D_SIGMA_THRESHOLD <- 5


#' @rdname mod
#'
#' @description 
#' \code{mod_D}
#'
#' @return
#' @export
#'
#' @examples
setClass("mod_D",
         contains = "mod",
         prototype = list(modType = "D")
)


#' @rdname convertReadsToPositions
#'
#' @description
#' \code{mod_D}: calls the default method
#'
#' @return
#' @export
#'
#' @examples
setMethod(
  f = "convertReadsToPositions",
  signature = signature(object = "mod_D",
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
#' \code{mod_D}
#' 
#' @return
#' @export
#' 
#' @importFrom stringr str_locate_all
#'
#' @examples
setMethod(
  f = "parseMod",
  signature = signature(object = "mod_D",
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


# returns the position data for D analysis.
# each entry in list a result for a gene of all replicates
.get_D_data <- function(ID,data){
  res <- lapply(data,function(x){
    return(x[[ID]][["default"]])
  })
  return(res)
}


# merge positions in one transcript
.analyze_transcript <- function(ID,data,gff,fafile){
  data <- .get_D_data(ID,data)
  # browser()
  
  # get sequence of transcript and subset gff for single transcript data
  gff <- gff[(is.na(S4Vectors::mcols(gff)$ID) & 
                S4Vectors::mcols(gff)$Name == ID) |
               (!is.na(S4Vectors::mcols(gff)$ID) & 
                  S4Vectors::mcols(gff)$ID == ID),]
  strand <- unique(as.character(strand(gff)))
  seq <- getSeq(fafile,gff)[[1]]
  
  # detect all G positions
  loc <- stringr::str_locate_all(as.character(seq), "U")
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
                   "_D_")
  }
  # Retrieve D positions
  modifications <- lapply(locations,
                          .check_for_D,
                          data,
                          strand,
                          name)
  
  # name the locations based on sequence position
  names(modifications) <- paste0("U_",loc)
  modifications <-  modifications[!vapply(modifications,is.null,logical(1))]
  return(modifications)
}

.check_for_D <- function(location, data, strand, name = NULL){
  # if( location == 456158){
  #   browser()
  # }
  
  # D expects a high number of reads at the N+1 position
  if(as.character(strand) == "-"){
    testLocation <- location-1
  } else {
    testLocation <- location+1
  }
  # number of replicates
  nbSamples <- length(data)
  # merge data for positions
  
  testData <- lapply(data,function(dataPerReplicate){
    dataPerReplicate <- dataPerReplicate[dataPerReplicate != 0]
    return(dataPerReplicate[names(dataPerReplicate) == testLocation])
  })
  testData <- unlist(testData)
  baseData <- lapply(data,function(dataPerReplicate){
    dataPerReplicate <- dataPerReplicate[dataPerReplicate != 0]
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
  p.value <- suppressWarnings(wilcox.test(baseData, testData)$p.value)
  
  useP <- getOption("RNAmod_use_p")
  if(!assertive::is_a_bool(useP)){
    useP <- as.logical(useP[[1]])
    warning("The option 'RNAmod_use_p' is not a single logical value. ",
            "Please set 'RNAmod_use_p' to TRUE or FALSE.",
            .call = FALSE)
  }
  if( (sigStrength.mean > RNAMOD_D_SIGMA_THRESHOLD &&
       p.value < RNAMOD_D_P_THRESHOLD) || 
      (sigStrength.mean > RNAMOD_D_SIGMA_THRESHOLD &&
       !useP ) ){
    # if( sigStrength.mean > RNAMOD_D_SIGMA_THRESHOLD  ){
    
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
#' \code{mod_D}
#'
#' @return
#' @export
#'
#' @examples
setMethod(
  f = "mergePositionsOfReplicates",
  signature = signature(object = "mod_D",
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
