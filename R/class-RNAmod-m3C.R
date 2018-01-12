#' @include class-RNAmod-type.R
NULL

RNAMOD_M3C_ROLLING_MEAN_WINDOW_WIDTH <- 9
RNAMOD_M3C_ARREST_RATE <- 0.95
RNAMOD_M3C_P_THRESHOLD <- 0.05
RNAMOD_M3C_SIGMA_THRESHOLD <- 3


#' @rdname mod
#'
#' @description 
#' \code{mod_m3C}
#'
#' @return
#' @export
#'
#' @examples
setClass("mod_m3C",
         contains = "mod",
         prototype = list(modType = "m3C")
)


#' @rdname convertReadsToPositions
#'
#' @description
#' \code{mod_m3C}: calls the default method
#'
#' @return
#' @export
#'
#' @examples
setMethod(
  f = "convertReadsToPositions",
  signature = signature(object = "mod_m3C",
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
#' \code{mod_m3C}
#' 
#' @return
#' @export
#' 
#' @importFrom stringr str_locate_all
#'
#' @examples
setMethod(
  f = "parseMod",
  signature = signature(object = "mod_m3C",
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
                  FUN = .analyze_m3C_transcript,
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


# returns the position data for m3C analysis.
# each entry in list a result for a gene of all replicates
.get_m3C_data <- function(ID,data){
  res <- lapply(data,function(x){
    return(x[[ID]][["default"]])
  })
  return(res)
}


# merge positions in one transcript
.analyze_m3C_transcript <- function(ID,data,gff,fafile){
  data <- .get_m3C_data(ID,data)
  # do not take into account position 1
  data <- lapply(data, function(x){
    x[as.numeric(names(x)) == 1] <- 0
    x
  })
  
  # get sequence of transcript and subset gff for single transcript data
  gff <- .subset_gff_for_unique_transcript(gff, ID)
  seq <- .get_seq_for_unique_transcript(gff,fafile,ID)
  
  # detect all C positions
  loc <- stringr::str_locate_all(as.character(seq), "C")
  loc <- loc[[1]][,"start"]
  if(length(loc) == 0) return(NULL)
  
  # Convert local G position to global positions
  locations <- .convert_local_to_global_locations(gff, loc)
  
  # if a transcript is encountered having ID present in the global variable
  # RNAmod_sample_transcript trigger sample plotting to access quality on
  # know locations
  # ID %in% options("RNAmod_sample_transcripts")
  
  name <- NULL
  if(ID %in% options("RNAmod_sample_transcripts")){
    name <- paste0(ID,
                   "_m3C_")
  }
  
  # Calculate the arrest rate per position
  arrestData <- lapply(data, .get_arrest_rate)
  
  # Retrieve m3C positions
  modifications <- lapply(locations,
                          .check_for_m3C,
                          data,
                          arrestData,
                          name)
  
  if(length(modifications) == 0) return(NULL)
  # name the locations based on sequence position
  names(modifications) <- paste0("C_",loc)
  modifications <-  modifications[!vapply(modifications,is.null,logical(1))]
  if(length(modifications) == 0) return(NULL)
  return(modifications)
}

.check_for_m3C <- function(location, 
                           data, 
                           arrestData,
                           name = NULL){
  # if( location == ){
  #   browser()
  # }
  
  locTest <- .calc_m3C_test_values(location,
                                   data,
                                   arrestData)
  # If insufficient data is present
  if(is.null(locTest)) return(NULL)
  
  # dynamic threshold based on the noise of the signal (high sd)
  if(!.validate_m3C_pos(RNAMOD_M3C_SIGMA_THRESHOLD, 
                        RNAMOD_M3C_P_THRESHOLD, 
                        locTest$sig.mean, 
                        locTest$p.value) ) return(NULL)
  
  # browser()
  # if location is among sample location (name is not null)
  # plot the data
  if(!is.null(name)){
    testData <- .aggregate_location_data(data, (location+1))
    baseData <- .aggregate_not_location_data(data, (location+1))
    .plot_sample_data(.create_m3C_plot_data(testData,
                                            baseData,
                                            paste0(name,location)), 
                      paste0(name,location))
  }
  
  # Return data
  return(list(location = location,
              signal = locTest$sig.mean,
              signal.sd = locTest$sig.sd,
              p.value = locTest$p.value,
              nbsamples = locTest$n))
}

.calc_m3C_test_values <- function(location,
                                  data,
                                  arrestData,
                                  locs){
  # do not take into account position 1
  if(location == 1) return(NULL)
  
  # merge data for positions
  # data on the N+1 location
  testData <- .aggregate_location_data(data, (location+1))
  testData <- testData[testData > 0]
  # data on the arrect direction
  testArrestData <- .aggregate_location_data(arrestData, location)
  # base data to compare against
  # use only G position
  baseData <- .aggregate_not_location_data(data, (location+1))
  
  # number of replicates
  n <- length(data)
  # if not enough data is present
  if(length(testData) == 0 | 
     length(baseData) < (3*n)) return(NULL)
  # No read arrest detectable
  if( sum(testArrestData) < 0 ) return(NULL)
  # To low arrest detectable
  testArrest <- length(testArrestData[testArrestData >= RNAMOD_M3C_ARREST_RATE])
  if( length(testArrestData) != testArrest ) {
    return(NULL)
  }
  
  # get test values
  # overall mean and sd
  mean <-  mean(baseData)
  sd <-  stats::sd(baseData)
  thresholdAddition <- abs(sd %/% mean)
  
  # Use the sigma level as value for signal strength
  sig <- (as.numeric(as.character(testData)) - mean) %/% sd
  sig.mean <- mean(sig)
  sig.sd <- stats::sd(sig)
  # Since normality of distribution can not be assumed use the MWU
  # generate p.value for single position
  p.value <- suppressWarnings(wilcox.test(baseData, testData)$p.value)
  
  return(list(thresholdAddition = thresholdAddition,
              sig = sig,
              sig.mean = sig.mean,
              sig.sd = sig.sd,
              p.value = p.value,
              n = n))
}

.validate_m3C_pos <- function(sig.threshold, 
                              p.threshold, 
                              sig, 
                              p.value){
  ((sig > sig.threshold &&
      p.value <= p.threshold) ||
     (sig > sig.threshold &&
        !.get_use_p()))
}

.create_m3C_plot_data <- function(testData, baseData, name){
  data.frame(x = c(rep(name,(length(testData) + 
                               length(baseData)))),
             y = c(testData, 
                   baseData),
             group = c(rep("Position",length(testData)), 
                       rep("Baseline", length(baseData))))
}

#' @rdname mergePositionsOfReplicates
#'
#' @description
#' \code{mod_m3C}
#'
#' @return
#' @export
#'
#' @examples
setMethod(
  f = "mergePositionsOfReplicates",
  signature = signature(object = "mod_m3C",
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
