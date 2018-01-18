#' @include class-RNAmod-mod-type.R
NULL

RNAMOD_M7G_ROLLING_MEAN_WINDOW_WIDTH <- 9
RNAMOD_M7G_ARREST_RATE <- 0.95
RNAMOD_M7G_P_THRESHOLD <- 0.05
RNAMOD_M7G_SIGMA_THRESHOLD <- 3


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
         contains = "mod",
         prototype = list(modType = "m7G")
)


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
#' @importFrom BiocParallel bplapply
#'
#' @examples
setMethod(
  f = "parseMod",
  signature = signature(object = "mod_m7G",
                        gff = "GRanges",
                        fafile = "FaFile",
                        data = "list"),
  definition = function(object,
                        gff,
                        fafile,
                        data) {
    message("Parsing data for m7G modifications...")
    # browser()
    # Process only genes found in all datasets
    IDs <- lapply(data,names)
    IDs <- Reduce(intersect, IDs)
    
    # detect modification per transcript
    # res <- lapply(IDs,
    res <- BiocParallel::bplapply(IDs,
                                  FUN = .analyze_M7G_transcript_prep,
                                  data = data,
                                  gff = gff,
                                  fafile = fafile)
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
.get_M7G_data <- function(ID,data){
  res <- lapply(data,function(x){
    return(x[[ID]][["default"]])
  })
  return(res)
}

# detect and merge modification positions in one transcript
.analyze_M7G_transcript_prep <- function(ID,
                                         data,
                                         gff,
                                         fafile){
  # debug
  if( getOption("RNAmod_debug") ){
    .print_transcript_info(paste(ID," prep"), "")
  }
  
  data <- .get_M7G_data(ID,data)
  # do not take into account position 1
  data <- lapply(data, function(x){
    x[as.numeric(names(x)) == 1] <- 0
    x
  })
  # get sequence of transcript and subset gff for single transcript data
  gff <- .subset_gff_for_unique_transcript(gff, ID)
  seq <- .get_seq_for_unique_transcript(gff,fafile,ID)
  # detect all G positions
  locations <- stringr::str_locate_all(as.character(seq), "G")
  locations <- locations[[1]][,"start"]
  if(length(locations) == 0) return(NULL)
  # Convert local G position to global positions
  globalLocations <- .convert_local_to_global_locations(gff, locations)
  
  res <- .analyze_M7G_transcript(ID = ID,
                                 data = data,
                                 globalLocations = globalLocations,
                                 iterationN = 1)
  res <- res[order(as.numeric(unlist(lapply(res, "[[", "location"))))]
  return(res)
}

# merge positions in one transcript
.analyze_M7G_transcript <- function(ID,
                                    data,
                                    globalLocations,
                                    iterationN){
  if( iterationN > .get_transcript_max_iteration()) return(NULL)
  # debug
  if( getOption("RNAmod_debug") ){
    .print_transcript_info(paste(ID," detect"), iterationN)
  }
  
  # Retrieve m7G positions
  modifications <- lapply(globalLocations,
                          .check_for_M7G,
                          data,
                          globalLocations)
  if(length(modifications) == 0) return(NULL)
  # name the locations based on sequence position
  names(modifications) <- paste0(ID,"_G_",globalLocations)
  modifications <-  modifications[!vapply(modifications,is.null,logical(1))]
  if(length(modifications) == 0) return(NULL)
  # check if by masking found modifications positions additional positions are
  # picked up - remember the +1 location data
  modLocations <- as.numeric(unlist(lapply(modifications, "[[","location")))
  # do not take into account positions found as modification
  data <- .mask_data(data, modLocations)
  globalLocations <- globalLocations[!(globalLocations %in% modLocations)]
  return(append(modifications,.analyze_M7G_transcript(ID = ID,
                                        data = data,
                                        globalLocations = globalLocations,
                                        (iterationN+1))))
}

# mask data for know modification locations
# extend to N+2,N+3 location?
.mask_data <- function(data, modLocations){
  lapply(data, function(x){
    x[as.numeric(names(x)) %in% (modLocations+1)] <- 
      x[as.numeric(names(x)) %in% (modLocations+1)]*(1-RNAMOD_M7G_ARREST_RATE)
    x
  })
}

# check for m7G at given position
.check_for_M7G <- function(location, 
                           data,
                           locs){
  # short cut if amount of data is not sufficient
  if( is.null(.do_M7G_pretest(location,
                            data))) return(NULL)
  # If potential modification right in front of current location
  if(length(locs[locs == (location-1)]) != 0) {
    locTestPre <- .check_for_M7G((location-1), 
                                 .mask_data(data, location),
                                 locs[locs != location])
    if(!is.null(locTestPre)){
      # udpate data accordingly
      data <- .mask_data(data, (location-1))
    }
  }
  # If potential modification right after of current location
  if(length(locs[locs == (location+1)]) != 0) {
    locTestPost <- .check_for_M7G((location+1), 
                                  .mask_data(data, location),
                                  locs[locs != location])
    if(!is.null(locTestPost)){
      # udpate data accordingly
      data <- .mask_data(data, (location+1))
    }
  }
  # Calculate the arrest rate per position
  arrestData <- lapply(data, .get_arrest_rate)
  # get test result for the current location
  locTest <- .calc_M7G_test_values(location,
                                   data,
                                   arrestData)
  # If insufficient data is present
  if(is.null(locTest)) return(NULL)
  # dynamic threshold based on the noise of the signal (high sd)
  if(!.validate_M7G_pos(RNAMOD_M7G_SIGMA_THRESHOLD, 
                        RNAMOD_M7G_P_THRESHOLD, 
                        locTest$sig.mean, 
                        locTest$p.value) ) {
    # debug
    if( getOption("RNAmod_debug") ){
      .print_location_info(paste(location,"_no"),locs)
    }
    return(NULL)
  }
  # debug
  if( getOption("RNAmod_debug") ){
    .print_location_info(paste(location,"_yes"), locs)
  }
  # Return data
  return(list(location = location,
              signal = locTest$sig.mean,
              signal.sd = locTest$sig.sd,
              p.value = locTest$p.value,
              nbsamples = locTest$n))
}

# check if any data is available to proceed with test
.do_M7G_pretest <- function(location,
                            data){
  # do not take into account position 1
  if(location == 1) return(NULL)
  
  # merge data for positions
  # data on the N+1 location
  testData <- .aggregate_location_data(data, (location+1))
  testData <- testData[testData > 0]
  # base data to compare against
  baseData <- .aggregate_not_location_data(data, (location+1))
  
  # number of replicates
  n <- length(data)
  # if not enough data is present
  if(length(testData) == 0 | 
     length(baseData) < (3*n)) return(NULL)
  return(list(n = n,
              testData = testData,
              baseData = baseData))
}

# test for m7G at current location
.calc_M7G_test_values <- function(location,
                                  data,
                                  arrestData){
  # short cut if amount of data is not sufficient
  pretestData <- .do_M7G_pretest(location,
                                 data)
  if(is.null(pretestData)) return(NULL)
  
  # data from pretest
  testData <- pretestData$testData
  baseData <- pretestData$baseData
  n <- pretestData$n
  
  # data on the arrect direction
  testArrestData <- .aggregate_location_data(arrestData, location)
  # No read arrest detectable
  if( sum(testArrestData) < 0 ) return(NULL)
  # To low arrest detectable
  testArrest <- length(testArrestData[testArrestData >= RNAMOD_M7G_ARREST_RATE])
  if( length(testArrestData) != testArrest ) {
    return(NULL)
  }
  # get test values
  # overall mean and sd
  mean <-  mean(baseData)
  sd <-  stats::sd(baseData)
  # Use the sigma level as value for signal strength
  sig <- (as.numeric(as.character(testData)) - mean) %/% sd
  sig.mean <- mean(sig)
  sig.sd <- stats::sd(sig)
  # Since normality of distribution can not be assumed use the MWU
  # generate p.value for single position
  p.value <- suppressWarnings(wilcox.test(baseData, testData)$p.value)
  return(list(sig = sig,
              sig.mean = sig.mean,
              sig.sd = sig.sd,
              p.value = p.value,
              n = n))
}

# call yes or nor position
.validate_M7G_pos <- function(sig.threshold, 
                              p.threshold, 
                              sig, 
                              p.value){
  ((sig > sig.threshold &&
      p.value <= p.threshold) ||
     (sig > sig.threshold &&
        !.get_use_p()))
}

# create data.frame for sample plotting
.create_M7G_plot_data <- function(testData, baseData, name){
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
                        fafile = "FaFile",
                        data = "list"),
  definition = function(object,
                        gff,
                        fafile,
                        data) {
    return(NA)
  }
)
