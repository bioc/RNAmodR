#' @include class-RNAmod-mod-type.R
NULL

RNAMOD_M7G_NUCLEOTIDE <- "G"
RNAMOD_M7G_SPM <- 3
RNAMOD_M7G_ARREST_RATE <- 0.95
RNAMOD_M7G_P_THRESHOLD <- 0.05
RNAMOD_M7G_SIGMA_THRESHOLD <- 15


#' @rdname mod
#'
#' @description 
#' \code{mod_m7G}
#' 
#' @export
setClass("mod_m7G",
         contains = "mod",
         prototype = list(modType = "m7G")
)

#' @rdname maskPositionData
#' 
#' @description 
#' \code{mod_m7G}
#' 
#' @export
setMethod(
  f = "maskPositionData",
  signature = signature(object = "mod_m7G",
                        data = "numeric",
                        modLocations = "numeric"),
  definition = function(object,
                        data,
                        modLocations) {
    data[as.numeric(names(data)) %in% (modLocations+1)] <- 
      data[as.numeric(names(data)) %in% (modLocations+1)] * 
      (1-RNAMOD_M7G_ARREST_RATE)
    return(data)
  }
)

#' @rdname preTest
#' 
#' @description 
#' \code{mod_m7G}
#' 
#' @export
setMethod(
  f = "preTest",
  signature = signature(object = "mod_m7G",
                        location = "numeric",
                        data = "list",
                        locations = "numeric"),
  definition = function(object,
                        location,
                        data,
                        locations) {
    # do pretest
    res <- .do_M7G_pretest(location,
                           locations,
                           data)
    return(res)
  }
)

# check if any data is available to proceed with test
# this is in a seperate function since it is also called by checkForModification
.do_M7G_pretest <- function(location,
                            locations,
                            data){
  # if non G position skip position
  if( names(locations[locations == location]) != RNAMOD_M7G_NUCLEOTIDE){
    return(NULL)
  }
  # do not take into account position 1
  if(location == 1) return(NULL)
  # browser()
  # number of replicates
  n <- length(data)
  # merge data for positions
  # data on the N+1 location
  testData <- .aggregate_location_data(data, (location+1))
  testData <- testData[testData > 0]
  # if spm is not high enough
  if(length(testData[testData >= RNAMOD_M7G_SPM]) < n) return(NULL)
  # base data to compare against
  baseData <- .aggregate_not_location_data(data, (location+1))
  # if not enough data is present
  if(length(testData) < n | 
     length(baseData) < (3*n)) return(NULL)
  # get test values
  # overall mean and sd
  mean <-  mean(baseData)
  sd <-  stats::sd(baseData)
  # Use the sigma level as value for signal strength
  if( mean((as.numeric(as.character(testData)) - mean) %/% (mean+sd)) 
      <= RNAMOD_M7G_SIGMA_THRESHOLD) {
    return(NULL)
  }
  return(list(n = n,
              testData = testData,
              baseData = baseData))
}

#' @rdname checkForModification
#' 
#' @description 
#' \code{mod_m7G}
#' 
#' @export
setMethod(
  f = "checkForModification",
  signature = signature(object = "mod_m7G",
                        location = "numeric",
                        locations = "numeric",
                        data = "list"),
  definition = function(object,
                        location,
                        locations,
                        data) {
    # if(location == 1575 | location == 599 | location == 1420) { browser() }
    # if(location == 1925) { browser() }
    # browser()
    # get test result for the current location
    locTest <- .calc_M7G_test_values(location,
                                     locations,
                                     data)
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
                type = getModType(object),
                signal = locTest$sig.mean,
                signal.sd = locTest$sig.sd,
                p.value = locTest$p.value,
                nbsamples = locTest$n))
  }
)

# test for m7G at current location
.calc_M7G_test_values <- function(location,
                                  locations,
                                  data){
  # short cut if amount of data is not sufficient
  pretestData <- .do_M7G_pretest(location,
                                 locations,
                                 data)
  if(is.null(pretestData)) return(NULL)
  # data from pretest
  testData <- pretestData$testData
  baseData <- pretestData$baseData
  n <- pretestData$n
  # Calculate the arrest rate per position
  arrestData <- lapply(data, .get_arrest_rate)
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
  sig <- (as.numeric(as.character(testData)) - mean) %/% (mean+sd)
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
