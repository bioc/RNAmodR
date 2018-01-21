#' @include class-RNAmod-mod-type.R
NULL

RNAMOD_D_NUCLEOTIDE <- "T"
RNAMOD_D_ROLLING_MEAN_WINDOW_WIDTH <- 9
RNAMOD_D_ARREST_RATE <- 0.95
RNAMOD_D_P_THRESHOLD <- 0.05
RNAMOD_D_SIGMA_THRESHOLD <- 3


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

#' @rdname maskPositionData
#' 
#' @description 
#' \code{mod_D}
#' 
#' @return
#' @export
#'
#' @examples
setMethod(
  f = "maskPositionData",
  signature = signature(object = "mod_D",
                        data = "numeric",
                        modLocations = "numeric"),
  definition = function(object,
                        data,
                        modLocations) {
    data[as.numeric(names(data)) %in% (modLocations+1)] <- 
      data[as.numeric(names(data)) %in% (modLocations+1)] * 
      (1-RNAMOD_D_ARREST_RATE)
    return(data)
  }
)

#' @rdname preTest
#' 
#' @description 
#' \code{mod_D}
#' 
#' @return
#' @export
#' 
#' @importFrom stringr str_locate_all
#' @importFrom BiocParallel bplapply
#'
#' @examples
setMethod(
  f = "preTest",
  signature = signature(object = "mod_D",
                        location = "numeric",
                        data = "list",
                        globalLocations = "numeric"),
  definition = function(object,
                        location,
                        data,
                        globalLocations) {
    # if non G position skip position
    if( names(globalLocations[globalLocations == location]) 
        != RNAMOD_D_NUCLEOTIDE){
      return(NULL)
    }
    # do pretest
    res <- .do_D_pretest(location,
                           data)
    return(res)
  }
)

# check if any data is available to proceed with test
# this is in a seperate function since it is also called by checkForModification
.do_D_pretest <- function(location,
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
  # get test values
  # overall mean and sd
  mean <-  mean(baseData)
  sd <-  stats::sd(baseData)
  # Use the sigma level as value for signal strength
  if( mean((as.numeric(as.character(testData)) - mean) %/% sd) 
      <= RNAMOD_D_SIGMA_THRESHOLD) {
    return(NULL)
  }
  return(list(n = n,
              testData = testData,
              baseData = baseData))
}

#' @rdname checkForModification
#' 
#' @description 
#' \code{mod_D}
#' 
#' @return
#' @export
#' 
#' @importFrom stringr str_locate_all
#' @importFrom BiocParallel bplapply
#'
#' @examples
setMethod(
  f = "checkForModification",
  signature = signature(object = "mod_D",
                        location = "numeric",
                        globalLocations = "numeric",
                        data = "list",
                        modClasses = "list"),
  definition = function(object,
                        location,
                        globalLocations,
                        data,
                        modClasses) {
    # browser()
    # get test result for the current location
    locTest <- .calc_D_test_values(location,
                                     data)
    # If insufficient data is present
    if(is.null(locTest)) return(NULL)
    # dynamic threshold based on the noise of the signal (high sd)
    if(!.validate_D_pos(RNAMOD_D_SIGMA_THRESHOLD, 
                          RNAMOD_D_P_THRESHOLD, 
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

# test for D at current location
.calc_D_test_values <- function(location,
                                  data){
  # short cut if amount of data is not sufficient
  pretestData <- .do_D_pretest(location,
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
  testArrest <- length(testArrestData[testArrestData >= RNAMOD_D_ARREST_RATE])
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
.validate_D_pos <- function(sig.threshold, 
                              p.threshold, 
                              sig, 
                              p.value){
  ((sig > sig.threshold &&
      p.value <= p.threshold) ||
     (sig > sig.threshold &&
        !.get_use_p()))
}

# create data.frame for sample plotting
.create_D_plot_data <- function(testData, baseData, name){
  data.frame(x = c(rep(name,(length(testData) + 
                               length(baseData)))),
             y = c(testData, 
                   baseData),
             group = c(rep("Position",length(testData)), 
                       rep("Baseline", length(baseData))))
}
