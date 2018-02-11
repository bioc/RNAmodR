#' @include class-RNAmodR-mod-type.R
NULL

RNAMODR_M7G_NUCLEOTIDE <- "G"
RNAMODR_M7G_ARREST_RATE <- 0.65
RNAMODR_M7G_ARREST_RATE_INV <- 1 - RNAMODR_M7G_ARREST_RATE
RNAMODR_M7G_Z_THRESHOLD <- 3
RNAMODR_M7G_SIG_THRESHOLD <- 5


#' @rdname RNAmodR-mod-class
#'
#' @description 
#' \code{mod_m7G}:
#' 
#' @export
setClass("mod_m7G",
         contains = "mod",
         prototype = list(modType = "m7G",
                          positionOffset = 1)
)

# check if any data is available to proceed with test
# this is in a seperate function since it is also called by checkForModification
.do_M7G_pretest <- function(location,
                            locations,
                            data,
                            coverage){
  # if non G position skip position
  if( names(locations[locations == location]) != RNAMODR_M7G_NUCLEOTIDE){
    return(NULL)
  }
  # do not take into account position 1 or the other end
  if(location == 1 || location > (max(locations) - 5) ) return(NULL)
  # split into conditions
  res <- mapply(FUN = .get_data_per_condition_M7G,
                split(data,names(data)),
                split(coverage,names(coverage)),
                MoreArgs = list(location = location),
                SIMPLIFY = FALSE)
  res <- res[!vapply(res,is.null,logical(1))]
  # check if Treated condition has valid data. Otherwise return NULL (Abort)
  if(is.null(res$Treated)){
    return(NULL)
  }
  return(res)
}

.get_data_per_condition_M7G <- function(data,
                                        coverage,
                                        location){
  # number of replicates
  n <- length(data)
  # data on the N location
  posData <- .aggregate_location_data(data,
                                      location)
  # data on the N+1 location
  testData <- .aggregate_location_data(data, 
                                       (location + 1))
  # if number of data points is not high enough or their are empty
  if(any(is.na(testData)) || length(unlist(testData)) == 0 ) return(NULL)
  # if minimal arrest rate requires to low coverage
  testCoverage <- .aggregate_location_data(coverage, 
                                           (location + 1))
  if(any( (unlist(testCoverage) * RNAMODR_M7G_ARREST_RATE_INV) < RNAMODR_DEFAULT_COVERAGE_MIN )) return(NULL)
  # base data to compare against
  baseData <- .aggregate_area_data(data, 
                                   (location + 1), 
                                   50)
  # if not enough data is present
  if(sum(vapply(lapply(baseData, 
                       function(x){x[!is.na(x)]}),
                length,
                numeric(1))) < n) {
    return(NULL)
  }
  return(list(n = n,
              posData = posData,
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
    # get test result for the current location
    locTest <- .calc_M7G_test_values(location,
                                     locations,
                                     data)
    # If insufficient data is present
    if(is.null(locTest)) return(NULL)
    # dynamic threshold based on the noise of the signal (high sd)
    if(!.validate_M7G_pos(RNAMODR_M7G_SIG_THRESHOLD, 
                          RNAMODR_M7G_Z_THRESHOLD, 
                          locTest$sig.mean, 
                          locTest$z) ) {
      return(NULL)
    }
    # Return data
    return(list(location = location,
                type = getModType(object),
                signal = locTest$sig.mean,
                signal.sd = locTest$sig.sd,
                z = locTest$z,
                nbsamples = locTest$n))
  }
)

# test for m7G at current location
.calc_M7G_test_values <- function(location,
                                  locations,
                                  data){
  # split data
  conditions <- names(data)
  coverage <- lapply(data,"[[","coverage")
  names(coverage) <- conditions
  data <- lapply(data,"[[","data")
  names(data) <- conditions
  
  # if(location == 1575) browser()
  # short cut if amount of data is not sufficient
  pretestData <- .do_M7G_pretest(location,
                                 locations,
                                 data,
                                 coverage)
  if(is.null(pretestData)) return(NULL)
  
  # If Control sample is available
  if(!is.null(pretestData$Control)){
    # data from pretest
    testData <- pretestData$Treated$testData
    baseData <- pretestData$Treated$baseData
    n <- pretestData$Treated$n
    testDataC <- pretestData$Control$testData
    baseDataC <- pretestData$Control$baseData
    nC <- pretestData$Control$n
    
    testDataControl <- mean(unlist(testDataC))
    testDataCorrected <- unlist(testData) - testDataControl
    # if stop coverage is to low
    if(any(testDataCorrected < (RNAMODR_M7G_ARREST_RATE - testDataControl))) return(NULL)
    # difference to threshold by relative percent
    # minimal value of 1 since y is one percent lower than arrest rate threshold
    sig <- floor((testDataCorrected - (RNAMODR_M7G_ARREST_RATE - testDataControl)) / 
                   (1 - (RNAMODR_M7G_ARREST_RATE - testDataControl)) * 100)
    sig.mean <- mean(sig)
    # approx. sd since covariance of testData and testDataC cannot not be 
    # assumed in case of unequal number of  observations
    sig.sd <- sqrt(stats::sd(unlist(testData))^2 + 
                     stats::sd(unlist(testDataC))^2)
    # generate z score
    baseData <- .merge_base_data_M7G(baseData,
                                     baseDataC)
    # z <- mean((testDataCorrected - mean(baseData, 
    #                                     na.rm = TRUE)) / sd(baseData, 
    #                                                         na.rm = TRUE))
    z <- mean((testDataCorrected - mean(baseData)) / sd(baseData))
  } else {
    # data from pretest
    testData <- unlist(pretestData$Treated$testData)
    baseData <- unlist(pretestData$Treated$baseData)
    n <- pretestData$Treated$n
    # if stop coverage is to low
    if(any(testData < RNAMODR_M7G_ARREST_RATE)) return(NULL)
    # difference to threshold by relative percent
    # minimal value of 1 since y is one percent lower than arrest rate threshold
    sig <- floor((testData - RNAMODR_M7G_ARREST_RATE) / 
                   (1 - RNAMODR_M7G_ARREST_RATE) * 100)
    sig.mean <- mean(sig)
    sig.sd <- stats::sd(sig)
    # generate z score
    z <- mean((testData - mean(baseData, 
                               na.rm = TRUE)) / sd(baseData, 
                                                   na.rm = TRUE))
  }
  return(list(sig = sig,
              sig.mean = sig.mean,
              sig.sd = sig.sd,
              z = z,
              n = n))
}
# iterates on every position and calculates the difference of the means
.merge_base_data_M7G <- function(treated,
                             control){
  df <- data.frame(append(treated,control))
  colnames(df) <- c(names(treated),names(control))
  df <- df[complete.cases(df),]
  treated <- df[,colnames(df) == "Treated"]
  control <- df[,colnames(df) == "Control"]
  res <- rowMeans(treated) - rowMeans(control)
  res[res < 0] <- 0
  return(res)
}

# call yes or nor position
.validate_M7G_pos <- function(sig.threshold, 
                              p.threshold, 
                              sig, 
                              z){
  ((sig > sig.threshold &&
      z >= p.threshold) ||
     (sig > sig.threshold &&
        !.get_use_p()))
}