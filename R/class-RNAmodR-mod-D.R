#' @include class-RNAmodR-mod-type.R
NULL

RNAMODR_D_NUCLEOTIDE <- "T"
RNAMODR_D_ARREST_RATE <- 0.6
RNAMODR_D_ARREST_RATE_INV <- 1 - RNAMODR_D_ARREST_RATE
RNAMODR_D_Z_THRESHOLD <- 3
RNAMODR_D_SIG_THRESHOLD <- 5


#' @rdname RNAmodR-mod-class
#'
#' @description 
#' \code{mod_D}:
#' 
#' @export
setClass("mod_D",
         contains = "mod",
         prototype = list(dataType = "5end",
                          modType = "D",
                          positionOffset = 1)
)

#' @rdname checkForModification
#' 
#' @description 
#' \code{mod_D}
setMethod(
  f = "checkForModification",
  signature = signature(x = "mod_D",
                        location = "numeric",
                        locations = "numeric",
                        data = "list"),
  definition = function(x,
                        location,
                        locations,
                        data) {
    # get test result for the current location
    locTest <- .calc_D_test_values(location,
                                   locations,
                                   data)
    # If insufficient data is present
    if(is.null(locTest)) return(NULL)
    # dynamic threshold based on the noise of the signal (high sd)
    if(!.validate_D_pos(RNAMODR_D_SIG_THRESHOLD, 
                        RNAMODR_D_Z_THRESHOLD, 
                        locTest$sig.mean, 
                        locTest$z) ) {
      return(NULL)
    }
    # Return data
    return(list(location = location,
                type = getModType(x),
                signal = locTest$sig.mean,
                signal.sd = locTest$sig.sd,
                z = locTest$z,
                nbsamples = locTest$n))
  }
)

# test for D at current location
.calc_D_test_values <- function(location,
                                locations,
                                data){
  # split data
  conditions <- names(data)
  coverage <- lapply(data,"[[","coverage")
  names(coverage) <- conditions
  data <- lapply(data,"[[","data")
  names(data) <- conditions
  
  # if(location == 20) browser()
  # short cut if amount of data is not sufficient
  pretestData <- .do_D_pretest(location,
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
    if(any(testDataCorrected < (RNAMODR_D_ARREST_RATE - testDataControl))) return(NULL)
    # difference to threshold by relative percent
    # minimal value of 1 since y is one percent lower than arrest rate threshold
    sig <- floor((testDataCorrected - (RNAMODR_D_ARREST_RATE - testDataControl)) / 
                   (1 - (RNAMODR_D_ARREST_RATE - testDataControl)) * 100)
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
    if(any(testData < RNAMODR_D_ARREST_RATE)) return(NULL)
    # difference to threshold by relative percent
    # minimal value of 1 since y is one percent lower than arrest rate threshold
    sig <- floor((testData - RNAMODR_D_ARREST_RATE) / 
                   (1 - RNAMODR_D_ARREST_RATE) * 100)
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

# check if any data is available to proceed with test
# this is in a seperate function since it is also called by checkForModification
.do_D_pretest <- function(location,
                          locations,
                          data,
                          coverage){
  # if non D position skip position
  if( names(locations[locations == location]) != RNAMODR_D_NUCLEOTIDE){
    return(NULL)
  }
  # do not take into account position 1 or the other end
  if(location == 1 || location > (max(locations) - 5) ) return(NULL)
  # split into conditions
  res <- mapply(FUN = .get_data_per_condition_D,
                split(data,names(data)),
                split(coverage,names(coverage)),
                MoreArgs = list(location = location),
                SIMPLIFY = FALSE)
  res <- res[!vapply(res,is.null,logical(1))]
  # check if Treated condition has valid data. Otherwise return NULL (Abort)
  # check if treated data ends at N+1 
  if(is.null(res$Treated) || any(is.na(res$Treated$posData))){
    return(NULL)
  }
  return(res)
}

.get_data_per_condition_D <- function(data, 
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
  if(any( (unlist(testCoverage) * RNAMODR_D_ARREST_RATE_INV) < RNAMODR_5END_COVERAGE_MIN )) return(NULL)
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

# iterates on every position and calculates the difference of the means
.merge_base_data_D <- function(treated,
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
.validate_D_pos <- function(sig.threshold,
                            p.threshold,
                            sig,
                            z){
  ((sig > sig.threshold &&
      z >= p.threshold) ||
     (sig > sig.threshold &&
        !.get_use_p()))
}