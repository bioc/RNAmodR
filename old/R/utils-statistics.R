#' @include RNAmodR.R
#' 



# RiboMeth scores --------------------------------------------------------------

# calculates score A according to Birkedal et al. 2015
.calculate_ribometh_score_A_c <- function(n,
                                        meanL,
                                        sdL,
                                        meanR,
                                        sdR){
  if(is.nan(meanL) || is.na(sdL) || is.na(meanR) || is.na(sdR)){
    return(NA)
  }
  dividend <- (2 * n  + 1)
  divisor <- (0.5 * abs(meanL - sdL) + n + 
                0.5 * abs(meanR - sdR) + 1 ) 
  # return(max(0, (1 - (dividend / divisor))))
  return((1 - (dividend / divisor)))
}
.calculate_ribometh_score_A <- compiler::cmpfun(.calculate_ribometh_score_A_c)

.calculate_ribometh_score_B_c <- function(n,
                                        areaL,
                                        areaR,
                                        weights){
  # split the weights
  weightsL <- weights[names(weights) < 0]
  weightsR <- weights[names(weights) > 0]
  if(length(weightsL) != length(areaL) || length(weightsR) != length(areaR)){
    return(NA)
  }
  # calculates score
  dividend <- abs(n - 0.5 * (stats::weighted.mean(areaL, weightsL) +
                              stats::weighted.mean(areaR, weightsR)))
  divisor <- (n + 1)
  return(dividend / divisor)
}
.calculate_ribometh_score_B <- compiler::cmpfun(.calculate_ribometh_score_B_c)

# calculates score C according to Birkedal et al. 2014
.calculate_ribometh_score_meth_c <- function(n,
                                           areaL,
                                           areaR,
                                           weights){
  # split the weights
  weightsL <- weights[names(weights) < 0]
  weightsR <- weights[names(weights) > 0]
  if(length(weightsL) != length(areaL) || length(weightsR) != length(areaR)){
    return(NA)
  }
  # calculates score
  return((1 - (n / 
                ( 0.5 * (stats::weighted.mean(areaL, weightsL) +
                           stats::weighted.mean(areaR, weightsR))
                  )
              )
         )
      )
}
.calculate_ribometh_score_meth <- compiler::cmpfun(.calculate_ribometh_score_meth_c)

# calculates score mean according to Marchand et al. 2016
.calculate_ribometh_score_mean_c <- function(locationData5,
                                             locationData3,
                                             mean5,
                                             mean3,
                                             sd5,
                                             sd3){
  if(length(locationData5) == 0 || length(locationData3) == 0){
    return(NA)
  }
  value5 <- (locationData5 - mean5)/sd5
  value3 <- (locationData3 - mean3)/sd3
  valueMean <- mean(c(value5,value3))
  return(max(abs(value5 - valueMean),abs(value3 - valueMean)))
}
.calculate_ribometh_score_mean <- compiler::cmpfun(.calculate_ribometh_score_mean_c)

# Matthew's correlation coefficient --------------------------------------------

.get_mcc_c <- function(tp,
                     tn,
                     fp,
                     fn){
  mcc <- ( (tp * tn) - (fp * fn) ) / 
    sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc
}
.get_mcc <- compiler::cmpfun(.get_mcc_c)

