#' @include RNAmodR.R
#' 



# RiboMeth scores --------------------------------------------------------------

# calculates score A according to Birkedal et al. 2014
.calculate_ribometh_score_A <- function(n,
                               meanL,
                               sdL,
                               meanR,
                               sdR){
  max(0, (1 - ( (2 * n  + 1) / 
                  (0.5 * abs(meanL - sdL) + n + 
                     0.5 * abs(meanR - sdR) + 1 ) 
                )
          )
      )
}

# calculates score B according to Birkedal et al. 2014
.calculate_ribometh_score_B <- function(n,
                               areaL,
                               areaR,
                               weights){
  weightsL <- weights[names(weights) < 0]
  weightsR <- weights[names(weights) > 0]
  x <- abs(n - 0.5 * (stats::weighted.mean(areaL, weightsL) / sum(weightsL) +
                        stats::weighted.mean(areaR, weightsR) / sum(weightsR)) 
           ) /  
    (n + 1)
  x
}

# calculates score C according to Birkedal et al. 2014
.calculate_ribometh_score_C <- function(n,
                               areaL,
                               areaR,
                               weights){
  weightsL <- weights[names(weights) < 0]
  weightsR <- weights[names(weights) > 0]
  max(0,(1 - (n / 
                ( 0.5 * (stats::weighted.mean(areaL, weightsL) / sum(weightsL) +
                           stats::weighted.mean(areaR, weightsR) / sum(weightsR))
                  )
              )
         )
      )
}

# calculates score MAX according to Marchand et al. 2016
.calculate_ribometh_score_MAX <- function(data,
                                 pos){
  
}

# Matthew's correlation coefficient --------------------------------------------

.get_mcc <- function(tp,
                     tn,
                     fp,
                     fn){
  mcc <- ( (tp * tn) - (fp * fn) ) / 
    sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc
}

