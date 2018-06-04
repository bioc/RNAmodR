#' @include class-RNAmodR-ident.R
NULL

#' @rdname RNAmodR-ident-class
#'
#' @description 
#' \code{RNAmodRident_RiboMeth}:
#' 
#' @export
setClass("RNAmodRident_RiboMeth",
         contains = "RNAmodRident",
         prototype = list(dataType = "methends",
                          modType = "RiboMeth",
                          param = list(nucleotide = "N")
                          )
)

#' @rdname scoreModificationsPerTranscript
#' 
#' @description 
#' \code{RNAmodRident_RiboMeth}
#' 
#' @export
setMethod(
  f = "scoreModificationsPerTranscript",
  signature = signature(x = "RNAmodRident_RiboMeth",
                        locations = "numeric",
                        data = "GRangesList",
                        args = "RNAmodRargs",
                        endLocation = "numeric"),
  definition = function(x,
                        locations,
                        data,
                        args,
                        endLocation) {
    # setup some data
    locations <- locations[locations > 2 & locations < (max(locations) - 2)]
    weights <- as.numeric(strsplit(getParam(args, "RiboMeth", "weights"),";")[[1]])
    names(weights) <- as.numeric(strsplit(getParam(args, "RiboMeth", "weights_rel_pos"),";")[[1]])
    # calculate scores
    n <- length(data)
    scorea <- lapply(data,
                     .aggregate_score_A,
                     locations,
                     weights)
    scoreb <- lapply(data,
                     .aggregate_score_B,
                     locations,
                     weights)
    scoremeth <- lapply(data,
                        .aggregate_score_meth,
                        locations,
                        weights)
    scoremean <- lapply(data,
                        .aggregate_score_mean,
                        locations,
                        as.numeric(getParam(args, "RiboMeth", "windows_size")))
    # return GPos vor valid RiboMeth positions
    gpos <- data[[1]][pos(data[[1]]) %in% locations]
    mcols(gpos)$counts <- rowMeans(data.frame(lapply(data, function(d){d[pos(d) %in% locations]$counts})))
    mcols(gpos)$Mod <- "RiboMeth"
    mcols(gpos)$n <- n
    mcols(gpos)$scorea <- rowMeans(data.frame(scorea))
    mcols(gpos)$scorea.sd <- matrixStats::rowSds(as.matrix(data.frame(scorea)))
    mcols(gpos)$scoreb <- rowMeans(data.frame(scoreb))
    mcols(gpos)$scoreb.sd <- matrixStats::rowSds(as.matrix(data.frame(scoreb)))
    mcols(gpos)$scoremeth <- rowMeans(data.frame(scoremeth))
    mcols(gpos)$scoremeth.sd <- matrixStats::rowSds(as.matrix(data.frame(scoremeth)))
    mcols(gpos)$scoremean <- rowMeans(data.frame(scoremean))
    mcols(gpos)$scoremean.sd <- matrixStats::rowSds(as.matrix(data.frame(scoremean)))
    return(gpos)
  }
)

# RiboMeth scores --------------------------------------------------------------

# calculates score A according to Birkedal et al. 2015
.aggregate_score_A <- function(data,
                               locations,
                               weights){
  min <- min(as.numeric(names(weights)))
  max <- max(as.numeric(names(weights)))
  # get data
  locationData <- data[pos(data) %in% locations]$counts
  area5 <- lapply(locations, function(l){data[pos(data) < l & pos(data) >= (l + min)]$counts})
  area3 <- lapply(locations, function(l){data[pos(data) <= (l + max) & pos(data) > l]$counts})
  # calculate means and sd
  mean5 <- lapply(area5, mean)
  mean3 <- lapply(area3, mean)
  sd5 <- lapply(area5, stats::sd)
  sd3 <- lapply(area3, stats::sd)
  # calc score per replicate
  scores <- mapply(FUN = .calculate_ribometh_score_A,
                   locationData,
                   mean5,
                   sd5,
                   mean3,
                   sd3)
  return(scores)
}

# calculates score B according to Birkedal et al. 2015
.aggregate_score_B <- function(data,
                               locations,
                               weights){
  min <- min(as.numeric(names(weights)))
  max <- max(as.numeric(names(weights)))
  # get data
  locationData <- data[pos(data) %in% locations]$counts
  area5 <- lapply(locations, function(l){data[pos(data) < l & pos(data) >= (l + min)]$counts})
  area3 <- lapply(locations, function(l){data[pos(data) <= (l + max) & pos(data) > l]$counts})
  # calc score per replicate
  scores <- mapply(FUN = .calculate_ribometh_score_B,
                   locationData,
                   area5,
                   area3,
                   MoreArgs = list(weights))
  return(scores)
}

# calculates score C according to Birkedal et al. 2015
.aggregate_score_meth <- function(data,
                                  locations,
                                  weights){
  min <- min(as.numeric(names(weights)))
  max <- max(as.numeric(names(weights)))
  # get data
  locationData <- data[pos(data) %in% locations]$counts
  area5 <- lapply(locations, function(l){data[pos(data) < l & pos(data) >= (l + min)]$counts})
  area3 <- lapply(locations, function(l){data[pos(data) <= (l + max) & pos(data) > l]$counts})
  # calc score per replicate
  scores <- mapply(FUN = .calculate_ribometh_score_meth,
                   locationData,
                   area5,
                   area3,
                   MoreArgs = list(weights))
  return(scores)
}

# calculates score C according to Birkedal et al. 2014
.aggregate_score_mean <- function(data,
                                  locations,
                                  windowSize){
  locationData5 <- 
    lapply(locations, 
           function(l){
             pos5 <- data[pos(data) == (l - 1)]$counts
             pos <- data[pos(data) == (l)]$counts
             return(1 - pos/pos5)
           })
  locationData3 <- 
    lapply(locations, 
           function(l){
             pos3 <- data[pos(data) == (l + 1)]$counts
             pos <- data[pos(data) == (l)]$counts
             return(1 - pos/pos3)
           })
  area5 <- 
    lapply(locations, 
           function(l){
             locs <- (l - windowSize):(l + windowSize)
             locs <- locs[locs > 2 & locs != l & locs < max(locations)]
             pos5 <- data[pos(data) %in% (locs - 1) ]$counts
             pos <- data[pos(data) %in% locs]$counts
             return(1 - pos/pos5)
           })
  area3 <- 
    lapply(locations, 
           function(l){
             locs <- (l - windowSize):(l + windowSize)
             locs <- locs[locs > 2 & locs != l & locs < max(locations)]
             pos3 <- data[pos(data) %in% (locs + 1) ]$counts
             pos <- data[pos(data) %in% locs]$counts
             return(1 - pos/pos3)
           })
  mean5 <- lapply(area5, mean)
  mean3 <- lapply(area3, mean)
  #
  locationData5[is.nan(unlist(locationData5)) | is.infinite(unlist(locationData5))] <- c(0)
  locationData3[is.nan(unlist(locationData3)) | is.infinite(unlist(locationData3))] <- c(0)
  mean5[is.nan(unlist(mean5)) | is.infinite(unlist(mean5))] <- c(1)
  mean3[is.nan(unlist(mean3)) | is.infinite(unlist(mean3))] <- c(1)
  # calc score per replicate
  scores <- mapply(FUN = .calculate_ribometh_score_mean,
                   locationData5,
                   locationData3,
                   mean5,
                   mean3)
  return(scores)
}


#' @rdname identifyModificationsPerTranscript
#' 
#' @description 
#' \code{RNAmodRident_RiboMeth}
#' 
#' @export
setMethod(
  f = "identifyModificationsPerTranscript",
  signature = signature(x = "RNAmodRident_RiboMeth",
                        data = "GPos",
                        args = "RNAmodRargs"),
  definition = function(x,
                        data,
                        args) {
    browser()
    
    
    return(return(gpos))
  }
)