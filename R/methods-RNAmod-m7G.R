#' @include RNAmod.R
NULL

RNAMOD_M7G_P_THRESHOLD <- 0.001
RNAMOD_M7G_SIGMA_THRESHOLD <- 5


#' @rdname parseMod
#' 
#' @description 
#' \code{mod_m7G}
#' 
#' @return
#' @export
#' 
#' @importFrom stringr str_locate_all
#'
#' @examples
setMethod(
  f = "parseMod",
  signature = signature(object = "mod_m7G",
                        counts = "numeric",
                        gff = "GRanges",
                        seq = "DNAString",
                        data = "DataFrame"),
  definition = function(object,
                        counts,
                        gff,
                        seq,
                        data) {
    
    strand <- as.character(strand(gff))
    loc <- stringr::str_locate_all(as.character(seq), "G")
    loc <- loc[[1]][,"start"]
    if(strand == "-"){
      pos <- data$pos + data$qwidth - 1
      locations <- (end(gff) - loc) + 1
    } else {
      pos <- data$pos
      locations <- (start(gff) + loc) - 1
    }
    # Normalize counts per positions against million of reads in BamFile
    posData <- table(pos)/(counts/10^6)
    
    sd <-  stats::sd(posData)
    mean <-  mean(posData)
    
    res <- lapply(locations,
                  .check_for_m7G,
                  posData,
                  strand,
                  mean,
                  sd)
    names(res) <- paste0("G_",loc)
    res <-  res[!vapply(res,is.null,logical(1))]
    
    return(res)
  }
)


.check_for_m7G <- function(location, posData, strand, mean, sd){
  # if(location < 456160){
  #   browser()
  # }
  subset <- posData[names(posData) < (location+50) & names(posData) > (location-50)]
  
  # Fill missing positions with 0?
  
  # m7G expects a high number of reads at the N+1 position
  if(as.character(strand) == "-"){
    loc <- location-1
  } else {
    loc <- location+1
  }
  bool <- names(subset) == loc
  testData <- subset[bool]
  baseData <- subset[!bool]
  
  if(length(testData) == 0){
    return(NULL)
  }
  if(length(baseData) < 3){
    return(NULL)
  }
  
  # Since normality of distribution can not be assumed use the MWU
  # Does not work very well since the distribution is half normal
  # test <- stats::wilcox.test(testData, baseData)
  # 
  # if( test$p.value < RNAMOD_M7G_P_THRESHOLD  ){
  #   return(list(location = location,
  #               value = as.numeric(as.character(testData)),
  #               p.value = test$p.value))
  # }
  
  # Use the sigma level as value for signal strength
  sigStrength <- (as.numeric(as.character(testData)) - mean) %/% sd
  # generate p.value
  p.value <- 0.01
  
  if( sigStrength > RNAMOD_M7G_SIGMA_THRESHOLD  ){
    return(list(location = location,
                value = as.numeric(as.character(testData)),
                signal = sigStrength,
                p.value = p.value))
  }
  return(NULL)
}


#' @rdname mergeModsOfReplicates
#'
#' @description
#' \code{mod_m7G}
#'
#' @return
#' @export
#'
#' @examples
setMethod(
  f = "mergeModsOfReplicates",
  signature = signature(object = "mod_m7G",
                        gff = "GRanges",
                        seq = "FaFile",
                        data = "list"),
  definition = function(object,
                        gff,
                        seq,
                        data) {
    # Process only genes found in all datasets
    IDs <- lapply(data,names)
    IDs <- Reduce(intersect, IDs)
    res <- lapply(IDs,
                  FUN = .merge_m7G_signal,
                  data)
    names(res) <- IDs
    res <- res[!is.null(res)]
    return(res)
  }
)

# merge signals for transcripts
.merge_m7G_signal <- function(ID,data){
  data <- lapply(data,function(x,id){
    return(x[[id]][["m7G"]])
  },ID)
  
  positions <- unique(unlist(lapply(data,names)))
  res <- lapply(positions,
                FUN = .merge_m7G_mod,
                data)
  names(res) <- positions
  res <- res[!is.null(res)]
  
  return(res)
}

# merge signals for one m7G position
.merge_m7G_mod <- function(pos,data){
  data <- lapply(data, "[[",pos)
  # Number of replicates not having the modified positions
  undetected <- length(data[is.null(data)])
  data <- data[!is.null(data)]
  if( length(unique(lapply(data,"[[","location"))) != 1 ){
    stop("Something went wrong. Locations for gene '",ID,"' do not match up.",
         call. = FALSE)
  }
  return(list(location = unique(unlist(lapply(data,"[[","location"))),
              value = lapply(data,"[[","value"),
              signal = mean(unlist(lapply(data,"[[","signal"))),
              signal.sd = stats::sd(unlist(lapply(data,"[[","signal"))),
              p.value = mean(unlist(lapply(data,"[[","p.value"))),
              p.value.sd = stats::sd(unlist(lapply(data,"[[","p.value"))),
              replicates = (length(data)+undetected),
              undetectedInNReplicate = undetected))
}


#' @rdname convertReadsToPositions
#'
#' @description
#' \code{mod_m7G}
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
                        seq = "DNAString",
                        data = "DataFrame"),
  definition = function(object,
                        counts,
                        gff,
                        seq,
                        data) {
    strand <- as.character(strand(gff))
    if(strand == "-"){
      pos <- data$pos + data$qwidth - 1
    } else {
      pos <- data$pos
    }
    # Normalize counts per positions against million of reads in BamFile
    posData <- table(pos)/(counts/10^6)

    return(posData)
  }
)



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
                        seq = "FaFile",
                        data = "list"),
  definition = function(object,
                        gff,
                        seq,
                        data) {
    # Process only genes found in all datasets
    IDs <- lapply(data,names)
    IDs <- Reduce(intersect, IDs)
    res <- lapply(IDs,
                  FUN = .merge_m7G_positions,
                  data)
    names(res) <- IDs
    res <- res[!is.null(res)]
    return(res)
  }
)

# merge positions in one transcript
.merge_m7G_positions <- function(ID,data){
  data <- lapply(data,function(x,id){
    return(x[[id]][["m7G"]])
  },ID)
  
  positions <- unique(unlist(lapply(data,names)))
  res <- lapply(positions,
                FUN = .merge_m7G_position,
                data)
  df <- data.frame(pos = unlist(lapply(res,"[[","pos")),
                   mean = unlist(lapply(res,"[[","mean")),
                   sd = unlist(lapply(res,"[[","sd")))
  return(df)
}

# merge position data for one position
.merge_m7G_position <- function(pos,data){
  data <- unlist(lapply(data, "[[",pos))
  data[is.na(data)] <- 0
  return(list(pos = pos,
              mean = mean(data),
              sd = stats::sd(data)))
}

