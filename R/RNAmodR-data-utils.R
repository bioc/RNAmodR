#' @include class-RNAmodR-analysis-type.R
NULL

# modificaton access for SE and GRL --------------------------------------------

.extract_modification_info_from_grl <- function(grl){
  res <- lapply(grl, function(gr){
    mcols <- S4Vectors::mcols(gr)
    DF <- cbind(S4Vectors::DataFrame(chrom = GenomeInfoDb::seqnames(gr),
                                     start = BiocGenerics::start(gr),
                                     end = BiocGenerics::end(gr),
                                     strand = as.character(BiocGenerics::strand(gr))),
                mcols[,!(colnames(mcols) %in% c("phase"))])
    DF$RNAmodR_signal <- as.numeric(DF$RNAmodR_signal)
    DF$RNAmodR_signal_sd <- as.numeric(DF$RNAmodR_signal_sd)
    DF$RNAmodR_p.value <- as.numeric(DF$RNAmodR_p.value)
    DF$RNAmodR_nbReplicates <- as.numeric(DF$RNAmodR_nbReplicates)
    rownames(DF) <- DF$ID
    DF
  })
  names(res) <- names(grl)
  res
}

.extract_modification_info_from_se <- function(ses){
  res <- lapply(ses, function(se){
    DF <- do.call(rbind,
                  seGetModifications(se))
    DF$RNAmodR_signal <- as.numeric(DF$RNAmodR_signal)
    DF$RNAmodR_signal_sd <- as.numeric(DF$RNAmodR_signal_sd)
    DF$RNAmodR_p.value <- as.numeric(DF$RNAmodR_p.value)
    DF$RNAmodR_nbReplicates <- as.numeric(DF$RNAmodR_nbReplicates)
    DF
  })
  names(res) <- names(ses)
  res
}


# common function for converting data ------------------------------------------

# Calculate FPKM values
.calculate_fpkm <- function(nReads,length,R){
  # Calculate fpkm
  fpkm <- nReads/((length/1000) * (R/1000000))
  return(fpkm)
}

# Calculate the arrest rate per position
.get_arrest_rate <- function(x){
  if(is.null(names(x))) stop("Unnamed position data.")
  y <- unlist(lapply(seq_along(x), function(i){
    a <- x[i]
    b <- x[(i+1)]
    if(is.na(b) || b == 0) {
      b <- 1
    }
    if(is.na(a) || a == 0) {
      a <- 1
    }
    # max = 1 in one direction
    if(a <= b) return((1-a/b))
    # min = -1 in other direction
    return((b/a-1))
  }))
  stats::setNames(y,names(x))
}

# calculates a rolling mean values
.get_rolling_mean <- function(x,
                              n=3){
  y <- stats::setNames(stats::filter(x,rep(1/n,n), sides=2),names(x))
  y[is.na(y)] <- x[is.na(y)]
  return(y)
}

# common function for subsetting data ------------------------------------------

.aggregate_location_data <- function(data, 
                                     location){
  unlist(lapply(data,function(dataPerReplicate){
    return(dataPerReplicate[as.numeric(names(dataPerReplicate)) == location])
  }))
}
.aggregate_not_location_data <- function(data,
                                         location){
  unlist(lapply(data,function(dataPerReplicate){
    return(dataPerReplicate[as.numeric(names(dataPerReplicate)) != location])
  }))
}
.aggregate_area_data <- function(data, 
                                 location, 
                                 width){
  unlist(.subset_area_data(data, 
                           location, 
                           width))
}
.subset_area_data <- function(data, 
                              location, 
                              width){
  lapply(data,function(dataPerReplicate){
    return(dataPerReplicate[as.numeric(names(dataPerReplicate)) < (location + width) &
                              as.numeric(names(dataPerReplicate)) > (location - width) &
                              as.numeric(names(dataPerReplicate)) != location])
  })
}
.subset_area_data2 <- function(data, 
                               location, 
                               width){
  lapply(data,function(dataPerReplicate){
    return(dataPerReplicate[as.numeric(names(dataPerReplicate)) < (location + width) &
                              as.numeric(names(dataPerReplicate)) > (location - width)])
  })
}