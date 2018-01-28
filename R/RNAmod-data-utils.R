#' @include class-RNAmod-analysis-type.R
NULL

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
    # infinite in one direction
    if(is.na(b) || b == 0) return(-1)
    # max = 1 in one direction
    if(a <= b) return((1-a/b))
    # infinite in other direction
    if(is.na(a) || a == 0) return(1)
    # min = -1 in other direction
    return((b/a-1))
  }))
  setNames(y,names(x))
}

# Calculate the arrest difference per position
.get_arrest_diff <- function(i, meanData){
  x <- data[[i]]
  y <- unlist(lapply(seq_along(x), function(j){
    x[as.numeric(names(x)) == j]-meanData[[i]][as.numeric(names(x)) == j]
  }))
  setNames(y,names(x))
}

# calculates a rolling mean values
.get_rolling_mean <- function(x,
                              n=3){
  y <- setNames(filter(x,rep(1/n,n), sides=2),names(x))
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
  unlist(lapply(data,function(dataPerReplicate){
    return(dataPerReplicate[as.numeric(names(dataPerReplicate)) < (location+width) &
                              as.numeric(names(dataPerReplicate)) > (location-width) &
                              as.numeric(names(dataPerReplicate)) != location])
  }))
}