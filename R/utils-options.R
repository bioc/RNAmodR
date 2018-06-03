#' @include RNAmodR.R
NULL

# option retrieval -------------------------------------------------------------

# whether to use a p value for detection or not
.get_min_map_quality <- function(){
  adviseText <- "Please set 'RNAmodR_min_mapquality' to a numeric value of >= 0."
  mapQ <- getOption("RNAmodR_min_mapquality")
  if(!assertive::is_a_number(mapQ)){
    mapQ <- RNAMODR_MIN_MAPQUALITY
    warning("The option 'RNAmodR_min_mapquality' is not a single numeric value. ",
            adviseText,
            call. = FALSE)
  }
  if(!(mapQ >= 0)){
    mapQ <- RNAMODR_MIN_MAPQUALITY
    warning("The option 'RNAmodR_min_mapquality' is not 0 or a positive value ",
            adviseText,
            call. = FALSE)
  }
  mapQ
}

# whether to use a p value for detection or not
.get_use_p <- function(){
  useP <- getOption("RNAmodR_use_p")
  if(!assertive::is_a_bool(useP)){
    useP <- as.logical(useP[[1]])
    warning("The option 'RNAmodR_use_p' is not a single logical value. ",
            "Please set 'RNAmodR_use_p' to TRUE or FALSE.",
            call. = FALSE)
  }
  useP
}

# returns the default colour palette
.get_color_palette <- function(){
  palette <- getOption("RNAmodR_palette")
  if(!assertive::is_a_string(palette)){
    palette <- RNAMODR_DEFAULT_PALETTE
    warning("The option 'RNAmodR_palette' is not a single string. ",
            "Please set 'RNAmodR_palette' to a valid palette identifier using ",
            "a single string.",
            call. = FALSE)
  }
  palette
}

