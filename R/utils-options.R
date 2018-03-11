#' @include RNAmodR.R
NULL

# option retrieval -------------------------------------------------------------

# whether to use a p value for detection or not
.get_map_quality <- function(){
  adviseText <- "Please set 'RNAmodR_map_quality' to a numeric value of >= 0."
  mapQ <- getOption("RNAmodR_map_quality")
  if(!assertive::is_a_number(mapQ)){
    mapQ <- RNAMODR_DEFAULT_MAPQ
    warning("The option 'RNAmodR_map_quality' is not a single numeric value. ",
            adviseText,
            call. = FALSE)
  }
  if(!(mapQ >= 0)){
    mapQ <- RNAMODR_DEFAULT_MAPQ
    warning("The option 'RNAmodR_map_quality' is not 0 or a positive value ",
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

# returns the width used for RiboMethScore calculate
.get_ribometh_score_width <- function(){
  width <- getOption("RNAmodR_RiboMethScore_width")
  if(!assertive::is_a_number(width) || !assertive::is_positive(width)){
    width <- RNAMODR_DEFAULT_RMS_WIDTH
    warning("The option 'RNAmodR_RiboMethScore_width' is not a number. ",
            "Please set 'RNAmodR_RiboMethScore_width' to a positive number. ",
            call. = FALSE)
  }
  width
}
# returns the width used for RiboMethScore calculate
.get_ribometh_score_weights <- function(){
  weights <- getOption("RNAmodR_RiboMethScore_weights")
  if(length(weights) != ( 2 * .get_ribometh_score_width() + 1 ) ){
    stop("The option 'RNAmodR_RiboMethScore_weights' is not a compatible with ",
         "the width of 'RNAmodR_RiboMethScore_width'. ",
         "Please set 'RNAmodR_RiboMethScore_weights' to a named list of ",
         "weights. ",
         call. = FALSE)
  }
  if(weights["0"] != 0){
    stop("Error in weighting list.")
  }
  ## more checks needed XXX
  weights
}
