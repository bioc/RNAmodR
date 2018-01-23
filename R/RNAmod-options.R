
# option retrieval -------------------------------------------------------------

# whether to use a p value for detection or not
.get_map_quality <- function(){
  adviseText <- "Please set 'RNAmod_map_quality' to a numeric value of >= 0."
  mapQ <- getOption("RNAmod_map_quality")
  if(!assertive::is_a_number(mapQ)){
    mapQ <- RNAMOD_DEFAULT_MAPQ
    warning("The option 'RNAmod_map_quality' is not a single numeric value. ",
            adviseText,
            call. = FALSE)
  }
  if(!(mapQ >= 0)){
    mapQ <- RNAMOD_DEFAULT_MAPQ
    warning("The option 'RNAmod_map_quality' is not 0 or a positive value ",
            adviseText,
            call. = FALSE)
  }
  mapQ
}

# whether to use a p value for detection or not
.get_use_p <- function(){
  useP <- getOption("RNAmod_use_p")
  if(!assertive::is_a_bool(useP)){
    useP <- as.logical(useP[[1]])
    warning("The option 'RNAmod_use_p' is not a single logical value. ",
            "Please set 'RNAmod_use_p' to TRUE or FALSE.",
            call. = FALSE)
  }
  useP
}


.get_color_palette <- function(){
  palette <- getOption("RNAmod_palette")
  if(!assertive::is_a_string(palette)){
    palette <- RNAMOD_DEFAULT_PALETTE
    warning("The option 'RNAmod_palette' is not a single string. ",
            "Please set 'RNAmod_palette' to a valid palette identifier using ",
            "a single string.",
            call. = FALSE)
  }
  palette
}


.get_transcript_max_iteration <- function(){
  iterations <- getOption("RNAmod_transcript_max_iteration")
  if(!assertive::is_a_number(iterations)){
    iterations <- RNAMOD_DEFAULT_TRANSCRIPT_MAX_ITERATIONS
    warning("The option 'RNAmod_transcript_max_iteration' is not a single ",
            "number. Please set 'RNAmod_palette' to a valid palette ",
            "identifier using a single string.",
            call. = FALSE)
  }
  iterations
}
