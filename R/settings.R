#' @include RNAmodR.R
NULL

# testthat

.test_test_TRUE <- function(x){TRUE}
.test_test_FALSE <- function(x){FALSE}

# tests ------------------------------------------------------------------------

.not_logical_operator <- function(x){
  .empty_character(x) | !(x %in% c("|","&"))
}

.not_single_numeric <- function(x){
  !is.numeric(x) | length(x) != 1 | is.na(x)
}
.not_numeric_between_0_100 <- function(x){
  .not_single_numeric(x) | x < 0 | x > 100
}
.not_numeric_between_0_1 <- function(x){
  .not_single_numeric(x) | x < 0 | x > 1
}
.not_numeric_bigger_zero <- function(x){
  .not_single_numeric(x) | x < 0
}

.not_single_integer <- function(x){
  !is.integer(x) | length(x) != 1 | is.na(x)
}
.not_integer_bigger_than_10 <- function(x){
  .not_single_integer(x) | x <= 10L
}
.not_integer_bigger_than_zero <- function(x){
  .not_single_integer(x) | x <= 0L
}
.not_integer_bigger_equal_than_zero <- function(x){
  .not_single_integer(x) | x < 0L
}
.not_integer_bigger_equal_than_one <- function(x){
  .not_single_integer(x) | x <= 1L
}

.is_not_GRanges_or_GRangesList <- function(x){
  !is(x,"GRanges") && !is(x,"GRangesList")
}

#' @importFrom grDevices col2rgb
.are_colours <- function(x) {
  vapply(x,
         function(z) {
           tryCatch(is.matrix(grDevices::col2rgb(z)),
                    error = function(e) FALSE)
         },
         logical(1))
}
.not_colours <- function(x){
  !is.character(x) | any(!.are_colours(x))
}

.empty_character <- function(x){
  if(!is.character(x) | is.na(x)){
    return(TRUE)
  }
  width(x) == 0L
}

# test import from assertive ---------------------------------------------------

.is_a_bool <- assertive::is_a_bool
.is_numeric_string <- assertive::is_numeric_string
.is_a_string <- assertive::is_a_string

# testing settings -------------------------------------------------------------

.get_name_in_parent_list <- function(...){
  xnames <- assertive::get_name_in_parent(list(...))
  xnames <- gsub("list\\(","",gsub("\\)","",xnames))
  xnames <- strsplit(xnames,", ")[[1]]
  xnames
}

.test_setting <- function(xname, settings, defaults, input){
  test <- settings$variable == xname
  FUN <- as.character(settings[test,"testFUN"])
  default <- defaults[[xname]]
  input <- input[[xname]]
  if(is.null(input)){
    return(default)
  }
  FUN <- get(FUN)
  if(FUN(input) == settings[test,"errorValue"]){
    stop(as.character(settings[test,"errorMessage"]), call. = FALSE)
  }
  input
}

.norm_settings <- function(input, settings, ...){
  if(!all(c("variable","testFUN","errorValue","errorMessage") %in% colnames(settings))){
    stop("Invalid columns in settings test definition.", call. = FALSE)
  }
  if(any(duplicated(settings$variable))){
    stop("Duplicated variable names in settings test definition.",
         call. = FALSE)
  }
  xnames <- .get_name_in_parent_list(...)
  defaults <- list(...)
  names(defaults) <- xnames
  f <- xnames %in% settings$variable
  if(!all(f)){
    stop("Test for variables '",
         paste(xnames[!f],collapse = "', '"),
         "' not found.", call. = FALSE)
  }
  args <- lapply(xnames, .test_setting, settings, defaults, input)
  names(args) <- xnames
  args
}
