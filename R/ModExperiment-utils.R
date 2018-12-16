#' @include RNAmodR.R
NULL



.get_ModExperiment_objects <- function(...){
  args <- list(...)
  ModExperiments <- args[lapply(args,is,"ModExperiment")]
  ModExperiments
}
