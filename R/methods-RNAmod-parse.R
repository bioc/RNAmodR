#' @include RNAmod.R
NULL


setMethod(
  f = "parse", 
  signature = signature(.Object = "RNAmod",
                        number = "numeric"),
  definition = function(.Object,
                        number){
    assertive::is_a_number(number)
    
    experiment <- getExperimentData(.Object,number)
    
    FUN <- function(file,
                    gff,
                    fasta){
      
      
    }
    
    data <- bplapply(experiment$BamFile,
                     .Object@.dataGFF,
                     .Object@.dataFasta)
    
    
  }
)