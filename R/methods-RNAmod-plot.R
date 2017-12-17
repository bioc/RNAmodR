#' @include RNAmod.R
NULL


#' @name visualizeModifications
#' 
#' @title visualizeModifications
#' 
#' @param .Object a RNAmod object. 
#' @param number a number defining an experiment. 
#' @param modifications name of modification to be used for analysis. 
#'
#' @return
#' @export
#'
#' @examples
setMethod(
  f = "visualizeModifications", 
  signature = signature(.Object = "RNAmod",
                        number = "numeric",
                        modifications = "character",
                        genes = "character"),
  definition = function(.Object,
                        number,
                        modifications,
                        genes){
    # Check input
    assertive::is_a_number(number)
    assertive::assert_all_are_non_empty_character(modifications)
    assertive::assert_all_are_non_empty_character(genes)
    
    folder <- paste0(getOutputFolder(.Object),"Visualizations/")
    
    # get intersection of requested and available genes
    genesAvail <- intersect(rownames(se), genes)
    if(length(genesAvail) != length(genes)){
      if(length(genesAvail) > 0){
        warning("No information available for the following genes: '",
                paste0(setdiff(genes,genesAvail), collapse = "', '"),"'",
                call. = FALSE)
      } else {
        stop("No information available for any genes given: '",
             paste0(genes, collapse = "', '"),"'",
             call. = FALSE)
      }
    }
    
    browser()
    positions <- SummarizedExperiment::rowData(se[rownames(se) %in% genesAvail,])$positions
    mods <- SummarizedExperiment::rowData(se[rownames(se) %in% genesAvail,])$mods
    gff <- .Object@.dataGFF
    fasta <- .Object@.dataFasta
    
    for(i in seq_along(genesAvail)){
      .plot_gene_with_modifications(positions,mods,gff,fasta)
    }
    
  }
)

.plot_gene_with_modifications <- function(positions,mods,gff,seq){
  
}

.save_gene_mod_plot <- function(){
  
}