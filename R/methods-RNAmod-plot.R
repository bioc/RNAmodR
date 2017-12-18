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
    
    positions <- SummarizedExperiment::rowData(se[rownames(se) %in% genesAvail,])$positions
    mods <- SummarizedExperiment::rowData(se[rownames(se) %in% genesAvail,])$mods
    gff <- .Object@.dataGFF
    fasta <- .Object@.dataFasta
    
    for(i in seq_along(genesAvail)){
      .plot_gene_with_modifications(genesAvail[[i]],positions,mods,gff,fasta)
    }
    
  }
)

.plot_gene_with_modifications <- function(geneName,positions,mods,gff,seq){
  
  
  gff <- gff[(is.na(S4Vectors::mcols(gff)$ID) & 
                S4Vectors::mcols(gff)$Name == geneName) |
               (!is.na(S4Vectors::mcols(gff)$ID) & 
                  S4Vectors::mcols(gff)$ID == geneName),]
  seq <- getSeq(fasta,gff)[[1]]
  browser()
  
  
  
  
  
}






.get_read_plot <- function(df, 
                           xlim, 
                           ylim, 
                           readType, 
                           title,
                           dfLayer){
  requireNamespace("ggplot2", quietly = TRUE)
  
  df <- as.data.frame(df)
  df$type <- rep(readType, nrow(df))
  
  plot <- ggplot()
  
  # plot CDS
  plot <- plot + geom_rect(aes_(xmin = dfLayer[dfLayer$type == "CDS","start"], 
                                ymin = 0, 
                                xmax = dfLayer[dfLayer$type == "CDS","end"], 
                                ymax = ~ylim, colour = "CDS", fill = "CDS"), 
                           alpha = 0.25)
  if( ("five_prime_UTR" %in% dfLayer$type) ){
    # plot CDS
    plot <- plot + geom_rect(aes_(xmin = dfLayer[dfLayer$type == "five_prime_UTR","start"], 
                                  ymin = 0, 
                                  xmax = dfLayer[dfLayer$type == "five_prime_UTR","end"], 
                                  ymax = ~ylim, colour = "five_prime_UTR", 
                                  fill = "five_prime_UTR"), 
                             alpha = 0.25)
  }
  if( ("three_prime_UTR" %in% dfLayer$type) ){
    # plot CDS
    plot <- plot + geom_rect(aes_(xmin = dfLayer[dfLayer$type == "three_prime_UTR","start"], 
                                  ymin = 0, 
                                  xmax = dfLayer[dfLayer$type == "three_prime_UTR","end"], 
                                  ymax = ~ylim, 
                                  colour = "three_prime_UTR",
                                  fill = "three_prime_UTR"), 
                             alpha = 0.25)
  }
  if( ("biological_region" %in% dfLayer$type) ){
    uORFdata <- dfLayer[dfLayer$type == "biological_region",]
    for(i in seq_len(nrow(uORFdata))){
      # plot CDS
      plot <- plot + geom_rect(aes_(xmin = uORFdata[i,"start"], 
                                    ymin = 0, 
                                    xmax = uORFdata[i,"end"], 
                                    ymax = ~ylim, 
                                    colour = "uORF", 
                                    fill = "uORF"), 
                               alpha = 0.25)
      
    }
  }
  if( ("intron" %in% dfLayer$type) ){
    intron_data <- dfLayer[dfLayer$type == "intron",]
    for(i in seq_len(nrow(intron_data))){
      # plot intron
      plot <- plot + geom_rect(aes_(xmin = intron_data[i,"start"], 
                                    ymin = 0, 
                                    xmax = intron_data[i,"end"], 
                                    ymax = ~ylim, 
                                    colour = "intron", 
                                    fill = "intron"), 
                               alpha = 0.25)
    }
  }
  
  
  plot <- plot +
    geom_col(mapping = aes_(~pos, ~mean, colour = ~type, fill = ~type), 
             data = df,
             size = 0,
             width = 1) +
    ylim(0,
         ylim) +
    xlim(xlim) +
    theme_minimal() +
    theme(legend.position = "none", 
          axis.line = element_blank()) +
    .get_read_line_colours() +
    .get_read_fill_colours() +
    labs(title = title, 
         y = paste0(unique(df$type), " [rpM]"), 
         x = "Genomic position")
  return(plot)
}

.get_read_line_colours <- function(){
  scale_colour_manual(values = c("Ribo" = "darkblue", 
                                 "RNA" = "darkgreen",  
                                 "gene" = "grey", 
                                 "mRNA" = "darkgrey", 
                                 "chromosome" = "black", 
                                 "CDS" = "orange",
                                 "intron" = "tan", 
                                 "uORF" = "red",
                                 "five_prime_UTR" = "violetred4",
                                 "three_prime_UTR" = "violetred"))
}
.get_read_fill_colours <- function(){
  scale_fill_manual(values = c("Ribo" = "darkblue", 
                               "RNA" = "darkgreen",  
                               "gene" = "grey", 
                               "mRNA" = "darkgrey", 
                               "chromosome" = "black", 
                               "CDS" = "orange", 
                               "intron" = "tan", 
                               "uORF" = "red",
                               "five_prime_UTR" = "violetred4",
                               "three_prime_UTR" = "violetred"))
}











.save_gene_mod_plot <- function(){
  
}



