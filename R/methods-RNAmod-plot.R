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
#' @import ggplot2
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
    
    se <- getSummarizedExperiment(.Object, number, modifications)
    
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
    
    positions <- SummarizedExperiment::rowData(se[rownames(se) %in% genesAvail,])$positions[[1]]
    mods <- SummarizedExperiment::rowData(se[rownames(se) %in% genesAvail,])$mods
    gff <- .Object@.dataGFF
    fasta <- .Object@.dataFasta
    
    modClasses <- .load_mod_classes(modifications)
    
    for(i in seq_along(genesAvail)){
      .plot_gene_with_modifications(genesAvail[[i]],positions,mods,gff,fasta,modClasses)
    }
    
  }
)

.plot_gene_with_modifications <- function(geneName,positions,mods,gff,seq,modClasses){
  gff_sub <- gff[(is.na(S4Vectors::mcols(gff)$ID) & 
                S4Vectors::mcols(gff)$Name == geneName) |
               (!is.na(S4Vectors::mcols(gff)$ID) & 
                  S4Vectors::mcols(gff)$ID == geneName),]
  seq <- Rsamtools::getSeq(fasta,gff_sub)[[1]]
  
  # create description layer plot
  layer <- .create_layer_data(gff,gff_sub)
  # layerPlot <- .get_gene_plot(geneName, layer)
  browser()
  
  # get data for different plot types
  posData <- .aggregate_pos_data(positions,modClasses)
  modData <- .aggregate_mod_data(mods,modClasses)
  
  dataPlots <- lapply(names(posData), function(type){
    pos <- posData[[type]]
    mods <- posData[[type]]
    
    plot <- .get_mod_plot(pos,mods)
    
  })
  
}


# returns a data.frame with the coordinates for plotting gene features
.create_layer_data <- function(gRanges,gRangeSelected){
  columns <- c("type","ID","Name","Parent")
  gRangeSelected <- gRangeSelected[,columns]
  gRanges <- gRanges[,columns]
  gRanges <- IRanges::subsetByOverlaps(gRanges, gRangeSelected)
  gRanges <- gRanges[as.character(S4Vectors::mcols(gRanges)$type) != "chromosome",]
  
  if( length(gRanges) == 0){
    return(NULL)
  }
  df <- S4Vectors::mcols(gRanges)[,c("type","ID","Name")]
  df$start <- BiocGenerics::start(gRanges)
  df$end <- BiocGenerics::end(gRanges)
  df$strand <- BiocGenerics::strand(gRanges)
  df <- as.data.frame(df)
  return(df)
}

# aggregates the position data for plotting of different plot types
.aggregate_pos_data <- function(positions, modClasses){
  plotTypes <- vapply(modClasses, getPlotType, character(1))
  plotTypes <- unique(plotTypes)
  posData <- lapply(plotTypes, function(type){
    if(is.null(positions[[type]])){
      stop("Positions data not available for type '",
           type,
           "'",
           .Call = FALSE)
    }
    positions[[type]]
  })
  names(posData) <- plotTypes
  return(posData)
}

# returns a plot showing all modifications on one type of position data
.get_read_plot <- function(pos,mods){
  requireNamespace("ggplot2", quietly = TRUE)
  
  plot <- ggplot(posData[[1]], aes_(x = ~pos, y = ~mean)) +
    geom_bar(stat = "identity")
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





#' @title .get_gene_plot
#' 
#' @description
#' returns a plot for the gene description
#'
#' @param geneName the name of the gene
#' @param dfLayer the data.frame containing the data for plotting the annotation
#' of gene, CDS and uORFS
#' @param xlim limits on the x axis
#' 
#' @return ggplot
#' 
#' @import ggplot2
.get_gene_plot <- function(geneName, 
                           dfLayer){
  requireNamespace("ggplot2", quietly = TRUE)
  ymin <- 0
  ymin2 <- 0
  
  arrowObj <- arrow(ends = "last", type = "closed")
  if( unique(dfLayer$strand) == "-" ){
    arrowObj <- arrow(ends = "first", type = "closed")
  }
  
  plot <- ggplot() +
    theme_minimal() +
    theme(legend.position = "none") +
    .get_read_line_colours() +
    scale_y_discrete(geneName) +
    labs(y = "none", x = "") + 
    xlim(xlim)
  # plot chromosome
  plot <- plot + geom_segment(aes_(x = xlim[1],
                                   y = ~ymin2, 
                                   xend = xlim[2], 
                                   yend = ~ymin2, 
                                   colour = "chromosome"), 
                              size = 1, 
                              arrow = arrowObj)
  # plot gene
  plot <- plot + geom_segment(aes_(x = dfLayer[dfLayer$type == "gene","start"], 
                                   y = ~ymin2,
                                   xend = dfLayer[dfLayer$type == "gene","end"], 
                                   yend = ~ymin2, 
                                   colour = "gene"), 
                              size = 3)
  # plot mRNA
  plot <- plot + geom_segment(aes_(x = dfLayer[dfLayer$type == "mRNA","start"], 
                                   y = ~ymin2, 
                                   xend = dfLayer[dfLayer$type == "mRNA","end"], 
                                   yend = ~ymin2, 
                                   colour = "mRNA"), 
                              size = 1.5)
  # plot CDS
  plot <- plot + geom_segment(aes_(x = dfLayer[dfLayer$type == "CDS","start"], 
                                   y = ~ymin2, xend = dfLayer[dfLayer$type == "CDS","end"], 
                                   yend = ~ymin2, 
                                   colour = "CDS"), 
                              size = 2.0)
  if( ("five_prime_UTR" %in% dfLayer$type) ){
    # plot CDS
    plot <- plot + geom_segment(aes_(x = dfLayer[dfLayer$type == "five_prime_UTR","start"], 
                                     y = ~ymin2, 
                                     xend = dfLayer[dfLayer$type == "five_prime_UTR","end"], 
                                     yend = ~ymin2, 
                                     colour = "five_prime_UTR"), 
                                size = 1.0)
  }
  if( ("three_prime_UTR" %in% dfLayer$type) ){
    # plot CDS
    plot <- plot + geom_segment(aes_(x = dfLayer[dfLayer$type == "three_prime_UTR","start"], 
                                     y = ~ymin2, 
                                     xend = dfLayer[dfLayer$type == "three_prime_UTR","end"], 
                                     yend = ~ymin2, 
                                     colour = "three_prime_UTR"), 
                                size = 1.0)
  }
  if( ("biological_region" %in% dfLayer$type) ){
    uORFdata <- dfLayer[dfLayer$type == "biological_region",]
    for(i in seq_len(nrow(uORFdata))){
      # plot CDS
      plot <- plot + geom_segment(aes_(x = uORFdata[i,"start"], 
                                       y = ~ymin2, 
                                       xend = uORFdata[i,"end"], 
                                       yend = ~ymin2, 
                                       colour = "uORF"), 
                                  size = 2.0)
    }
  }
  
  return(plot)
}





.save_gene_mod_plot <- function(){
  
}



