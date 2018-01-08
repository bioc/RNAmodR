#' @include RNAmod.R
NULL

RNAMOD_PLOT_POS_WIDTH <- 0.5
RNAMOD_PLOT_LAYER_HEIGHT <- 50
RNAMOD_PLOT_DATA_HEIGHT <- 100

#' @name getModPlot
#' 
#' @aliases getModPlot saveModPlots
#' 
#' @title visualization of modifications
#' 
#' @param .Object a RNAmod object. 
#' @param number a number defining an experiment. 
#' @param modifications name of modification to be used for analysis. 
#' @param gene a single gene name
#'
#' @return a grid of layer plot and any mod plot
#' @export
#' 
#' @import ggplot2
#'
#' @examples
#' \donttest{
#' getModPlot(mod,1,"m7G","RDN18-1")
#' saveModPlot(mod,1,"m7G",c("RDN18-1","tS\(CGA\)P"))
#' }
setMethod(
  f = "getModPlot", 
  signature = signature(.Object = "RNAmod",
                        number = "numeric",
                        modifications = "character",
                        gene = "character"),
  definition = function(.Object,
                        number,
                        modifications,
                        gene,
                        focus){
    # Check input
    assertive::is_a_number(number)
    assertive::assert_all_are_non_empty_character(modifications)
    assertive::assert_all_are_non_empty_character(genes)
    assertive::assert_is_a_bool(focus)
    
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
    genesAvail <- unlist(genesAvail)
    assertive::assert_is_a_non_empty_string(genesAvail)
    
    positions <- SummarizedExperiment::rowData(se[rownames(se) %in% genesAvail,])$positions[[1]]
    mods <- SummarizedExperiment::rowData(se[rownames(se) %in% genesAvail,])$mods
    gff <- .Object@.dataGFF
    fasta <- .Object@.dataFasta
    
    modClasses <- .load_mod_classes(modifications)
    if(!focus){
      plot <- .plot_gene_with_modifications(genesAvail[[i]],
                                            positions,
                                            mods,
                                            gff,
                                            fasta,
                                            modClasses)
    } else {
      plot <- .plot_per_modifications(genesAvail[[i]],
                                      positions,
                                      mods,
                                      gff,
                                      fasta,
                                      modClasses)
    }
    return(plot$plot)
  }
)


#' @rdname getModPlot
#'
#' @return TRUE if successful
#' @export
#' 
#' @import ggplot2
setMethod(
  f = "saveModPlot", 
  signature = signature(.Object = "RNAmod",
                        number = "numeric",
                        modifications = "character",
                        genes = "character"),
  definition = function(.Object,
                        number,
                        modifications,
                        genes,
                        focus){
    # Check input
    assertive::is_a_number(number)
    assertive::assert_all_are_non_empty_character(modifications)
    assertive::assert_all_are_non_empty_character(genes)
    assertive::assert_is_a_bool(focus)
    
    # create folder
    experiment <- getExperimentData(.Object,
                                    number)
    folder <- paste0(getOutputFolder(.Object),
                     "ModPlots/",
                     unique(experiment$SampleName),
                     "/")
    if(!assertive::is_dir(folder)){
      dir.create(folder, recursive = TRUE)
    }
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
    mods <- SummarizedExperiment::rowData(se[rownames(se) %in% genesAvail,])$mods[[1]]
    gff <- .Object@.dataGFF
    fasta <- .Object@.dataFasta
    
    modClasses <- .load_mod_classes(modifications)
    
    if(!focus){
      for(i in seq_along(genesAvail)){
        plot <- .plot_gene_with_modifications(genesAvail[[i]],
                                              positions,
                                              mods,
                                              gff,
                                              fasta,
                                              modClasses)
        names <- paste0(paste(modifications, collapse = "-"),
                        "_",
                        genesAvail[[i]])
        .save_gene_mod_plot(folder,
                            names,
                            plot$plot,
                            plot$width,
                            plot$height)
      }
    } else {
      for(i in seq_along(genesAvail)){
        plot <- .plot_per_modifications(genesAvail[[i]],
                                        positions,
                                        mods,
                                        gff,
                                        fasta,
                                        modClasses)
        names <- paste0(genesAvail[[i]],
                        "_",
                        plot$names)
        .save_gene_mod_plot(folder,
                            names,
                            plot$plot,
                            plot$width,
                            plot$height)
      }
    }
    message("done.")
    return(invisible(TRUE))
  }
)

# plotting a complete gene -----------------------------------------------------

.plot_gene_with_modifications <- function(geneName,
                                          positions,
                                          mods,
                                          gff,
                                          seq,
                                          modClasses){
  gff_sub <- gff[(is.na(S4Vectors::mcols(gff)$ID) & 
                    S4Vectors::mcols(gff)$Name == geneName) |
                   (!is.na(S4Vectors::mcols(gff)$ID) & 
                      S4Vectors::mcols(gff)$ID == geneName),]
  seq <- Rsamtools::getSeq(fasta,gff_sub)[[1]]
  
  # create description layer plot
  layer <- .create_layer_data(gff,gff_sub)
  
  # get data for different plot types
  posData <- .aggregate_pos_data(gff_sub,positions,modClasses,lim)
  modData <- .aggregate_mod_data(mods,modClasses,lim)
  
  # to inhibit plotting by grid.arrange
  grDevices::pdf(file = NULL)
  
  dataPlots <- lapply(names(posData), function(type){
    pos <- posData[[type]]
    mods <- modData[[type]]
    
    plot <- .get_mod_plot(pos,mods)
    return(plot)
  })
  # layerPlot <- .get_gene_plot(geneName, layer)
  layerPlot <- list()
  
  nrow <- length(layerPlot) + length(dataPlots)
  width <- 40 + nrow(posData[[1]]) * RNAMOD_PLOT_POS_WIDTH
  height <- length(layerPlot) * RNAMOD_PLOT_LAYER_HEIGHT + 
    length(dataPlots) * RNAMOD_PLOT_DATA_HEIGHT
  
  # grid <- do.call(gridExtra::grid.arrange, append(dataPlots,
  #                                                 layerPlot,
  #                                                 list(nrow = nrow)))
  grid <- do.call(gridExtra::grid.arrange, append(dataPlots,
                                                  list(nrow = nrow)))
  
  # to inhibit plotting by grid.arrange
  grDevices::dev.off()
  
  return(list(plot = grid,
              width = width,
              height = height))
}

# plot per modification --------------------------------------------------------

.plot_per_modifications <- function(geneName,
                                    positions,
                                    mods,
                                    gff,
                                    seq,
                                    modClasses){
  gff_sub <- gff[(is.na(S4Vectors::mcols(gff)$ID) & 
                    S4Vectors::mcols(gff)$Name == geneName) |
                   (!is.na(S4Vectors::mcols(gff)$ID) & 
                      S4Vectors::mcols(gff)$ID == geneName),]
  seq <- Rsamtools::getSeq(fasta,gff_sub)[[1]]
  
  # create description layer plot
  layer <- .create_layer_data(gff,gff_sub)
  
  # get data for different plot types
  posData <- .aggregate_pos_data(gff_sub,positions,modClasses,lim)
  modData <- .aggregate_mod_data(mods,modClasses,lim)
  modData <- lapply(names(modData), function(x){
    y <- modData[[x]]
    y$plotType <- x
    y
  })
  modData <- do.call(rbind,modData)
  browser()
  
  
  # to inhibit plotting by grid.arrange
  grDevices::pdf(file = NULL)
  
  plots <- lapply(1:nrow(modData), function(i){
    mod <- modData[i,]
    pos <- posData[[mod$plotType]]
    
    
    
    # layerPlot <- .get_gene_plot(geneName, layer)
    layerPlot <- list()
    plot <- .get_mod_plot(pos,mod)
    
    nrow <- length(layerPlot) + 1
    width <- 40 + nrow(pos[[1]]) * RNAMOD_PLOT_POS_WIDTH
    height <- RNAMOD_PLOT_LAYER_HEIGHT + 
      RNAMOD_PLOT_DATA_HEIGHT
    
    # grid <- gridExtra::grid.arrange(plot, layerPlot, nrow = nrow)
    grid <- gridExtra::grid.arrange(plot, nrow = nrow)
    return(list(plot = grid,
                width = width,
                height = height))
  })
  plot <- lapply(plots, "[[","plot")
  width <- lapply(plots, "[[","width")
  height <- lapply(plots, "[[","height")
  
  # to inhibit plotting by grid.arrange
  grDevices::dev.off()
  
  return(list(plot = plot,
              width = width,
              height = height))
}



# data aggregation function ----------------------------------------------------

# aggregates the position data for plotting of different plot types
.aggregate_pos_data <- function(gff,data, modClasses,lim){
  plotTypes <- vapply(modClasses, getPlotType, character(1))
  plotTypes <- unique(plotTypes)
  data <- lapply(plotTypes, function(type, data){
    if(is.null(data[[type]])){
      stop("Positions data not available for type '",
           type,
           "'",
           .Call = FALSE)
    }
    res <- data[[type]]
    if( as.character(BiocGenerics::strand(gff)) == "-"){
      res <- res[order(res$pos, decreasing = TRUE),]
    }
    res$localPos <- 1:nrow(res)
    rownames(res) <- res$localPos
    return(res)
  }, data)
  names(data) <- plotTypes
  return(data)
}

# aggregates the modification data
.aggregate_mod_data <- function(data, modClasses,lim){
  modTypes <- vapply(modClasses, getModType, character(1))
  modTypes <- unique(modTypes)
  data <- lapply(modTypes, function(type, data){
    if(nrow(data[data$RNAmod_type == type,]) == 0){
      return(NULL)
    }
    data[data$RNAmod_type == type,]
  }, data)
  names(data) <- vapply(modClasses, getPlotType, character(1))
  return(data)
}


# modification visualization ---------------------------------------------------

# returns a plot showing all modifications on one type of position data
.get_mod_plot <- function(pos,mods){
  requireNamespace("ggplot2", quietly = TRUE)
  
  break_FUN <- function(lim){
    if(length(lim[is.na(lim)])>0) return(lim)
    x <- 0
    y <- floor(abs(max(lim)/100))*100
    diff <- y - x
    n <- floor(diff / 100)
    as.numeric(c(x,unlist(lapply(1:n,function(z){100*z})),round(max(lim))))
  }
  
  plot <- ggplot(pos, aes_(x = ~localPos, y = ~mean)) +
    scale_x_continuous(name = "position of transcript [nt]") +
    scale_y_continuous(name = "mean(number of read ends)",
                       labels = scales::scientific,
                       limits = c(NA,max(pos$mean)*1.25)) +
    theme_minimal()
  
  # prepare mod data
  mods$localStart <- pos[pos$pos == mods$start,"localPos"]
  mods$localEnd <- pos[pos$pos == mods$end,"localPos"]
  mods$vStart <- max(pos[pos$localPos < mods$localStart+10 &
                       pos$localPos > mods$localStart-10,"mean"])*1.01
  modsPositions <- mods[mods$start == mods$end,]
  modsArea<- mods[mods$start != mods$end,]
  
  
  if( nrow(modsArea) > 0){
  }
  plot <- geom_bar(stat = "identity")
  if( nrow(modsPositions) > 0){
    plot <- plot + geom_segment(data = as.data.frame(modsPositions),
                                mapping = aes_(x = ~localStart,
                                               y = ~vStart,
                                               xend = ~localStart,
                                               yend = Inf,
                                               colour = ~RNAmod_type)) +
      scale_colour_brewer(name = "modification\ntype",
                          palette = "Set1") 
    
    p_text <- unlist(lapply( modsPositions$RNAmod_p.value, function(p){
      if( p < 0.0001){
        return(paste0("< ",0.0001))
      } else {
        return(paste0(": ",round(p, digits = 4)))
      }
    }))
    
    label <- paste0(modsPositions$RNAmod_type,
                    "\n",
                    rownames(modsPositions),
                    "\n \u03C3: ",
                    modsPositions$RNAmod_signal,
                    " (p ",
                    p_text,
                    ")")
    plot <- plot + annotate("label",
                            x = modsPositions$localStart,
                            y = Inf,
                            label = label,
                            vjust = "inward",
                            size = 3)
  }
  
  return(plot)
}


# gene feature visualization ---------------------------------------------------

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

# returns a plot showing all the gene annotation layers
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







# plot saving ------------------------------------------------------------------

.save_gene_mod_plot <- function(folder,names,plot,width,height){
  # check that for all plot all values are available
  x <- c(length(names),
         length(plot),
         length(width),
         length(height))
  if(!all( abs(x - mean(x)) < .Machine$double.eps ^ 0.5 ))
    stop("Not the same number of names, plots, height and width values ",
         "provided.",
         .call = TRUE)
  assertive::assert_all_are_non_empty_character(names)
  assertive::assert_all_are_whole_numbers(width)
  assertive::assert_all_are_whole_numbers(height)
  
  fileNames <- paste0(folder,
                     "RNAmod_",
                     names,
                     ".pdf")
  for(i in 1:seq_along(length(fileNames))){
    
    
    ggplot2::ggsave(plot[[i]],
                    filename = fileNames[[i]],
                    units = "mm",
                    width = width[[i]],
                    height = height[[i]])
    message("Saving mod plot for gene '",names[[i]],"'")
  }
}



