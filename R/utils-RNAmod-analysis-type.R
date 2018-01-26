#' @include class-RNAmod-analysis-type.R
NULL

# sample plotting --------------------------------------------------------------

.plot_sample_data <- function(df,name){
  requireNamespace("ggplot2", quietly = TRUE)
  palette <- .get_color_palette()
  # Fix zero values for log display
  df2 <- df[df$y != 0,]
  # plot data and create grid
  plot <- ggplot(df, aes_(x = ~x, y = ~y, colour = ~group)) +
    scale_x_discrete(name = "position") + 
    scale_colour_brewer(palette = palette,
                        name = "position\ndata")
  # Log display
  plot2 <- plot + 
    geom_violin(data = df2,
                trim = FALSE) + 
    geom_jitter(data = df2,
                width = 0.1) +
    scale_y_log10(name = "counts per position", 
                  labels = scales::scientific)
  #Linear display
  plot <- plot + 
    geom_violin(trim = FALSE) + 
    geom_jitter(width = 0.1) +
    scale_y_continuous(name = "counts per position", 
                       labels = scales::scientific)
  grid <- gridExtra::grid.arrange(plot, plot2, nrow = 1, ncol = 2)
  # Resurrect calling RNAmod object to get the output folder
  object <- get(".Object", envir = .where(".Object"))
  # create file path 
  folder <- paste0(getOutputFolder(object),
                   "SamplePos/")
  if(!assertive::is_dir(folder)){
    dir.create(folder, recursive = TRUE)
  }
  fileName <- paste0(folder,
                     name,
                     ".pdf")
  # save plots as pdf
  ggsave(plot = grid,
         filename = fileName,
         width = 6,
         height = 4)
}

# common function handling positions -------------------------------------------

# converts the position on a transcript from start to end, to the correct
# genomic position
.convert_local_to_global_locations <- function(gff,
                                               loc){
  
  
  
  # browser()
  # Intron handling needs to be added
  strand <- unique(as.character(strand(gff)))
  if(strand == "-"){
    locations <- (end(gff) - loc) + 1
  } else {
    locations <- (start(gff) + loc) - 1
  }
  
  
  return(locations)
}

# converts global to local positions and modifies data accoringly
.convert_global_to_local_position <- function(gff,
                                              gr,
                                              data){
  # get a list of introns and the position which are void
  posToBeRemoved <- .get_intron_positions(gff,
                                          gr$ID)
  # get the genomic distance
  length <- width(gr) - length(unlist(posToBeRemoved))
  #
  #
  # GA specific handling
  # move position based on strand
  strand <- .get_unique_strand(gr)
  data <- data[.is_on_correct_strand(data,strand)]
  # interest in read's 5' position
  if(.is_on_minus_strand(gr)){
    positions <- BiocGenerics::end(data)
  } else {
    positions <- BiocGenerics::start(data)
  }
  # offset positions based on how many positions the read has passed from
  # transcription start
  positions <- unlist(lapply(positions, 
                             FUN = .move_postion,
                             posToBeRemoved, 
                             strand))
  # reset to relative positions to gene start
  if(.is_on_minus_strand(gr)){
    positions <- BiocGenerics::end(gr) - positions + 1
  } else {
    positions <- positions - BiocGenerics::start(gr) + 1
  }
  # discard reads out of boundaries
  positions <- positions[positions>0]
  positions <- positions[positions<=length]
  return(positions)
}

# offset positions based on how many positions the read has passed from
# transcription start
.move_postion <- function(positions,
                          posToBeRemoved,
                          strand){
  x <- unlist(posToBeRemoved)
  unlist(lapply(positions, function(position){
    if(.is_minus_strand(strand)){
      position <- position + length(x[x>position])
    } else {
      position <- position - length(x[x<position])
    }
    position
  }))
}

# aggregate the positions occupied by introns
.get_intron_positions <- function(gff, 
                                  ID){
  # get childs
  childs <- .subset_gff_for_unique_transcript(gff, 
                                              ID, 
                                              wo.childs = FALSE)
  lapply(childs[grepl("intron",childs$type),], 
         function(intron){
           Biostrings::start(intron):Biostrings::end(intron)
         })
}

.get_transcript_sequence <- function(gff,ID,seq){3
  # get gr
  gr <- .subset_gff_for_unique_transcript(gff, ID)
  # get a list of introns and the position which are void
  posToBeRemoved <- .get_intron_positions(gff,
                                          gr$ID)
  # convert DNAString into character vector
  seq <- .get_single_position_letters(as.character(seq))
  # set positions as names
  if(.is_on_minus_strand(gr)){
    names(seq) <- BiocGenerics::end(gr):BiocGenerics::start(gr)
  } else{
    names(seq) <- BiocGenerics::start(gr):BiocGenerics::end(gr)
  }
  # remove intron sequence
  seq <- seq[!(names(seq) %in% unlist(posToBeRemoved))]
  # set local position
  names(seq) <- 1:length(seq)
  return(seq)
}

# common function for converting data ------------------------------------------

# Calculate the arrest rate per position
.get_arrest_rate <- function(x){
  if(is.null(names(x))) stop("Unnamed position data.")
  y <- unlist(lapply(seq_along(x), function(i){
    a <- x[i]
    b <- x[(i+1)]
    # infinite in one direction
    if(is.na(b) || b == 0) return(-1)
    # max = 1 in one direction
    if(a <= b) return((1-a/b))
    # infinite in other direction
    if(is.na(a) || a == 0) return(1)
    # min = -1 in other direction
    return((b/a-1))
  }))
  setNames(y,names(x))
}

# Calculate the arrest difference per position
.get_arrest_diff <- function(i, meanData){
  x <- data[[i]]
  y <- unlist(lapply(seq_along(x), function(j){
    x[as.numeric(names(x)) == j]-meanData[[i]][as.numeric(names(x)) == j]
  }))
  setNames(y,names(x))
}

# calculates a rolling mean values
.get_rolling_mean <- function(x,
                              n=3){
  y <- setNames(filter(x,rep(1/n,n), sides=2),names(x))
  y[is.na(y)] <- x[is.na(y)]
  return(y)
}

# common function for subsetting data ------------------------------------------

.aggregate_location_data <- function(data, 
                                     location){
  unlist(lapply(data,function(dataPerReplicate){
    return(dataPerReplicate[as.numeric(names(dataPerReplicate)) == location])
  }))
}
.aggregate_not_location_data <- function(data,
                                         location){
  unlist(lapply(data,function(dataPerReplicate){
    return(dataPerReplicate[as.numeric(names(dataPerReplicate)) != location])
  }))
}
.aggregate_area_data <- function(data, 
                                 location, 
                                 width){
  unlist(lapply(data,function(dataPerReplicate){
    return(dataPerReplicate[as.numeric(names(dataPerReplicate)) < (location+width) &
                              as.numeric(names(dataPerReplicate)) > (location-width) &
                              as.numeric(names(dataPerReplicate)) != location])
  }))
}
