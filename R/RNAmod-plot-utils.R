#' @include methods-RNAmod-plot.R
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

