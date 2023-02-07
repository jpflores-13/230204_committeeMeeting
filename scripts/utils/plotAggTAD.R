#' @param AggTAD Aggregated TAD object to plot
#' @param maxval integer representing the maximum value to plot (essentialy sets the top of the zrange)
#' @param cols color palette for plotting
#' @param title title of the plot

### Define a function that builds aggregate TADS
plotAggTAD <- function(AggTAD,maxval = 10000,cols = RColorBrewer::brewer.pal(6,"YlGnBu"),title="")
{
  # Convert to long format for ggplot
  AggTAD_long = setNames(melt(AggTAD), c('x', 'y', 'counts'))
  AggTAD_long$counts = AggTAD_long$counts
  
  ggplot(data=AggTAD_long,mapping=aes(x=x,y=y,fill=counts)) + 
    geom_tile() + 
    theme_void() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(aspect.ratio = 1) + 
    ggtitle(title) +
    scale_fill_gradientn(colours = cols,
                         na.value=cols[maxval],
                         limits=c(0,maxval),
                         oob = scales::squish) 
}