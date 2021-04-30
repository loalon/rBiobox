#' Plot Genes
#'
#' Given a matrix of gene expression data and a set of genes
#' obtain the plot for the individual genes of that set.
#'
#' @param data Gene expression matrix, samples in rows, genes in columns
#' @param genes List of genes to be represented
#' @param xVector Vector with the sample common characteristics (time, stage, ...)
#' @param xLabel Label of the x axis e.g. Time
#' @param yLabel Label of the y axis e.g. VST
#' @param title Optional plot title
#' @param grid Shows grid on the plot
#' @param legend Shows legend
#' @param legendTitle Adds a title to the legend
#' @param interactive Generates an interactive plot
#' @param scale center and scale the gene expression matrix
#'
#' @examples
#' plotGenes(toyData$expressionData, toyData$geneCluster[1:2], toyData$time, scale=F, xLabel="time", yLabel="VST", grid=F)
#'
#' @export
plotGenes <- function(data, genes, xVector, xLabel = "", yLabel="", title = "",
                      grid = F, legend = T, legendTitle = "Genes",
                      interactive = T, scale = T) {

  # Checks
  if(nrow(data) != length(xVector) )
    stop("xVector must contain the same number of items as rows in the expression data")

  if (length(genes) < 1 )
    stop("This function needs at least one gene")

  if (!all(genes %in% colnames(data)))
    stop("Gene expression data does not contain all query genes")

  if(!is.numeric(as.matrix(data)))
    stop("The expression data is not numeric")
    ## add option to different lines shapes

  # subset the data only with genes of interest
  localData <- data[, genes]

  # Scale if needed
  if(scale){
    localData <- scale(localData)
  }

  # join data with meta and melt
  localData <- as.data.frame(
    cbind(localData, xVector = xVector, sampleID = rownames(localData))
  )

  localData <- reshape2::melt(localData, measure.vars = genes)

  # prepare the minimal plot
  finalPlot <-
    ggplot2::ggplot(localData, ggplot2::aes(x = xVector, y = value,
                                            group = variable, col = variable)) +
    ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "ribbon",
                          fill = "lightgrey", alpha = 0.75) +
    ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "line", lwd = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y=ggplot2::element_blank())


  #TODO add orientation
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))

  #TODO check if time is factor and advise the user to manually stablish levels
  #if they are not numeric

  # If no grid
  if(!grid) {
    finalPlot <- finalPlot +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            )
  }

  # legend is TRUE
  finalPlot <- if(legend) {
    finalPlot +
      ggplot2::theme(legend.position = "right") +
      ggplot2::labs(color = legendTitle)
  } else {
    finalPlot +
      ggplot2::theme(legend.position = "none")
  }

  # Add general attributes
  finalPlot <- finalPlot +
    ggplot2::xlab(xLabel) +
    ggplot2::ylab(yLabel) +
    ggplot2::labs(title = title)

  # Return interactive or normal plot
  if(interactive){
    return(plotly::ggplotly(finalPlot))
  } else {
    return(finalPlot)
  }

}



