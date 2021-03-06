#' Plot Eigengene
#'
#' Given a matrix of gene expression data and a set of genes
#' obtain the plot for the eigengene of that set.
#'
#' @param data Gene expression matrix, samples in rows, genes in columns
#' @param genes List of genes to be represented
#' @param xVector Vector with the samples common characteristics (time, stage, ...)
#' @param yVector Vector with the sample conditions (control, fertilised, ...)
#' @param xLabel Label of the x axis e.g. Time
#' @param yLabel Label of the y axis e.g. VST
#' @param inverse Plot the reverse profile
#' @param title Optional plot title
#' @param grid Shows grid on the plot
#' @param colors Color to use with the y-axis data
#' @param legend Shows legend
#' @param multiline Plot individual genes, instead of the eigengene
#' @param legendTitle Legend title
#' @param plotBothProfiles Plot main eigengene profile and the reverse one
#'
#' @examples
#' plotEigengene(toyData$expressionData, toyData$geneCluster, toyData$time, toyData$conditions)
#'
#' @export
plotEigengene <- function(data, genes, xVector, yVector, xLabel = "", yLabel="z-score",
                          inverse = F, title = "", grid = F, colors = c("#009444", "#BE9230"),
                          legend=T, legendTitle="", plotBothProfiles=F) {

  if(nrow(data) != length(xVector) || nrow(data) != length(yVector)){
    stop("xVector and yVector must contain the same number of items as rows in the expression data")
  }

  if(plotBothProfiles) {
    if(length(unique(yVector)) > 1)
      stop("Plotting both expression profiles can only be done with a one condition dataset, e.g. Control.")
    if(length(genes) == 1 )
      stop("Can't plot both expression profiles with only one gene.")
    if(inverse)
      stop("Can't plot both expression profiles with the inverse parameter activated.")
  }

  if (length(genes) < 1 )
    stop("This function needs at least one gene")


  ## add option to different lines shapes

  localData <- data[, genes]

  expr <- if (length(genes) == 1 ) { #1 gene, directly the data
    scale(localData)

  } else if (length(genes > 1)) { #cluster, eigengene
    pca <- prcomp(localData)
    pc1 <- pca$x[, 1] # extract PC1

    # get the main profile
    if(sum(sign(cor(pca$x[,1,drop = FALSE], localData))) < 0) {
      pc1 <- pc1 * -1
    }

    # reverse profile
    if(inverse && !plotBothProfiles)
      pc1 <- pc1 * -1

    pc1

  }


  myplot <- ggplot2::ggplot(data.frame(x = xVector, y = scale(expr), g = yVector),
                            ggplot2::aes(x = x, y = y, group = g)) +
    ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "ribbon", fill = "lightgrey", alpha = 0.75) +
    ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "line", ggplot2::aes(col = g), lwd = 2) +
    ggplot2::labs(color = legendTitle)

  # print both profiles, related to their %
  if(plotBothProfiles) {

    pcaPos<- rownames(pca$rotation[which (pca$rotation[,'PC1'] > 0),])
    pcaNeg<- rownames(pca$rotation[which (pca$rotation[,'PC1'] < 0),])

    positive <- NA
    negative <- NA
    if (length(pcaPos) > length(pcaNeg)) {
      positive <- colnames(data[,which(colnames(data) %in% pcaPos)])
      negative <- colnames(data[,which(colnames(data) %in% pcaNeg)])
    } else {
      positive <- colnames(data[,which(colnames(data) %in% pcaNeg)])
      negative <- colnames(data[,which(colnames(data) %in% pcaPos)])
    }

    ratio <- round(length(negative)/(length(positive)+length(negative)), digits=2)
    mData <- data.frame(PC1 = expr, PC1i = expr * -1)

    mData <- reshape2::melt(mData, measure.vars = c("PC1", "PC1i"))
    mData$xVector <- xVector

    myplot <- ggplot2::ggplot(mData, ggplot2::aes(x = xVector, y =scale(value), group = variable, col = variable, alpha=variable)) +
      #stat_summary(fun.data = mean_se, alpha = 0.25, geom = "ribbon", col = "grey90") +

      ggplot2::stat_summary(fun.data = ggplot2::mean_se, lwd = 2, geom = "line") +
      ggplot2::scale_alpha_manual(values=c(1,ratio))
    #scale_alpha_discrete(range = c(0.35, 0.9))
      #theme_bw() +
      #facet_wrap(~Treatment, nrow = 2, scales = "free_y") +
      #theme(text = element_text(size = 10)) +
      #scale_color_discrete(name = "")

  }


  #TODO add orientation
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))

  #TODO check if time is factor and advise the user to manually establish levels if they are not numeric

  if(!grid) {
    myplot <- myplot +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank())
  }

  if(is.vector(colors) && length(colors) >= length(unique(yVector))) {
    myplot <- myplot + ggplot2::scale_color_manual(values=colors)
  }

  if(!legend) {
    myplot <- myplot +
      ggplot2::theme(legend.position = "none")
  }

  ##add general stuff
  myplot <- myplot +
    ggplot2::xlab(xLabel) +
    ggplot2::ylab(yLabel) +
    #ggtitle(title)
    ggplot2::labs(title = title)

  return (myplot)

}



