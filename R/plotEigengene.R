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
plotEigengene <- function(data, genes, xVector, yVector, xLabel = "", yLabel="",
                          inverse = F, title = "", grid = F, colors = c("#009444", "#BE9230"),
                          legend=T, multiline = F, legendTitle="", plotBothProfiles=F) {

  if(nrow(data) != length(xVector) || nrow(data) != length(yVector)){
    stop("xVector and yVector must contain the same number of items as rows in the expression data")
  }

  if(plotBothProfiles) {
    if(length(unique(yVector)) > 1)
      stop("Plotting both expression profiles can only be done with a one condition dataset, e.g. Control.")
    if(multiline)
      stop("Can't plot both expression profiles in multiline mode.")
    if(length(genes) == 1 )
      stop("Can't plot both expression profiles with only one gene.")
    if(inverse)
      stop("Can't plot both expression profiles with the inverse parameter activated.")
  }

  if (length(genes) < 1 )
    stop("This function needs at least one gene")


  ## add option to different lines shapes

  ## maybe not scale, or scale optional
  expr <- NA

  #d <- data[,which(colnames(data) %in% gene)]
  d <- data[, genes]


  if (length(genes) == 1 ) { #1 gene, directly the data
    expr <- scale(d)

  } else if (length(genes > 1)) { #cluster, eigengene
    pca <- prcomp((d))
    pc1 <- pca$x[, 1]
    #print(percents <- round(summary(pca)$importance[2,]*100))

    if(sum(sign(cor(pca$x[,1,drop = FALSE], d))) < 0) {
      pc1 <- pc1 * -1
    }
    if(inverse && !plotBothProfiles)
      pc1 <- pc1 * -1


    expr <- pc1

  }

  myplot <- NA

  if(multiline && length(genes) > 1 ) {
    d <- as.data.frame(scale(d))
    d <- as.data.frame(cbind(d, time = time, sampleID=rownames(d)))
    d <- reshape2::melt(d, measure.vars = genes)

    myplot <- ggplot2::ggplot(d, aes(x = time, y = value, group = variable, col = variable)) +
      stat_summary(fun.data = mean_se, geom = "line", lwd = 1) +
      theme_bw()

  } else {
    myplot <- ggplot2::ggplot(data.frame(x = xVector, y = scale(expr), g = yVector),
                     aes(x = x, y = y, group = g)) +
      stat_summary(fun.data = mean_se, geom = "ribbon", fill = "lightgrey", alpha = 0.75) +
      stat_summary(fun.data = mean_se, geom = "line", aes(col = g), lwd = 2) + #      plot_output_list <- lapply(shiftedColors, function(color) {

      labs(color = legendTitle)

  }

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

    ratio <- round(length(negative)/(length(positive)+length(negative)) ,digits=2 )
    mData <- data.frame(PC1 = expr, PC1i = expr * -1)
    #mData$PC1i <- expr * -1

    mData <- reshape2::melt(mData, measure.vars = c("PC1", "PC1i"))
    mData$xVector <- xVector
    #mData$Treatment <- type

    myplot <- ggplot2::ggplot(mData, aes(x = xVector, y =scale(value), group = variable, col = variable, alpha=variable)) +
      #stat_summary(fun.data = mean_se, alpha = 0.25, geom = "ribbon", col = "grey90") +

      stat_summary(fun.data = mean_se, lwd = 2, geom = "line") +
      scale_alpha_manual(values=c(1,ratio))
    #scale_alpha_discrete(range = c(0.35, 0.9))
      #theme_bw() +
      #facet_wrap(~Treatment, nrow = 2, scales = "free_y") +
      #theme(text = element_text(size = 10)) +
      #scale_color_discrete(name = "")

  }


  #TODO add orientation
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))

  #TODO check if time is factor and advise the user to manually stablish levels if they are not numeric

  if(!grid) {
    myplot <- myplot +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())
  }

  if(multiline == F && is.vector(colors) && length(colors) >= length(unique(yVector))) {
    myplot <- myplot + scale_color_manual(values=colors)
  }

  if(legend) {
    myplot <- myplot +
    theme(legend.position = "none")
  }

  ##add general stuff
  myplot <- myplot +
    xlab(xLabel) +
    ylab("z-score") +
    #ggtitle(title)
    labs(title = title)

  return (myplot)

}



