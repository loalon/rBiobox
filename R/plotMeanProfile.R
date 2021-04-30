plotMeanProfile <- function(data, genes, condition, time, timeUnits = "Time",
         inverse = F, title = "", noGrid = T, colors = c("#009444", "#BE9230"),
         noLegend=T, multiline = F, legendTitle="",
         plotBothProfiles=F) {
  require(ggplot2)
  require(RColorBrewer)
  require(reshape2)

  if (length(genes) < 1 )
    stop("This function needs at least one gene")


  ## add option to different lines shapes

  ## maybe not scale, or scale optional
  expr <- NA

  #d <- data[,which(colnames(data) %in% gene)]
  d <- data[, genes]

  expr <- scale(d)


  myplot <- NA

    d <- as.data.frame(scale(d))
    d <- as.data.frame(cbind(d, time = time, condition = condition, sampleID=rownames(d)))
    d <- melt(d, measure.vars = genes)

    # myplot <- ggplot(d, aes(x = time, y = value, group = variable, col = variable)) +
    #   stat_summary(fun.data = mean_se, geom = "line", lwd = 1) +
    #   theme_bw()

    print(d)

    myplot <- ggplot(d,
                     aes(x = time, y = value, group = condition)) +
      stat_summary(fun.data = mean_se, geom = "ribbon", fill = "lightgrey", alpha = 0.75) +
      stat_summary(fun.data = mean_se, geom = "line", aes(col = condition), lwd = 2) + #      plot_output_list <- lapply(shiftedColors, function(color) {

      labs(color = legendTitle)



  #TODO add orientation
  #theme(axis.text.x = element_text(angle = 90, hjust = 1))

  #TODO check if time is factor and advise the user to manually stablish levels if they are not numeric

  if(noGrid) {
    myplot <- myplot +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())
  }

    myplot <- myplot + scale_color_manual(values=colors)

  if(noLegend) {
    myplot <- myplot +
      theme(legend.position = "none")
  }



  ##add general stuff
  myplot <- myplot +
    xlab(timeUnits) +
    ylab("z-score") +
    #ggtitle(title)
    labs(title = title)
  print("here")
  return (myplot)

}

# plotMeanProfile(rawData,
#                 sweetGenes,
#                 substr(rownames(rawData), 5, 5),
#                 as.numeric(substr(rownames(rawData), 2, 3)))
#
# plotMeanProfile(newTable,
#                 sweetGenes,
#                 substr(rownames(newTable), 4, 4),
#                 as.numeric(substr(rownames(newTable), 2, 3)))
