% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotEigengene.R
\name{plotEigengene}
\alias{plotEigengene}
\title{Plot Eigengene}
\usage{
plotEigengene(
  data,
  genes,
  xVector,
  yVector,
  xLabel = "",
  yLabel = "",
  inverse = F,
  title = "",
  grid = F,
  colors = c("#009444", "#BE9230"),
  legend = T,
  multiline = F,
  legendTitle = "",
  plotBothProfiles = F
)
}
\arguments{
\item{data}{Gene expression matrix, samples in rows, genes in columns}

\item{genes}{List of genes to be represented}

\item{xVector}{Vector with the samples common characteristics (time, stage, ...)}

\item{yVector}{Vector with the sample conditions (control, fertilised, ...)}

\item{xLabel}{Label of the x axis e.g. Time}

\item{yLabel}{Label of the y axis e.g. VST}

\item{inverse}{Plot the reverse profile}

\item{title}{Optional plot title}

\item{grid}{Shows grid on the plot}

\item{colors}{Color to use with the y-axis data}

\item{legend}{Shows legend}

\item{multiline}{Plot individual genes, instead of the eigengene}

\item{legendTitle}{Legend title}

\item{plotBothProfiles}{Plot main eigengene profile and the reverse one}
}
\description{
Given a matrix of gene expression data and a set of genes
obtain the plot for the eigengene of that set.
}
\examples{
plotEigengene(toyData$expressionData, toyData$geneCluster, toyData$time, toyData$conditions)

}
