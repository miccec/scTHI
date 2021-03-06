% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TME_plot.R
\name{TME_plot}
\alias{TME_plot}
\title{TME_plot}
\usage{
TME_plot(tsneData, Class, cexPoint = 0.8)
}
\arguments{
\item{tsneData}{X and y coordinates of points in the plot.}

\item{Class}{Object returned by TME_classification function.}

\item{cexPoint}{Set the point size.}
}
\value{
None
}
\description{
Generates a plot on the t-SNE coordinates, labeling cells by TME
classification.
}
\examples{
library(scTHI.data)
data(scExample)
result <-  scTHI_score(scExample,
           cellCusterA = colnames(scExample)[1:30],
           cellCusterB = colnames(scExample)[31:100],
           cellCusterAName = "ClusterA",
           cellCusterBName = "ClusterB", filterCutoff = 0,
           pvalueCutoff = 1, nPermu = 100, ncore = 8)
result <- scTHI_runTsne(result)
Class <- TME_classification(scExample)
TME_plot(tsneData = result$tsneData, Class)
}
