\name{plotFit}
\alias{plotFit}
\title{Network Fit Plot}
\description{
  This function plots the graph corresponding to the minimum scoring network.
}

\usage{
 plotFit(ternaryFit, type="interactive", ...)
}
\arguments{
  \item{ternaryFit}{a ternaryFit object}
  \item{type}{the type of plot to produce. "interactive" produces a plot
    that can be altered in the plotting window using the tkplot
    function from the igraph package. "static" produces a standard plot
    in any R graphics device.}
  \item{...}{additional parameters passed to the plotting function}
}
\value{
  A plot of the network corresponding to the minimum score (stored in the graphObjMin slot) is plotted.
}
\author{Matthew N. McCall and Anthony Almudevar}
\seealso{
  Almudevar A, McCall MN, McMurray H, Land H (2011). Fitting Boolean
  Networks from Steady State Perturbation Data, Statistical Applications
  in Genetics and Molecular Biology, 10(1): Article 47.
  }
\examples{
ssObj <- matrix(c(1,1,1,0,1,1,0,0,1),nrow=3)
pObj <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3)
tnfitObj <- tnetfit(ssObj, pObj)
plotFit(tnfitObj, type="static")
}
\keyword{manip}

