\name{plotPost}
\alias{plotPost}
\title{Network Posterior Plot}
\description{
  This function plots the graph consisting of all edges with a marginal posterior probability greater than the selected threshold. 
}

\usage{
plotPost(ternaryPost, threshold=0.5, type="interactive", ...) 
}
\arguments{
  \item{ternaryPost}{a ternaryPost object}
  \item{type}{the type of plot to produce. "interactive" produces a plot
    that can be altered in the plotting window using the tkplot
    function from the igraph package. "static" produces a standard plot
    in any R graphics device.}
  \item{threshold}{the marginal posterior probability required for an
    edge to be included in the plot.}
  \item{...}{additional parameters passed to the plotting function}
}
\value{
  A plot of the network consisting of all edges with a marginal posterior probability greater than the selected threshold. 
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
tnpostObj <- tnetpost(tnfitObj, mdelta=10, msample=10)
plotPost(tnpostObj, type="static")
}
\keyword{manip}

