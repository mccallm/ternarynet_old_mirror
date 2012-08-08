\name{plot-methods}
\alias{plot, ternaryFit-method}
\alias{plot, ternaryPost-method}
\title{Plotting functions for ternaryFit and ternaryPost objects}
\description{
  These functions plot visual representations of networks based on the
  type of object -- ternaryFit or ternaryPost.
}
\usage{
plot(ternaryFit, type="interactive", ...)
plot(ternaryPost, threshold=0.5, type="interactive", ...) 
}
\arguments{
  \item{ternaryFit}{a ternaryFit object}
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
  For a ternaryFit object, the network corresponding to the minimum
  score (stored in the graphObjMin slot) is plotted.
  For a ternaryPost object, a network consisting of all edges with a
  marginal posterior probability greater than the selected threshold. 
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
plot(tnfitObj)
}
\keyword{manip}

