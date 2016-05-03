\name{attractorSummary}
\alias{attractorSummary}
\title{Summarize Attractors}
\description{
  This function summarizes the posterior probability of possible attractors.
}
\usage{
attractorSummary(tpost, post.prob.limit = 0.01, wildtype = TRUE)
}
\arguments{
  \item{tpost}{a ternaryPost object}
  \item{post.prob.limit}{the minimum posterior probability for an
    attractor to be listed}
  \item{wildtype}{if TRUE, the wildtype attractors are summarized; if
    FALSE, the perturbed attractors are summarized.}
}
\value{
  The function returns a matrix of attractors and posterior
  probabilities for each perturbation.
}
\author{Matthew N. McCall and Anthony Almudevar}
\seealso{Almudevar A, McCall MN, McMurray H, Land H (2011). Fitting
  Boolean Networks from Steady State Perturbation Data, Statistical
  Applications in Genetics and Molecular Biology, 10(1): Article 47.}
\examples{
ssObj <- matrix(c(1,1,1,0,1,1,0,0,1),nrow=3)
pObj <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3)
rownames(ssObj) <- rownames(pObj) <- colnames(ssObj) <- colnames(pObj) <- c("Gene1","Gene2","Gene3")
tnfitObj <- tnetfit(ssObj, pObj)
tnpostObj <- tnetpost(tnfitObj, mdelta=10, msample=10)
attractorSummary(tnpostObj)
}
\keyword{manip}

