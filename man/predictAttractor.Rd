\name{predictAttractor}
\alias{predictAttractor}
\title{Predict the attractor(s) resulting from a given perturbation}
\description{
  This function computes the posterior probabilities of attractors reached for a
  given perturbation using the networks from a ternaryPost object.
}
\usage{
predictAttractor(tpost, perturbations, wildtype = TRUE, verbose = FALSE)
}
\arguments{
  \item{tpost}{a ternaryPost object}
  \item{perturbations}{a list with two elements: perturbed.genes and
    forced.states}
  \item{wildtype}{if TRUE, the wildtype attractors are summarized; if
    FALSE, the perturbed attractors are summarized.}
  \item{verbose}{if TRUE, periodic reports on progress are printed.}
}
\value{
  The function returns a list with two elements: \\
  post.prob: the posterior probability of each attractor \\
  attractor.summary: a single vector of steady states based on the
  resulting attractor
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
predictAttractor(tnpostObj, list(perturbed.genes=c(1,2),forced.states=c(1,1)))
}
\keyword{manip}

