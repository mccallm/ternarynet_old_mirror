\name{ternaryFit-class}
\docType{class}
\alias{class:ternaryFit}
\alias{ternaryFit-class}
\alias{ternaryFit}
\alias{ternaryFit-methods}
\alias{initialize,ternaryFit-method}
\alias{dim,ternaryFit-method}
\alias{show,ternaryFit-method}
\alias{perturbationObj}
\alias{perturbationObj,ternaryFit-method}
\alias{steadyStateObj}
\alias{steadyStateObj,ternaryFit-method}
\alias{degreeObjMin}
\alias{degreeObjMin,ternaryFit-method}
\alias{graphObjMin}
\alias{graphObjMin,ternaryFit-method}
\alias{tableObjMin}
\alias{tableObjMin,ternaryFit-method}
\alias{newScore}
\alias{newScore,ternaryFit-method}
\alias{minScore}
\alias{minScore,ternaryFit-method}
\alias{finalTemperature}
\alias{finalTemperature,ternaryFit-method}
\alias{traces}
\alias{traces,ternaryFit-method}
\alias{stageCount}
\alias{stageCount,ternaryFit-method}
\alias{xSeed}
\alias{xSeed,ternaryFit-method}
\alias{inputParams}
\alias{inputParams,ternaryFit-method}
\alias{perturbationObj<-}
\alias{perturbationObj<-,ternaryFit-method}
\alias{steadyStateObj<-}
\alias{steadyStateObj<-,ternaryFit-method}
\alias{degreeObjMin<-}
\alias{degreeObjMin<-,ternaryFit-method}
\alias{graphObjMin<-}
\alias{graphObjMin<-,ternaryFit-method}
\alias{tableObjMin<-}
\alias{tableObjMin<-,ternaryFit-method}
\alias{newScore<-}
\alias{newScore<-,ternaryFit-method}
\alias{minScore<-}
\alias{minScore<-,ternaryFit-method}
\alias{finalTemperature<-}
\alias{finalTemperature<-,ternaryFit-method}
\alias{traces<-}
\alias{traces<-,ternaryFit-method}
\alias{stageCount<-}
\alias{stageCount<-,ternaryFit-method}
\alias{xSeed<-}
\alias{xSeed<-,ternaryFit-method}
\alias{inputParams<-}
\alias{inputParams<-,ternaryFit-method}

\title{Ternary Network Fit}
\description{This is a class representation of the output of the ternary
  network fitting algorithm implemented in the function \code{tnetfit}.}

\section{Creating Objects}{
  While one can create their own objects using the function
  \code{ternaryFit()}, this is highly discouraged. Typically this class
  is created by running the \code{tnetfit} function.
}
  
\section{Slots}{
  \describe{
    \item{\code{perturbationObj}:}{a matrix of perturbation
  experiments. Rows are genes and columns are experiments.}
    \item{\code{steadyStateObj}:}{a matrix of steady gene expression
  observations from a perturbation experiment. Rows are genes and
  columns are experiments.}
    \item{\code{degreeObjMin}:}{a vector containing the in-degree of
  each node in the fit achieving the minimum score}
    \item{\code{graphObjMin}:}{a matrix containing the parents of
  each node in the fit achieving the minimum score}
    \item{\code{tableObjMin}:}{a matrix containing the table in the fit
  achieving the minimum score} 
    \item{\code{newScore}:}{the most recent score}
    \item{\code{minScore}:}{the minimum score}
    \item{\code{finalTemperature}:}{the final value of the temperature
  parameter} 
    \item{\code{traces}:}{a dataframe contain the traces for 4 parameters}
    \item{\code{stageCount}:}{the number of stages}
    \item{\code{xSeed}:}{the random seed.}
    \item{\code{inputParams}:}{the ternaryFitParameters object used.}
  }
}

\section{Methods}{
  All named elements can be accessed and set in the standard way
  (e.g. \code{xSeed(object)} and \code{xSeed(object)<-}).
}
\author{Matthew N. McCall and Anthony Almudevar}
\seealso{
  \code{tnetpost}, \code{ternaryFitParameters-class}, \code{ternaryPost-class}.
  Almudevar A, McCall MN, McMurray H, Land H (2011). Fitting
  Boolean Networks from Steady State Perturbation Data, Statistical
  Applications in Genetics and Molecular Biology, 10(1): Article 47.
} 

\examples{
ssObj <- matrix(c(1,1,1,0,1,1,0,0,1),nrow=3)
pObj <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3)
rownames(ssObj) <- rownames(pObj) <- colnames(ssObj) <- colnames(pObj) <- c("Gene1","Gene2","Gene3")
tnfitObj <- tnetfit(ssObj, pObj)
class(tnfitObj)
}

\keyword{classes}