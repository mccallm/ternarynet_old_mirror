tnetpost <- function(tfit, mdelta=as.integer(2000), msample=as.integer(2000), temperatureScale=1.0, xSeed=NA){

  steadyStateObj <- as.matrix(steadyStateObj(tfit))+2
  perturbationObj <- matrix(as.integer(as.matrix(perturbationObj(tfit))+2), ncol=ncol(steadyStateObj))

  nGene <- as.integer(nrow(steadyStateObj))
  nExperiment <- as.integer(ncol(steadyStateObj))

  if(is.na(xSeed)) xSeed <- as.integer(floor(runif(1)*10000))
  set.seed(xSeed)
  
  inputParams <- inputParams(tfit)

  obj <- .Call("tnetpost", msample, mdelta, temperatureScale, xSeed(tfit), scoreType(inputParams), perturbationType(inputParams), backupStage(inputParams), maxStage(inputParams), maxTransition(inputParams), epsilon(inputParams), beta0(inputParams), chi0(inputParams), delta(inputParams), rho(inputParams), ne(inputParams), maxDegree(inputParams), pAddParent(inputParams), pExchangeParent(inputParams), neighborDegree(inputParams), pNeighborhood(inputParams), nGene, nExperiment, edgePenalty(inputParams), m0(inputParams), steadyStateObj, perturbationObj, newScore(tfit), minScore(tfit), degreeObjMin(tfit), graphObjMin(tfit), tableObjMin(tfit), finalTemperature(tfit))

  ternaryPost(perturbationObj=perturbationObj-2, steadyStateObj=steadyStateObj-2, scores=obj$scores, degreeObjs=matrix(obj$degreeObjs, nrow=length(obj$scores), byrow=TRUE), graphObjs=array(obj$graphObjs, dim=c(nGene,nGene,length(obj$scores))), tableObjs=array(obj$tableObjs, dim=c(length(obj$tableObjs)/(nGene*length(obj$scores)),nGene,length(obj$scores))), inputParams=inputParams)

}
