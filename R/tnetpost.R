tnetpost <- function(tfit, mdelta=as.integer(10000), msample=as.integer(2000), temperatureScale=1.0, xSeed=NA){

  ssObj <- as.vector(t(steadyStateObj(tfit))+2)
  pObj <- as.integer(t(perturbationObj(tfit))+2)

  nGene <- as.integer(nrow(steadyStateObj(tfit)))
  nExperiment <- as.integer(ncol(steadyStateObj(tfit)))

  if(is.na(xSeed)) xSeed <- as.integer(floor(runif(1)*10000))
  set.seed(xSeed)
  
  inputParams <- inputParams(tfit)

  obj <- .Call("tnetpost", msample, mdelta, temperatureScale, xSeed(tfit), scoreType(inputParams), perturbationType(inputParams), backupStage(inputParams), maxStage(inputParams), maxTransition(inputParams), epsilon(inputParams), beta0(inputParams), chi0(inputParams), delta(inputParams), rho(inputParams), ne(inputParams), maxDegree(inputParams), pAddParent(inputParams), pExchangeParent(inputParams), neighborDegree(inputParams), pNeighborhood(inputParams), nGene, nExperiment, edgePenalty(inputParams), m0(inputParams), ssObj, pObj, newScore(tfit), minScore(tfit), degreeObjMin(tfit), graphObjMin(tfit), tableObjMin(tfit), finalTemperature(tfit))

  ternaryPost(perturbationObj=perturbationObj(tfit), steadyStateObj=steadyStateObj(tfit), scores=obj$scores, degreeObjs=matrix(obj$degreeObjs, nrow=length(obj$scores), byrow=TRUE), graphObjs=array(obj$graphObjs, dim=c(nGene,nGene,length(obj$scores))), tableObjs=array(obj$tableObjs, dim=c(length(obj$tableObjs)/(nGene*length(obj$scores)),nGene,length(obj$scores))), inputParams=inputParams)

}
