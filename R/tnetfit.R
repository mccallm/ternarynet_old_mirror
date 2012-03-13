tnetfit <- function(steadyStateObj, perturbationObj, params=ternaryFitParameters(), xSeed=NA){

  ssObj <- as.vector(t(steadyStateObj))+2
  pObj <- as.integer(t(perturbationObj)+2)
  
  nGene <- as.integer(nrow(steadyStateObj))
  nExperiment <- as.integer(ncol(steadyStateObj))

  if(is.na(xSeed)) xSeed <- as.integer(floor(runif(1)*10000))
  set.seed(xSeed)
    
  tmp <- .Call("tnetfit", as.integer(perturbationType(params)), as.integer(scoreType(params)), as.integer(backupStage(params)), as.integer(maxStage(params)), as.integer(maxTransition(params)), epsilon(params), as.integer(beta0(params)), chi0(params), delta(params), as.integer(ne(params)), as.integer(m0(params)), as.integer(maxDegree(params)), pAddParent(params), pExchangeParent(params), as.integer(neighborDegree(params)), pNeighborhood(params), rho(params), nGene, nExperiment, edgePenalty(params), ssObj, pObj)

  traces <- data.frame(muTrace=tmp$muTrace, sigmaTrace=tmp$sigmaTrace, sigmaRhoTrace=tmp$sigmaRhoTrace, temperatureTrace=tmp$temperatureTrace)
  
  ternaryFit(perturbationObj, steadyStateObj, tmp$degreeObjMin, matrix(tmp$graphObjMin, nrow=nGene), matrix(tmp$tableObjMin, ncol=nGene), tmp$newScore, tmp$minScore, tmp$scTemperature, traces, tmp$stageCount, as.integer(xSeed), params)
  
}
