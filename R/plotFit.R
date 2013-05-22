plotFit <- function(ternaryFit, type="interactive", ...){
  grph=graphObjMin(ternaryFit)
  deg=apply(grph,2,function(x) sum(x>0))
  el=cbind(as.vector(grph)[as.vector(grph)>0],rep(1:ncol(grph),deg))
  el2=el
  el2[,1] <- colnames(grph)[el[,1]]
  el2[,2] <- colnames(grph)[el[,2]]
  if(type=="static") plot(graph.edgelist(el2), ...)
  if(type=="interactive") tkplot(graph.edgelist(el2), ...)
}

