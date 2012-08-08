plotFit <- function(ternaryFit, type="interactive", ...){
  grph=graphObjMin(ternaryFit)
  deg=apply(grph,2,function(x) sum(x>0))
  el=cbind(rep(1:ncol(grph),deg),as.vector(grph)[as.vector(grph)>0])
  if(type=="static") plot(graph.edgelist(el), ...)
  if(type=="interactive") tkplot(graph.edgelist(el), ...)
}

