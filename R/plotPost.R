plotPost <- function(ternaryPost, threshold=0.5, type="interactive", ...){
  grph=graphPosterior(ternaryPost)
  el=which(grph[,-1]>threshold, arr.ind=TRUE)[,c(2,1)]
  if(type=="static") plot(graph.edgelist(el), ...)
  if(type=="interactive") tkplot(graph.edgelist(el), ...)
}
