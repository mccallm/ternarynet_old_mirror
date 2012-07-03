getAttractor <- function(tpost, perturbations, wildtype=TRUE, verbose=FALSE){

  fg <- perturbations$perturbed.genes
  fs <- perturbations$forced.states+2
  
  n.post <- dim(tableObjs(tpost))[3]
  n.gene <- dim(perturbationObj(tpost))[1] 
  n.sample <- dim(perturbationObj(tpost))[2] 

  attractor.list <- matrix(nrow=n.post, ncol=n.gene)
  for (i in c(1:n.post)) {
    tn.model <- list(table=tableObjs(tpost)[,,i], graph=graphObjs(tpost)[,,i], degree=degreeObjs(tpost)[i,])

    vec0 <- rep(2,length(tn.model$degree))
    vec0[fg] <- fs		
    junk <- attractor.distance.summary(vec0, tn.model, wildtype, fg, fs) 
    nc <- dim(junk$attractor)[2]	
    attractor.list[i,] <- summarize.attractor.2(junk$attractor[,2:nc])
    if(verbose) if(i %% 1000 == 0) cat(paste(i,"/",n.post,"\n"))
  }
  
  tmp <- condense.attractor(attractor.list)
  return(list(post.prob=tmp$p.attr, attractor.summary=tmp$uni.attr))
}
