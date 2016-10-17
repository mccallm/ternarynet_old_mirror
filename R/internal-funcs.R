
########## for R library

####################################################################

read.cc.app<-function(fname) {

	junk<- scan(fname,nlines=5,what="character",sep="\n")
	
	StartDate<-junk[1]
	EndDate<-junk[2]
	InputFile<-junk[3]
	OutputFile<-junk[4]
	Memo<-junk[5]

	junk<-scan(fname,skip=5)
	
	PerturbType <-junk[1]
	ScoreType <- junk[2]
	nBackupStage <- junk[3]
	MaxStage <- junk[4]
	MaxTransition <- junk[5]
	epsilon <- junk[6]
	beta <- junk[7]
	chi0 <- junk[8]
	delta <- junk[9]
	ne <- junk[10]
	m0 <- junk[11]
	MaxDegree <- junk[12]
	PAddParent <- junk[13]
	PExchangeParent <- junk[14]
	NeighborDegree <- junk[15]
	icount<-15	
	pNeighborhood <- NULL
	if (NeighborDegree > 1) {
		for (i in c(1:(NeighborDegree-1))) {
			icount<-icount+1
			pNeighborhood<-c(pNeighborhood, junk[icount])
		}
	}
	icount <- icount+1
	rho <- junk[icount]
	icount <- icount+1
	xSeed <- junk[icount]
	icount <- icount+1
	nGene <- junk[icount]
	icount <- icount+1
	nExp <- junk[icount]
	icount <- icount+1
	edgePenalty <- junk[icount]		

	DataResponse <- matrix(junk[c(1:(nGene*nExp))+icount],nGene,nExp,byrow=T)
	icount <- icount + (nGene*nExp)
	DataPerturbation <- matrix(junk[c(1:(nGene*nExp))+icount],nGene,nExp,byrow=T)
	icount <- icount + (nGene*nExp)

	icount <- icount+1
	NewScore <- junk[icount]
	icount <- icount+1	
	MinScore <- junk[icount]
	
	Fit.InDegree<-junk[(1:nGene)+icount]
	icount<-icount+nGene	
	net.graph<-NULL
	for (j in c(1:nGene)) {
		if (Fit.InDegree[j] > 0) {
			incr<-Fit.InDegree[j]
			net.graph<-append(net.graph,list(junk[c(1:incr)+icount]))
			icount<-icount+incr
		}
		if (Fit.InDegree[j] == 0) {
			net.graph<-append(net.graph,list(NULL))
		}
	}
	net.op<-NULL
	for (j in c(1:nGene)) {
		incr <- 3^Fit.InDegree[j]
		net.op<-append(net.op,list(junk[c(1:incr)+icount]))
		icount<-icount+incr
	}	
	net.table<-NULL
	
	MinFit.InDegree<-junk[(1:nGene)+icount]
	icount<-icount+nGene	
	min.net.graph<-NULL
	for (j in c(1:nGene)) {
		if (MinFit.InDegree[j] > 0) {
			incr<-MinFit.InDegree[j]
			min.net.graph<-append(min.net.graph,list(junk[c(1:incr)+icount]))
			icount<-icount+incr
		}
		if (MinFit.InDegree[j] == 0) {
			min.net.graph<-append(min.net.graph,list(NULL))
		}
	}
	min.net.op<-NULL
	for (j in c(1:nGene)) {
		incr <- 3^MinFit.InDegree[j]
		min.net.op<-append(min.net.op,list(junk[c(1:incr)+icount]))
		icount<-icount+incr
	}	


	icount<-icount+1
	FinalTemperature <- junk[icount]
	icount<-icount+1
	nStage <- junk[icount]	
	
	MuTrace<-junk[c(1:nStage)+icount]
	icount<-icount+nStage	
	SigmaTrace<-junk[c(1:nStage)+icount]
	icount<-icount+nStage	
	SigmaRhoTrace<-junk[c(1:nStage)+icount]
	icount<-icount+nStage	
	TemperatureTrace<-junk[c(1:nStage)+icount]
	icount<-icount+nStage	
	AcceptanceTrace<-junk[c(1:nStage)+icount]
	icount<-icount+nStage	

	fit.obj<-list(
		xSeed=xSeed,
        	nStage=nStage, 
        	PerturbType=PerturbType,
        	ScoreType=ScoreType,
        	nBackupStage=nBackupStage,
		nExp=nExp,
		nGene=nGene, 
		edgePenalty=edgePenalty,
        	MaxStage=MaxStage,
        	MaxTransition=MaxTransition,
        	epsilon=epsilon,
        	beta=beta,
        	chi0=chi0,
        	delta=delta,
        	ne=ne,
        	m0=m0,
        	m1=NULL,
        	m2=NULL,
		MaxDegree=MaxDegree,
		PAddParent=PAddParent,
		PExchangeParent=PExchangeParent,
		NeighborDegree=NeighborDegree,
		pNeighborhood=pNeighborhood,
        	rho=rho,
		FinalTemperature=FinalTemperature,		
		NewScore=NewScore,
        	MinScore=MinScore,
        	MinStage=NULL,
        	MinTransition=NULL,

		MinFit.InDegree=MinFit.InDegree,
		Fit.InDegree=Fit.InDegree,
		min.net.graph=min.net.graph, 
		min.net.op=min.net.op,
		net.graph=net.graph, 
		net.op=net.op,
		MuTrace=MuTrace,
	        SigmaTrace=SigmaTrace,
	        SigmaRhoTrace=SigmaRhoTrace,
	        TemperatureTrace=TemperatureTrace,
		AcceptanceTrace=AcceptanceTrace,
		DataResponse=DataResponse,
		DataPerturbation=DataPerturbation,

		Memo=Memo,
		StartDate=StartDate,
		EndDate=EndDate,
		InputFile=InputFile,
		OutputFile=OutputFile
	)

	return(fit.obj)

}

###################################################


read.cc.app.table.version<-function(fname) {

	junk<- scan(fname,nlines=5,what="character",sep="\n")
	
	StartDate<-junk[1]
	EndDate<-junk[2]
	InputFile<-junk[3]
	OutputFile<-junk[4]
	Memo<-junk[5]

	junk<-scan(fname,skip=5)
	
	PerturbType <-junk[1]
	ScoreType <- junk[2]
	nBackupStage <- junk[3]
	MaxStage <- junk[4]
	MaxTransition <- junk[5]
	epsilon <- junk[6]
	beta <- junk[7]
	chi0 <- junk[8]
	delta <- junk[9]
	ne <- junk[10]
	m0 <- junk[11]
	MaxDegree <- junk[12]
	PAddParent <- junk[13]
	PRemoveParent <- junk[14]
	PExchangeParent <- junk[15]
	NeighborDegree <- junk[16]
	icount<-16	
	pNeighborhood <- NULL
	if (NeighborDegree > 1) {
		for (i in c(1:(NeighborDegree-1))) {
			icount<-icount+1
			pNeighborhood<-c(pNeighborhood, junk[icount])
		}
	}
	icount <- icount+1
	rho <- junk[icount]
	icount <- icount+1
	xSeed <- junk[icount]
	icount <- icount+1
	nGene <- junk[icount]
	icount <- icount+1
	nExp <- junk[icount]
	icount <- icount+1
	edgePenalty <- junk[icount]		

	DataResponse <- matrix(junk[c(1:(nGene*nExp))+icount],nGene,nExp,byrow=T)
	icount <- icount + (nGene*nExp)
	DataPerturbation <- matrix(junk[c(1:(nGene*nExp))+icount],nGene,nExp,byrow=T)
	icount <- icount + (nGene*nExp)

	icount <- icount+1
	NewScore <- junk[icount]
	icount <- icount+1	
	MinScore <- junk[icount]
	
	Fit.InDegree<-junk[(1:nGene)+icount]
	icount<-icount+nGene	
	net.graph<-NULL
	for (j in c(1:nGene)) {
		if (Fit.InDegree[j] > 0) {
			incr<-Fit.InDegree[j]
			net.graph<-append(net.graph,list(junk[c(1:incr)+icount]))
			icount<-icount+incr
		}
		if (Fit.InDegree[j] == 0) {
			net.graph<-append(net.graph,list(NULL))
		}
	}
	net.op<-NULL
	for (j in c(1:nGene)) {
		incr <- 3^Fit.InDegree[j]
		net.op<-append(net.op,list(junk[c(1:incr)+icount]))
		icount<-icount+incr
	}	
	net.table<-NULL
	
	MinFit.InDegree<-junk[(1:nGene)+icount]
	icount<-icount+nGene	
	min.net.graph<-NULL
	for (j in c(1:nGene)) {
		if (MinFit.InDegree[j] > 0) {
			incr<-MinFit.InDegree[j]
			min.net.graph<-append(min.net.graph,list(junk[c(1:incr)+icount]))
			icount<-icount+incr
		}
		if (MinFit.InDegree[j] == 0) {
			min.net.graph<-append(min.net.graph,list(NULL))
		}
	}
	min.net.op<-NULL
	for (j in c(1:nGene)) {
		incr <- 3^MinFit.InDegree[j]
		min.net.op<-append(min.net.op,list(junk[c(1:incr)+icount]))
		icount<-icount+incr
	}	


	icount<-icount+1
	FinalTemperature <- junk[icount]
	icount<-icount+1
	nStage <- junk[icount]	
	
	MuTrace<-junk[c(1:nStage)+icount]
	icount<-icount+nStage	
	SigmaTrace<-junk[c(1:nStage)+icount]
	icount<-icount+nStage	
	SigmaRhoTrace<-junk[c(1:nStage)+icount]
	icount<-icount+nStage	
	TemperatureTrace<-junk[c(1:nStage)+icount]
	icount<-icount+nStage	
	
	fit.obj<-list(
		xSeed=xSeed,
        	nStage=nStage, 
        	PerturbType=PerturbType,
        	ScoreType=ScoreType,
        	nBackupStage=nBackupStage,
		nExp=nExp,
		nGene=nGene, 
		edgePenalty=edgePenalty,
        	MaxStage=MaxStage,
        	MaxTransition=MaxTransition,
        	epsilon=epsilon,
        	beta=beta,
        	chi0=chi0,
        	delta=delta,
        	ne=ne,
        	m0=m0,
        	m1=NULL,
        	m2=NULL,
		MaxDegree=MaxDegree,
		PAddParent=PAddParent,
		PRemoveParent=PRemoveParent,
		PExchangeParent=PExchangeParent,
		NeighborDegree=NeighborDegree,
		pNeighborhood=pNeighborhood,
        	rho=rho,
		FinalTemperature=FinalTemperature,		
		NewScore=NewScore,
        	MinScore=MinScore,
        	MinStage=NULL,
        	MinTransition=NULL,

		MinFit.InDegree=MinFit.InDegree,
		Fit.InDegree=Fit.InDegree,
		min.net.graph=min.net.graph, 
		min.net.op=min.net.op,
		net.graph=net.graph, 
		net.op=net.op,
		MuTrace=MuTrace,
	        SigmaTrace=SigmaTrace,
	        SigmaRhoTrace=SigmaRhoTrace,
	        TemperatureTrace=TemperatureTrace,
		AcceptanceTrace=NULL,
		DataResponse=DataResponse,
		DataPerturbation=DataPerturbation,

		Memo=Memo,
		StartDate=StartDate,
		EndDate=EndDate,
		InputFile=InputFile,
		OutputFile=OutputFile
	)

	return(fit.obj)

}



###################################################

read.cc.app.bool<-function(fname) {

	junk<- scan(fname,nlines=5,what="character",sep="\n")
	
	StartDate<-junk[1]
	EndDate<-junk[2]
	InputFile<-junk[3]
	OutputFile<-junk[4]
	Memo<-junk[5]

	junk<-scan(fname,skip=5)
	
	PerturbType <-junk[1]
	ScoreType <- junk[2]
	nBackupStage <- junk[3]
	MaxStage <- junk[4]
	MaxTransition <- junk[5]
	epsilon <- junk[6]
	beta <- junk[7]
	chi0 <- junk[8]
	delta <- junk[9]
	ne <- junk[10]
	m0 <- junk[11]
	MaxDegree <- junk[12]
	PAddParent <- junk[13]
	PExchangeParent <- junk[14]
	NeighborDegree <- junk[15]
	icount<-15	
	pNeighborhood <- NULL
	if (NeighborDegree > 1) {
		for (i in c(1:(NeighborDegree-1))) {
			icount<-icount+1
			pNeighborhood<-c(pNeighborhood, junk[icount])
		}
	}
	icount <- icount+1
	rho <- junk[icount]
	icount <- icount+1
	xSeed <- junk[icount]
	icount <- icount+1
	nGene <- junk[icount]
	icount <- icount+1
	nExp <- junk[icount]
	icount <- icount+1
	edgePenalty <- junk[icount]		
	icount <- icount+1
	cyclePenaltyN <- junk[icount]
	cyclePenaltyKnots <- junk[icount + c(1:(2*cyclePenaltyN+2))] 
	icount<-icount + (2*cyclePenaltyN+2)

	DataResponse <- matrix(junk[c(1:(nGene*nExp))+icount],nGene,nExp,byrow=T)
	icount <- icount + (nGene*nExp)
	DataPerturbation <- matrix(junk[c(1:(nGene*nExp))+icount],nGene,nExp,byrow=T)
	icount <- icount + (nGene*nExp)

	icount <- icount+1
	NewScore <- junk[icount]
	icount <- icount+1	
	MinScore <- junk[icount]
	
	Fit.InDegree<-junk[(1:nGene)+icount]
	icount<-icount+nGene	
	net.graph<-NULL
	for (j in c(1:nGene)) {
		if (Fit.InDegree[j] > 0) {
			incr<-Fit.InDegree[j]
			net.graph<-append(net.graph,list(junk[c(1:incr)+icount]))
			icount<-icount+incr
		}
		if (Fit.InDegree[j] == 0) {
			net.graph<-append(net.graph,list(NULL))
		}
	}
	net.op<-NULL
	for (j in c(1:nGene)) {
		incr <- 2^Fit.InDegree[j]
		net.op<-append(net.op,list(junk[c(1:incr)+icount]))
		icount<-icount+incr
	}	
	net.table<-NULL
	
	MinFit.InDegree<-junk[(1:nGene)+icount]
	icount<-icount+nGene	
	min.net.graph<-NULL
	for (j in c(1:nGene)) {
		if (MinFit.InDegree[j] > 0) {
			incr<-MinFit.InDegree[j]
			min.net.graph<-append(min.net.graph,list(junk[c(1:incr)+icount]))
			icount<-icount+incr
		}
		if (MinFit.InDegree[j] == 0) {
			min.net.graph<-append(min.net.graph,list(NULL))
		}
	}
	min.net.op<-NULL
	for (j in c(1:nGene)) {
		incr <- 2^MinFit.InDegree[j]
		min.net.op<-append(min.net.op,list(junk[c(1:incr)+icount]))
		icount<-icount+incr
	}	


	icount<-icount+1
	FinalTemperature <- junk[icount]
	icount<-icount+1
	nStage <- junk[icount]	
	
	MuTrace<-junk[c(1:nStage)+icount]
	icount<-icount+nStage	
	SigmaTrace<-junk[c(1:nStage)+icount]
	icount<-icount+nStage	
	SigmaRhoTrace<-junk[c(1:nStage)+icount]
	icount<-icount+nStage	
	TemperatureTrace<-junk[c(1:nStage)+icount]
	icount<-icount+nStage	
	
	fit.obj<-list(
		xSeed=xSeed,
        	nStage=nStage, 
        	PerturbType=PerturbType,
        	ScoreType=ScoreType,
        	nBackupStage=nBackupStage,
		nExp=nExp,
		nGene=nGene, 
		edgePenalty=edgePenalty,
		cyclePenaltyN=cyclePenaltyN,
		cyclePenaltyKnots=cyclePenaltyKnots,
        	MaxStage=MaxStage,
        	MaxTransition=MaxTransition,
        	epsilon=epsilon,
        	beta=beta,
        	chi0=chi0,
        	delta=delta,
        	ne=ne,
        	m0=m0,
        	m1=NULL,
        	m2=NULL,
		MaxDegree=MaxDegree,
		PAddParent=PAddParent,
		PExchangeParent=PExchangeParent,
		NeighborDegree=NeighborDegree,
		pNeighborhood=pNeighborhood,
        	rho=rho,
		FinalTemperature=FinalTemperature,		
		NewScore=NewScore,
        	MinScore=MinScore,
        	MinStage=NULL,
        	MinTransition=NULL,

		MinFit.InDegree=MinFit.InDegree,
		Fit.InDegree=Fit.InDegree,
		min.net.graph=min.net.graph, 
		min.net.op=min.net.op,
		net.graph=net.graph, 
		net.op=net.op,
		MuTrace=MuTrace,
	        SigmaTrace=SigmaTrace,
	        SigmaRhoTrace=SigmaRhoTrace,
	        TemperatureTrace=TemperatureTrace,
		AcceptanceTrace=NULL,
		DataResponse=DataResponse,
		DataPerturbation=DataPerturbation,

		Memo=Memo,
		StartDate=StartDate,
		EndDate=EndDate,
		InputFile=InputFile,
		OutputFile=OutputFile
	)

	return(fit.obj)

}



########################################################

read.cc.app.posterior <- function(fname) {
	
	junk<- scan(fname,nlines=5,what="character",sep="\n")
	
	StartDate<-junk[1]
	EndDate<-junk[2]
	FitFile<-junk[3]
	OutputPosteriorFile<-junk[4]
	OutputPosteriorMemo<-junk[5]

	junk<-scan(fname,skip=5)
	
	m0 <-junk[1]
	mdelta <- junk[2]
	parameter <- junk[3]
	xSeed <- junk[4]
	FinalTemperature <- junk[5]
	nGene <- junk[6]
	icount <- 6
	
	fit.score<-rep(NA, m0)
	fit.indegree.list<-NULL
	for (i in c(1:m0)) {fit.indegree.list<-append(fit.indegree.list,list(NULL))}
	fit.graph.list<-fit.indegree.list
	fit.op.list<-fit.indegree.list
	
	for (iii in 1:m0) {

		icount<-icount+1
		fit.score[iii] <- junk[icount]

		fit.indegree.list[[iii]]<-junk[(1:nGene)+icount]
		icount<-icount+nGene	

		fit.graph<-NULL
		for (j in c(1:nGene)) {
			if (fit.indegree.list[[iii]][j] > 0) {
				incr<- fit.indegree.list[[iii]][j]
				fit.graph<-append(fit.graph,list(junk[c(1:incr)+icount]))
				icount<-icount+incr
			}
			if (fit.indegree.list[[iii]][j] == 0) {
				fit.graph<-append(fit.graph,list(NULL))
			}
		}
		fit.graph.list[[iii]] <- fit.graph
		fit.op<-NULL
		for (j in c(1:nGene)) {
			incr <- 3^fit.indegree.list[[iii]][j]
			fit.op<-append(fit.op,list(junk[c(1:incr)+icount]))
			icount<-icount+incr
		}	
		fit.op.list[[iii]]<-fit.op
	}	
	
	posterior.obj<-list(

		StartDate=StartDate,
		EndDate=EndDate,
		FitFile=FitFile,
		OutputPosteriorFile=OutputPosteriorFile,
		OutputPosteriorMemo=OutputPosteriorMemo,
		m0 = m0,
		mdelta = mdelta,
		parameter =parameter,
		xSeed = xSeed,
		FinalTemperature = FinalTemperature,
		nGene = nGene,
		fit.score = fit.score,
		fit.indegree.list = fit.indegree.list,
		fit.graph.list = fit.graph.list,
		fit.op.list = fit.op.list
	)

	return(posterior.obj)
}

########################################################

read.cc.app.posterior.bool <- function(fname) {
	
	junk<- scan(fname,nlines=5,what="character",sep="\n")
	
	StartDate<-junk[1]
	EndDate<-junk[2]
	FitFile<-junk[3]
	OutputPosteriorFile<-junk[4]
	OutputPosteriorMemo<-junk[5]

	junk<-scan(fname,skip=5)
	
	m0 <-junk[1]
	mdelta <- junk[2]
	parameter <- junk[3]
	xSeed <- junk[4]
	FinalTemperature <- junk[5]
	nGene <- junk[6]
	icount <- 6
	
	fit.score<-rep(NA, m0)
	fit.indegree.list<-NULL
	for (i in c(1:m0)) {fit.indegree.list<-append(fit.indegree.list,list(NULL))}
	fit.graph.list<-fit.indegree.list
	fit.op.list<-fit.indegree.list
	
	for (iii in 1:m0) {

		icount<-icount+1
		fit.score[iii] <- junk[icount]

		fit.indegree.list[[iii]]<-junk[(1:nGene)+icount]
		icount<-icount+nGene	

		fit.graph<-NULL
		for (j in c(1:nGene)) {
			if (fit.indegree.list[[iii]][j] > 0) {
				incr<- fit.indegree.list[[iii]][j]
				fit.graph<-append(fit.graph,list(junk[c(1:incr)+icount]))
				icount<-icount+incr
			}
			if (fit.indegree.list[[iii]][j] == 0) {
				fit.graph<-append(fit.graph,list(NULL))
			}
		}
		fit.graph.list[[iii]] <- fit.graph
		fit.op<-NULL
		for (j in c(1:nGene)) {
			incr <- 2^fit.indegree.list[[iii]][j]
			fit.op<-append(fit.op,list(junk[c(1:incr)+icount]))
			icount<-icount+incr
		}	
		fit.op.list[[iii]]<-fit.op
	}	
	
	posterior.obj<-list(

		StartDate=StartDate,
		EndDate=EndDate,
		FitFile=FitFile,
		OutputPosteriorFile=OutputPosteriorFile,
		OutputPosteriorMemo=OutputPosteriorMemo,
		m0 = m0,
		mdelta = mdelta,
		parameter =parameter,
		xSeed = xSeed,
		FinalTemperature = FinalTemperature,
		nGene = nGene,
		fit.score = fit.score,
		fit.indegree.list = fit.indegree.list,
		fit.graph.list = fit.graph.list,
		fit.op.list = fit.op.list
	)

	return(posterior.obj)
}

####################################################################

ArrayToHash<-function(pIA, nOutcomes=3) {


	n<-length(pIA)
        tempi<-0

        for ( i in c(1:n) ) {tempi<-tempi+(nOutcomes^(i-1))*(pIA[i]-1)}

        tempi<-tempi+1

	return(tempi)

}

############################################################

f4p<-function(x) {sum(4^x)}

############################################################

condense.attractor<-function(m,thr=0) {

	nr<-dim(m)[1]
	nc<-dim(m)[2]
	mm<-m


	for (i in c(1:nc)) {mm[,i][m[,i]==9]<-2}
	mm<-mm+1
	code<-apply(mm,1,f4p)
	
	unicode<-unique(code) 

	p.attr<-rep(0,length(unicode)) 
	uni.attr<-matrix(0,length(unicode),nc)

	for (i in c(1:length(unicode))) {
		p.attr[i]<-mean(code==unicode[i])
		ind<-min(c(1:nr)[code==unicode[i]])
		uni.attr[i,]<-m[ind,]
	}

	ind<-rev(sort.list(p.attr))
	uni.attr<-as.matrix(uni.attr[ind,])
	p.attr<-p.attr[ind]
	uni.attr<-uni.attr[p.attr >= thr, ]
	p.attr<-p.attr[p.attr >= thr]


	return(list(p.attr=p.attr,uni.attr=uni.attr)) 

}

############################################################

summarize.attractor.calc<-function(m01) {

	m0<-m01[1]
	m1<-m01[2]

	ans<-0
	if ( (m0 == -1) & (m1 <= 0) ) {ans <- -1}
	if ( (m0 == -1) & (m1 == 1) ) {ans <- 9} 
	if ( (m0 >= 0) & (m1 == 1) ) {ans <- 1} 	

	return(ans)
}

############################################################


summarize.attractor.1<-function(m) {

	m<-as.matrix(m)
	attractor<-apply(m,1,mean)
	return(attractor) 

}

############################################################


summarize.attractor.2<-function(m) {

	m<-as.matrix(m)
	m0<-apply(m,1,min)
	m1<-apply(m,1,max)
	attractor<-apply(cbind(m0,m1),1,summarize.attractor.calc)
	return(attractor)	
	
}


############################################################



summarize.attractor.bool<-function(m) {

	m<-as.matrix(m)
	m0<-apply(m,1,min)
	m1<-apply(m,1,max)
	attractor<-1*(m0==1 & m1==1) + 9*(m0 < m1) 
	return(attractor)	
	
}


############################################################


score.attractor<-function(x) {
	
	obs<-x[1]
	mn<-x[2]
	mx<-x[3]

	ans<-1
	if ( (obs==1) & (mn==1) & (mx <= 2) ) {ans<-0}
	if ( (obs==2) & (mn==2) & (mx == 2) ) {ans<-0}
	if ( (obs==3) & (mn >= 2) & (mx == 3) ) {ans<-0}
	
	return(ans)
}


############################################################

apply.op<-function(vec,tn.model) {
	n<-length(vec)
	new.vec<-vec
	for (i in c(1:n)) {
		if (tn.model$degree[i]==0) {new.vec[i]<-2}
		if (tn.model$degree[i] > 0) {
			in.state<-vec[tn.model$graph[1:tn.model$degree[i],i]]
			tempi<-ArrayToHash(in.state, 3)
			new.vec[i]<-tn.model$table[tempi,i]
 		}
	}
	return(new.vec)
}

############################################################

attractor.distance.summary<-function(vec0, tn.model, wt, forced.genes=NULL,forced.states=NULL,nOutcomes=3) {

	vec.list<-cbind(NULL,vec0)

	list.pos<-0
	list.size<-1
	while (list.pos==0) {
		new.vec<-apply.op(vec.list[,list.size], tn.model)
		if (!wt) {new.vec[forced.genes]<-forced.states}
		for (i in c(1:list.size)) {
			if (sum(new.vec != vec.list[,i])==0) {list.pos<-i}
		}
		vec.list<-cbind(vec.list,new.vec)
		list.size<-list.size+1
	}
	
	
	attractor<-cbind(vec0, vec.list[,(list.pos+1):list.size]) - nOutcomes+1
	return(list(attractor=attractor,path=vec.list-nOutcomes+1))
}


#############################################################################################
##
## visible functions
##
#############################################################################################


attractor.2<-function(tn.model, model.obj, wt, nOutcomes=3) {

	n.gene<-dim(model.obj$initial.states)[1] 
	n.sample<-dim(model.obj$initial.states)[2] 

	attractor.matrix<-matrix(NA,n.gene,n.sample)	

	for (i in c(1:n.sample)) {
		fg<-model.obj$forced.gene.list[[i]]$forced.genes
		fs<-model.obj$forced.gene.list[[i]]$forced.states
		#vec0<-model.obj$initial.states[,i]
		vec0<-rep(2,n.gene)
		vec0[fg]<-fs		
		junk<-attractor.distance.summary(vec0, tn.model, wt, fg, fs) 
		nc<-dim(junk$attractor)[2]	
		if (nOutcomes==3) {attractor.matrix[,i]<-summarize.attractor.2(junk$attractor[,2:nc])}
		if (nOutcomes==2) {attractor.matrix[,i]<-summarize.attractor.bool(junk$attractor[,2:nc])}			
	}

	return(attractor.matrix)

}		

#######################################

attractor.2.detail<-function(net.op, net.graph, model.obj, wt) {

	n.gene<-dim(model.obj$initial.states)[1] 
	n.sample<-dim(model.obj$initial.states)[2] 

	junk.list<-NULL	

	for (i in c(1:n.sample)) {
		fg<-model.obj$forced.gene.list[[i]]$forced.genes
		fs<-model.obj$forced.gene.list[[i]]$forced.states
		vec0<-model.obj$initial.states[,i]
		junk<-attractor.distance.summary(vec0, model.obj, wt, fg, fs, model.obj$nOutcomes) 
		junk.list<-append(junk.list, list(junk))			
	}
	return(junk.list)
}		

#####################################

attractor.posterior.2<-function(posteriorObj, model.obj, wt ) {
	

	n.post <- dim(tableObjs(posteriorObj))[3]
	n.gene <- dim(perturbationObj(posteriorObj))[1] 
	n.sample <- dim(perturbationObj(posteriorObj))[2] 

	attractor.list<-NULL
	for (i in c(1:n.sample)) {attractor.list<-append(attractor.list,list(NULL))}

	for (i in c(1:n.sample)) {attractor.list[[i]]<-matrix(0,n.post,n.gene)}
	
	for (i in c(1:n.post)) {

		tn.model<-list(table=tableObjs(posteriorObj)[,,i], 
				graph=graphObjs(posteriorObj)[,,i], 
				degree=degreeObjs(posteriorObj)[i,]) 

		junk<-attractor.2(tn.model, model.obj, wt)
		for (j in c(1:n.sample)) {
			attractor.list[[j]][i,]<-junk[,j]
		}
	}
	
	return(attractor.list)

}

####################################################################


graph.posterior.marginal <- function( posteriorObj ) {


	n.post <- dim(tableObjs(posteriorObj))[3]
	n.gene <- dim(perturbationObj(posteriorObj))[1] 
	n.sample <- dim(perturbationObj(posteriorObj))[2] 


	parent.post.table<-matrix(0,n.gene,n.gene+1)	

	for (i in c(1:n.gene)) {
		for (j in c(1:n.post)) {

			gr<-graphObjs(posteriorObj)[,,j]			

			pp<-c(gr[,i],0,0)[1:2]

			if (pp[1]==0 & pp[2]==0) {
				parent.post.table[i,pp[1]+1]<-parent.post.table[i,pp[1]+1]+1	
			}
			if (pp[1]+pp[2] > 0) {
				parent.post.table[i,pp[1]+1]<-parent.post.table[i,pp[1]+1]+1	
				parent.post.table[i,pp[2]+1]<-parent.post.table[i,pp[2]+1]+1
			}
		}
	}
	
	parent.post.table <- parent.post.table/n.post		

	return(parent.post.table)

}
	
####################################################################

