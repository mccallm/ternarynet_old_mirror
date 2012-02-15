#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "math.h"
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "tnetfuncs.h"

using namespace std;

extern "C" {
  SEXP tnetfit(SEXP RperturbationType, SEXP RscoreType, SEXP RbackupStage, SEXP RmaxStage, SEXP RmaxTransition, SEXP Repsilon, SEXP Rbeta, SEXP Rchi0, SEXP Rdelta, SEXP Rne, SEXP Rm0, SEXP RmaxDegree, SEXP RpAddParent, SEXP RpExchangeParent, SEXP RneighborDegree, SEXP RpNeighborhood, SEXP Rrho, SEXP RnGene, SEXP RnExperiment, SEXP RedgePenalty, SEXP RsteadyStateObj, SEXP RperturbationObj){

    SEXP Routput, Rnames;
    int perturbationType, scoreType, backupStage, maxStage, maxTransition;
    double epsilon;
    int beta;
    double chi0, delta;
    int ne, m0, maxDegree;
    double pAddParent, pExchangeParent;
    int neighborDegree;
    double rho;
    int nGene, nExperiment;
    double edgePenalty;

    perturbationType = INTEGER(coerceVector(RperturbationType, INTSXP))[0];
    scoreType = INTEGER(coerceVector(RscoreType, INTSXP))[0];
    backupStage = INTEGER(coerceVector(RbackupStage, INTSXP))[0];
    maxStage = INTEGER(coerceVector(RmaxStage, INTSXP))[0];
    maxTransition = INTEGER(coerceVector(RmaxTransition, INTSXP))[0];
    epsilon = REAL(Repsilon)[0];
    beta = INTEGER(coerceVector(Rbeta, INTSXP))[0];
    chi0 = REAL(Rchi0)[0];
    delta = REAL(Rdelta)[0];
    ne = INTEGER(coerceVector(Rne, INTSXP))[0];
    m0 = INTEGER(coerceVector(Rm0, INTSXP))[0];
    maxDegree = INTEGER(coerceVector(RmaxDegree, INTSXP))[0];
    pAddParent = REAL(RpAddParent)[0];
    pExchangeParent = REAL(RpExchangeParent)[0];
    neighborDegree = INTEGER(coerceVector(RneighborDegree, INTSXP))[0];
    rho = REAL(Rrho)[0];
    nGene = INTEGER(coerceVector(RnGene, INTSXP))[0];
    nExperiment = INTEGER(coerceVector(RnExperiment, INTSXP))[0];
    edgePenalty = REAL(RedgePenalty)[0];
    
    const double *pNeighborhood = REAL(RpNeighborhood);

    const double *steadyStateObj = REAL(RsteadyStateObj);

    /*    
	  cout << "steadyStateObj" << endl;
	  for (int i = 0; i < (I*J); i++) cout << steadyStateObj[i] << " ";
	  cout << endl;
    */    

    const int *perturbationObj = INTEGER(RperturbationObj);
    
    /*
      cout << "perturbationObj" << endl;
      for (int i = 0; i < (I*J); i++) cout << perturbationObj[i] << " ";
      cout << endl;
    */

    int nOutcomes=3, l0, stage_count, tempi, jtrans, nTransition;
    int tableWidth = powi(nOutcomes,maxDegree);
    bool sigma_zero_flag;
    double sigmaci_toler = 0.0000001;
    double newScore, oldScore, difference, minScore = 999999999;
    double temperature, acceptProb;
    int proposalDegree, iProposal;
    int *degreeObj = new int[nGene+1];
    int *graphObj = new int[nGene*nGene];
    int *tableObj = new int[nGene*tableWidth]; 
    int *degreeObjMin = new int[nGene+1];
    int *graphObjMin = new int[nGene*nGene];
    int *tableObjMin = new int[nGene*tableWidth];
    int *proposalGene = new int[neighborDegree+1];
    int *oldDegree = new int[neighborDegree+1];
    int *newDegree = new int[1];
    int *oldGraphComponent = new int[neighborDegree*nGene];
    int *newGraphComponent = new int[nGene+1];
    int *oldTableComponent = new int[neighborDegree*tableWidth]; 
    int *newTableComponent = new int[tableWidth+1]; 
    double *currentValues = new double[ne+1];  
    double *muTrace = new double[maxStage+1];
    double *sigmaTrace = new double[maxStage+1];
    double *sigmaRhoTrace = new double[maxStage+1];
    double *temperatureTrace = new double[maxStage+1];
    double sigma0, sigmaci, sigmaRaw, muci, muci_old, convergence;
    int acceptCount;
    double u;

    // cout << "set intial state" << endl;
    
    // set random seed
    GetRNGstate();

    // set initial state
    for (int i=1; i <= nGene; ++i) degreeObj[i]=0;
    for (int i=0; i < (nGene*tableWidth); ++i) tableObj[i]=2;
    for (int i=0; i < (nGene*nGene); ++i) graphObj[i]=0;

    // cout << "estimate initial temp" << endl;

    // estimate initial temperature
    temperature = initialTemp(chi0, m0, nOutcomes, maxDegree, nGene, nExperiment, edgePenalty, pAddParent, pExchangeParent, neighborDegree, pNeighborhood, degreeObj, graphObj, tableObj, steadyStateObj, perturbationObj); 
    
    // cout << "oldScore" << endl;

    oldScore = AttractorDistanceForced(nGene, nOutcomes, maxDegree, nExperiment, edgePenalty, tableObj, graphObj, degreeObj, steadyStateObj, perturbationObj); 
    l0 = 0;
    sigma_zero_flag = false;	
    stage_count = 0;	
    
    // cout << "begin main loop" << endl;

    // being main loop
    while (l0 < beta && !sigma_zero_flag && stage_count < maxStage) {	
      nTransition = 0;
      jtrans = 0;
      acceptCount=0;
      while (jtrans < ne && nTransition < maxTransition) {
	// generate proposal 
	proposalDegree = randomNeighborDegree(neighborDegree, pNeighborhood);
	for (iProposal = 1; iProposal <= proposalDegree; ++iProposal) {
	  proposalGene[iProposal] = (int)floor(unif_rand() * nGene) + 1; //add (int)floor() 01-18-2012

	  oldDegree[iProposal] = degreeObj[proposalGene[iProposal]];
	  newDegree[0] = degreeObj[proposalGene[iProposal]];
	  for (int j=1; j <= oldDegree[iProposal]; ++j) {
	    oldGraphComponent[(iProposal-1)*nGene+j-1]=graphObj[nGene*(proposalGene[iProposal]-1)+j-1];
	    newGraphComponent[j]=graphObj[nGene*(proposalGene[iProposal]-1)+j-1];
	  }				
	  for (int j=1; j <= powi(nOutcomes, oldDegree[iProposal]); ++j) {
	    oldTableComponent[(iProposal-1)*tableWidth+j-1]=tableObj[(proposalGene[iProposal]-1)*tableWidth+j-1];
	    newTableComponent[j]=tableObj[(proposalGene[iProposal]-1)*tableWidth+j-1]; 		
	  }
	  PerturbGene(proposalGene[iProposal], nGene, maxDegree, nOutcomes, 
		      pAddParent, pExchangeParent, newDegree, newTableComponent, newGraphComponent);
	  degreeObj[proposalGene[iProposal]] = newDegree[0];
	  for (int j=1; j <= newDegree[0]; ++j) {
	    graphObj[nGene*(proposalGene[iProposal]-1)+j-1]=newGraphComponent[j];}
	  for (int j=1; j <= powi(nOutcomes, newDegree[0]); ++j) {
	    tableObj[(proposalGene[iProposal]-1)*tableWidth+j-1] = newTableComponent[j];}
	}
	// calculate new score
	newScore = AttractorDistanceForced(nGene, nOutcomes, maxDegree, nExperiment, edgePenalty, tableObj, graphObj, degreeObj, steadyStateObj, perturbationObj); 
	// accept or reject proposal 
	difference = newScore-oldScore;
	if (difference <= 0) {
	  ++acceptCount;      				
	}
	else{
	  if (temperature != 0) {
	    acceptProb = exp((-difference)/temperature); 
	    u = uniformdist();				
	    if (acceptProb < u){
	      newScore = oldScore;
	      for (iProposal = proposalDegree; iProposal >= 1; --iProposal) {
		degreeObj[proposalGene[iProposal]] = oldDegree[iProposal];
		for (int j=1; j <= oldDegree[iProposal]; ++j) {
		  graphObj[(proposalGene[iProposal]-1)*nGene+j-1] 
		    = oldGraphComponent[(iProposal-1)*nGene+j-1];}
		for (int j=1; j <= powi(nOutcomes, oldDegree[iProposal]); ++j) { 
		  tableObj[(proposalGene[iProposal]-1)*tableWidth+j-1] 
		    = oldTableComponent[(iProposal-1)*tableWidth+j-1];}
	      }
	    } else
	      {
		++acceptCount;
	      }
	  }    		
	}
	++jtrans;
	currentValues[jtrans] = newScore;
	if (newScore > oldScore) ++nTransition;
	// capture minimum, if needed 			
	if (newScore < minScore) {
	  minScore = newScore;
	  for (int j=1; j <= nGene; ++j) degreeObjMin[j] = degreeObj[j];
	  for (int j=0; j < (nGene*nGene); ++j) graphObjMin[j] = graphObj[j];
	  for (int j=0; j < (nGene*tableWidth); ++j) tableObjMin[j] = tableObj[j];
	}
	oldScore = newScore;
      } // while jtrans
      sigmaRaw = stagevariance(currentValues, jtrans);
      muci = stagemean(currentValues, jtrans);
      if (stage_count == 1) {
	sigma0 = sqrt(sigmaRaw);
	muci_old = muci;
      } else
	{
	  sigmaci = rho*sigmaRaw + (1-rho)*sigmaci;
	  if (sigmaci < sigmaci_toler) { 
	    sigma_zero_flag = true;      
	  }
	  else{
	    convergence = (muci-muci_old)/sigma0; 
	    if ( convergence < epsilon && convergence > -epsilon ) ++l0;
	    muci_old = muci;
	    temperature = temperature/(1+(temperature*log(1+delta)/(3*sigmaci)));
	  }
	}
      ++stage_count;
      muTrace[stage_count]=muci;
      sigmaTrace[stage_count]=sigmaRaw;
      sigmaRhoTrace[stage_count]=sigmaci;
      temperatureTrace[stage_count]=temperature;
    } // while
    
    PROTECT(Routput = allocVector(VECSXP, 14));
    
    SEXP RdegreeObj;
    PROTECT(RdegreeObj = allocVector(INTSXP, nGene));
    for (int i = 0; i < nGene; i++)
      INTEGER(RdegreeObj)[i] = degreeObj[i+1];
    SET_VECTOR_ELT(Routput, 0, RdegreeObj);
    
    SEXP RgraphObj;
    PROTECT(RgraphObj = allocVector(INTSXP, nGene*nGene));
    for (int i = 0; i < nGene; i++)
      for (int j = 0; j < nGene; j++){
	if(j < degreeObj[i+1]) { 
	  INTEGER(RgraphObj)[j + nGene * i] = graphObj[j + nGene * i]; 
	} else 
	  {
	    INTEGER(RgraphObj)[j + nGene * i] = 0; 
	}
      } 
    SET_VECTOR_ELT(Routput, 1, RgraphObj);
    
    SEXP RtableObj;
    PROTECT(RtableObj = allocVector(INTSXP, tableWidth*nGene));
    for (int i = 0; i < nGene; i++){
      for (int j = 0; j < tableWidth; j++){
	if(j < powi(nOutcomes, degreeObj[i+1])){
	  INTEGER(RtableObj)[j + tableWidth*i] = tableObj[j + tableWidth * i];
	} else
	  {
	    INTEGER(RtableObj)[j + tableWidth*i] = 0;
	  }
      }
    }
    SET_VECTOR_ELT(Routput, 2, RtableObj);
    
    SEXP RdegreeObjMin;
    PROTECT(RdegreeObjMin = allocVector(INTSXP, nGene));
    for (int i = 0; i < nGene; i++)
      INTEGER(RdegreeObjMin)[i] = degreeObjMin[i+1];
    SET_VECTOR_ELT(Routput, 3, RdegreeObjMin);
    
    SEXP RgraphObjMin;
    PROTECT(RgraphObjMin = allocVector(INTSXP, nGene*nGene));
    for (int i = 0; i < nGene; i++)
      for (int j = 0; j < nGene; j++){
	if(j < degreeObjMin[i+1]) { 
	  INTEGER(RgraphObjMin)[j + nGene * i] = graphObjMin[j + nGene * i]; 
	} else 
	  {
	  INTEGER(RgraphObjMin)[j + nGene * i] = 0; 
	}
      } 
    SET_VECTOR_ELT(Routput, 4, RgraphObjMin);

    SEXP RtableObjMin;
    PROTECT(RtableObjMin = allocVector(INTSXP, tableWidth*nGene));
    for (int i = 0; i < nGene; i++)
      for (int j = 0; j < tableWidth; j++){
	if(j < powi(nOutcomes, degreeObjMin[i+1])){
	  INTEGER(RtableObjMin)[j + tableWidth*i] = tableObjMin[j + tableWidth * i];
	} else
	  {
	    INTEGER(RtableObjMin)[j + tableWidth*i] = 0;
	  }
      }
    SET_VECTOR_ELT(Routput, 5, RtableObjMin);
    
    SEXP RmuTrace;
    PROTECT(RmuTrace = allocVector(REALSXP, stage_count));
    for (int i = 0; i < stage_count; i++)
      REAL(RmuTrace)[i] = muTrace[i+1];
    SET_VECTOR_ELT(Routput, 6, RmuTrace);
    
    SEXP RsigmaTrace;
    PROTECT(RsigmaTrace = allocVector(REALSXP, stage_count));
    for (int i = 0; i < stage_count; i++)
      REAL(RsigmaTrace)[i] = sigmaTrace[i+1];
    SET_VECTOR_ELT(Routput, 7, RsigmaTrace);
    
    SEXP RsigmaRhoTrace;
    PROTECT(RsigmaRhoTrace = allocVector(REALSXP, stage_count));
    for (int i = 0; i < stage_count; i++)
      REAL(RsigmaRhoTrace)[i] = sigmaRhoTrace[i+1];
    SET_VECTOR_ELT(Routput, 8, RsigmaRhoTrace);
    
    SEXP RtemperatureTrace;
    PROTECT(RtemperatureTrace = allocVector(REALSXP, stage_count));
    for (int i = 0; i < stage_count; i++)
      REAL(RtemperatureTrace)[i] = temperatureTrace[i+1];
    SET_VECTOR_ELT(Routput, 9, RtemperatureTrace);
    
    SET_VECTOR_ELT(Routput, 10, ScalarReal(newScore));

    SET_VECTOR_ELT(Routput, 11, ScalarReal(minScore));

    SET_VECTOR_ELT(Routput, 12, ScalarInteger(stage_count));

    SET_VECTOR_ELT(Routput, 13, ScalarReal(temperatureTrace[stage_count]));

    PROTECT(Rnames= allocVector(STRSXP,14));
    SET_STRING_ELT(Rnames,0,mkChar("degreeObj"));
    SET_STRING_ELT(Rnames,1,mkChar("graphObj"));
    SET_STRING_ELT(Rnames,2,mkChar("tableObj"));
    SET_STRING_ELT(Rnames,3,mkChar("degreeObjMin"));
    SET_STRING_ELT(Rnames,4,mkChar("graphObjMin"));
    SET_STRING_ELT(Rnames,5,mkChar("tableObjMin"));
    SET_STRING_ELT(Rnames,6,mkChar("muTrace"));
    SET_STRING_ELT(Rnames,7,mkChar("sigmaTrace"));
    SET_STRING_ELT(Rnames,8,mkChar("sigmaRhoTrace"));
    SET_STRING_ELT(Rnames,9,mkChar("temperatureTrace"));
    SET_STRING_ELT(Rnames,10,mkChar("newScore"));
    SET_STRING_ELT(Rnames,11,mkChar("minScore"));
    SET_STRING_ELT(Rnames,12,mkChar("stageCount"));
    SET_STRING_ELT(Rnames,13,mkChar("scTemperature"));
    setAttrib(Routput, R_NamesSymbol, Rnames);
    
    UNPROTECT(12);

    delete[] degreeObj;
    delete[] graphObj;
    delete[] tableObj;
    delete[] degreeObjMin;
    delete[] graphObjMin;
    delete[] tableObjMin;
    delete[] proposalGene;
    delete[] oldDegree;
    delete[] newDegree;
    delete[] oldGraphComponent;
    delete[] newGraphComponent;
    delete[] oldTableComponent;
    delete[] newTableComponent;
    delete[] currentValues;
    delete[] muTrace;
    delete[] sigmaTrace;
    delete[] sigmaRhoTrace;
    delete[] temperatureTrace;

    PutRNGstate();

    return Routput;
  }
}
