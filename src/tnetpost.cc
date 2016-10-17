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
  SEXP tnetpost(SEXP Rm0, SEXP Rmdelta, SEXP Rparameter, SEXP RxSeedFit, SEXP RscoreType, SEXP RperturbationType, SEXP RbackupStage, SEXP RmaxStage, SEXP RmaxTransition, SEXP Repsilon, SEXP Rbeta, SEXP Rchi0, SEXP Rdelta, SEXP Rrho, SEXP Rne, SEXP RmaxDegree, SEXP RpAddParent, SEXP RpExchangeParent, SEXP RneighborDegree, SEXP RpNeighborhood, SEXP RnGene, SEXP RnExperiment, SEXP RedgePenalty, SEXP Rm0Fit, SEXP RsteadyStateObj, SEXP RperturbationObj, SEXP RnewScore, SEXP RminScore, SEXP RdegreeObjMin, SEXP RgraphObjMin, SEXP RtableObjMin, SEXP RfinalTemperature){

    SEXP Routput, Rnames;

    int ipause;
	
    int m0fit;
    double newScore, oldScore, minScore, finalScore;
    int nOutcomes = 3;
    int l0, stage_count, tempi, jtrans, nTransition;
    double difference;
    double temperature, acceptProb;
    int acceptCount;
    int proposalDegree, iProposal;
    double u;
    int m0, mdelta, xSeedFit;
    double parameter, finalTemperature; 
    int perturbationType, scoreType, backupStage;
    int maxStage, maxTransition;
    double epsilon, chi0, delta, rho; 
    int beta;
    int nGene, nExperiment;
    double edgePenalty; 
    int ne, maxDegree;
    double pAddParent, pExchangeParent;
    int neighborDegree; 
    
    m0 = INTEGER(coerceVector(Rm0, INTSXP))[0];
    mdelta = INTEGER(coerceVector(Rmdelta, INTSXP))[0];
    parameter = REAL(coerceVector(Rparameter, REALSXP))[0];
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
    newScore = REAL(coerceVector(RnewScore, REALSXP))[0];
    minScore = REAL(coerceVector(RminScore, REALSXP))[0];
    finalTemperature = REAL(coerceVector(RfinalTemperature, REALSXP))[0];
    xSeedFit = INTEGER(coerceVector(RxSeedFit, INTSXP))[0];
    
    const double *pNeighborhood = REAL(RpNeighborhood);

    Rprintf("Single value parameters loaded \n");

    // set random seed
    GetRNGstate();

    int tableWidth = powi(nOutcomes,maxDegree);
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

    const double *steadyStateObj = REAL(RsteadyStateObj);
    const int *perturbationObj = INTEGER(RperturbationObj);

    Rprintf("Perturbation matrices loaded \n");

    int *degreeObj = new int[nGene+1];
    RdegreeObjMin = coerceVector(RdegreeObjMin,INTSXP);
    for(int i=0; i < nGene; i++) degreeObj[i+1] = INTEGER(RdegreeObjMin)[i];

    int *graphObj = new int[nGene*nGene];
    RgraphObjMin = coerceVector(RgraphObjMin,INTSXP);
    for(int i=0; i < nGene; i++)
      for(int j=0; j < nGene; j++)
	graphObj[i*nGene + j] = INTEGER(RgraphObjMin)[i*nGene + j]; 

    int *tableObj = new int[nGene*tableWidth];
    RtableObjMin = coerceVector(RtableObjMin,INTSXP);
    for(int i=0; i < nGene; i++)
      for(int j=0; j < tableWidth; j++) 
	tableObj[i*tableWidth + j] = INTEGER(RtableObjMin)[i*tableWidth + j];

    Rprintf("Fit parameters loaded \n");

    //cout << "BP1" << endl;
    //cin >> ipause >> endl;


    double *myscores = new double[m0]; // debug 01-28-2012 m0*mdelta -> m0
    int myscores_count=0;

    //cout << "BP2" << endl;
    //cin >> ipause >> endl;

    int *degreeObjs = new int[m0*nGene]; // debug 01-28-2012 m0*mdelta -> m0

    //cout << "BP3" << endl;
    //cin >> ipause >> endl;

    int *graphObjs = new int[m0*nGene*nGene]; // debug 01-28-2012 m0*mdelta -> m0

    //cout << "BP4" << endl;
    //cin >> ipause >> endl;

    int *tableObjs = new int[m0*nGene*tableWidth]; // debug 01-28-2012 m0*mdelta -> m0

    //cout << "BP5" << endl;
    //cin >> ipause >> endl;

    Rprintf("Output parameters created \n");

    // start sampling

    oldScore = AttractorDistanceForced(nGene, nOutcomes, maxDegree, nExperiment, edgePenalty, tableObj, graphObj, degreeObj, steadyStateObj, perturbationObj); 

    // begin main loop

    bool minDump = true;
	
    temperature = parameter*finalTemperature;

    Rprintf("Starting to sample \n");

    for (int outerLoop = 1; outerLoop <= m0; ++outerLoop) {
      		
      acceptCount=0;
      for (int innerLoop = 1; innerLoop <= mdelta; ++innerLoop) {

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
	  PerturbGene(proposalGene[iProposal], nGene, maxDegree, nOutcomes, pAddParent, pExchangeParent, newDegree, newTableComponent, newGraphComponent);
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
	  minDump = true; 				
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
		minDump = true; 	
	      }
	  }    		
	}
	oldScore = newScore;
		
	// check for new minmimum
	/*		
	if ( (newScore < minScore) && minDump) {
	  myscores[myscores_count] = newScore;
	  
	  for (int i=0; i < nGene; ++i)
	    degreeObjs[myscores_count*nGene + i] = degreeObj[i+1];
	  
	  for (int i=0; i < nGene; ++i) {
	    for (int j=0; j < degreeObj[i+1]; ++j) {
	      graphObjs[myscores_count*nGene*nGene + nGene*i + j] = graphObj[nGene*i+j];
	    }
	  }
	  
	  for (int i=0; i < nGene; ++i) {
	    for (int j=0; j < powi(nOutcomes, degreeObj[i+1]); ++j) { 
	      tableObjs[myscores_count*nGene*tableWidth + tableWidth*i + j] = tableObj[tableWidth*i+j];
	    }
	  }
		
	  myscores_count++;
	  minDump = false;			
	}
	*/		
      } // for innerLoop

      if ((outerLoop % 10) == 0) Rprintf("Sample %d of %d \n", outerLoop, m0);
      
      // output model
      myscores[myscores_count] = newScore;
      for (int i=0; i < nGene; ++i)
      	degreeObjs[myscores_count*nGene + i] = degreeObj[i+1];
      
      for (int i=0; i < nGene; ++i) {
      	for (int j=0; j < degreeObj[i+1]; ++j) {
      	  graphObjs[myscores_count*nGene*nGene + nGene*i + j] = graphObj[nGene*i+j];
      	}
      }
      
      for (int i=0; i < nGene; ++i) {
      	for (int j=0; j < powi(nOutcomes, degreeObj[i+1]); ++j) { 
      	  tableObjs[myscores_count*nGene*tableWidth + tableWidth*i + j] = tableObj[tableWidth*i+j];
      	}
      }
      myscores_count++;
    } // for outerLoop

    PROTECT(Routput = allocVector(VECSXP, 4));
    
    SEXP Rmyscores;
    PROTECT(Rmyscores = allocVector(REALSXP, myscores_count));
    for (int i=0; i < myscores_count; i++)
      REAL(Rmyscores)[i] = myscores[i];
    SET_VECTOR_ELT(Routput, 0, Rmyscores);

    SEXP RdegreeObjs;
    PROTECT(RdegreeObjs = allocVector(INTSXP, myscores_count*nGene));
    for (int k = 0; k < myscores_count; k++)
      for (int i = 0; i < nGene; i++)
	INTEGER(RdegreeObjs)[k*nGene + i] = degreeObjs[k*nGene + i];
    SET_VECTOR_ELT(Routput, 1, RdegreeObjs);
    
    SEXP RgraphObjs;
    PROTECT(RgraphObjs = allocVector(INTSXP, myscores_count*nGene*nGene));
    for (int k = 0; k < myscores_count; k++)
      for (int i = 0; i < nGene; i++)
	for (int j = 0; j < nGene; j++){
	  if(j < degreeObjs[k*nGene + i]) { 
	    INTEGER(RgraphObjs)[k*nGene*nGene + i*nGene + j] = graphObjs[k*nGene*nGene + i*nGene + j]; 
	  } else 
	    {
	      INTEGER(RgraphObjs)[k*nGene*nGene + i*nGene + j] = 0; 
	    }
	}
    SET_VECTOR_ELT(Routput, 2, RgraphObjs);
    
    SEXP RtableObjs;
    PROTECT(RtableObjs = allocVector(INTSXP, myscores_count*nGene*tableWidth));
    for (int k = 0; k < myscores_count; k++)
      for (int i = 0; i < nGene; i++)
	for (int j = 0; j < tableWidth; j++){
	  if(j < powi(nOutcomes, degreeObjs[k*nGene + i])){
	    INTEGER(RtableObjs)[k*nGene*tableWidth + i*tableWidth + j] = tableObjs[k*nGene*tableWidth + i*tableWidth + j];
	  } else
	    {
	      INTEGER(RtableObjs)[k*nGene*tableWidth + i*tableWidth + j] = 0;
	    }
	}
    SET_VECTOR_ELT(Routput, 3, RtableObjs);

    PROTECT(Rnames= allocVector(STRSXP,4));
    SET_STRING_ELT(Rnames,0,mkChar("scores"));
    SET_STRING_ELT(Rnames,1,mkChar("degreeObjs"));
    SET_STRING_ELT(Rnames,2,mkChar("graphObjs"));
    SET_STRING_ELT(Rnames,3,mkChar("tableObjs"));
    setAttrib(Routput, R_NamesSymbol, Rnames);
    
    UNPROTECT(6);

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
    delete[] degreeObj;
    delete[] graphObj;
    delete[] tableObj;
    delete[] myscores;
    delete[] degreeObjs;
    delete[] graphObjs;
    delete[] tableObjs;

    PutRNGstate();

    return Routput;
  }
}
