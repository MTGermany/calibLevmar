////////////////////////////////////////////////////////////////////////////////////
//  Copyright (C) 2008  Manolis Lourakis (lourakis at ics forth gr)
//  Institute of Computer Science, Foundation for Research & Technology - Hellas
//  Heraklion, Crete, Greece.
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
// c++ 
#include <iostream>
using namespace std;

#include "levmar.h"

#include "InOut.h"
#include "Math.h"
#include "Statistics.h"
#include "RandomUtils.h"
#include "general.h"

#ifndef LM_DBL_PREC
#error Example program assumes that levmar has been compiled with double precision, see LM_DBL_PREC!
#endif


//######################################
// user defined global settings
//######################################

const int NDATA=10000; // n <= NDATA measurements 

// Selector for alternatives
int delta(int i, int j){return (i==j) ? 1 : 0;}

// some nonlinear functions
double cutofflin(double x, double a){return (fabs(x)>a) ? x/fabs(x)*(fabs(x)-a) : 0;}
double powfun(double x, double p){return (fabs(x)>1e-6) ? x/fabs(x)*pow(fabs(x),p) : 0;}

// function for deterministic utilities
#include "Vfunc.cpp"



//######################################
// i.i.d. Progit
//######################################

// estimated endogeneous variable hatz= array of sqrt(y_kn*ln(1/P_kn)), 
// so -ln L in "levmar":
// -ln L= sum_{kn} y_kn ln P_kn=sum_kn(hatz-0)^2=min


void calcLSEsummandsProbit(double *beta, 
			   //double (*Vfunc)(int i, const double CdataPers[], 
			  //const double SdataPers[],const double beta[]),
			     double *hatz,
			     int Mparam, int n_est, double *data
			     ){


  // extract exogeneous and endogeneous input data from data container


  // n_est=nlines*I=number of estimated prob =size of hatz
  // nlines=n_est/NI number of data lines/persons/decisions
  // NI=number of alternatives (external const)
  // NC=number of characteristica (external const)
  // NS=number of sociodemographic vars (external const)
  // func: preliminary external function
  //==========================================================================
  // data: data container:
  //      * elements 0.. n_input-1: data exog. and endog. var 
  //        (n_input=nlines*(NC*NI+NS+NI)
  //        - one data line= NC*NI+NS+NI elements:
  //          - NC*NI generic vars C_{jni} (not paramlin. factors!)
  //          - NS socioeconomic or external vars (like weather) S_{jn}
  //          - NI observed endogeneous vars y_{ni}
  //      * elements n_input .. n_input+n_endog-1 (n_endog=nlines*NI):
  //        modelled endogeneous vars haty_{ni}=y_nP_{ni}]
  //      * elements n_input+n_endog,n_input+n_endog+1,...: Ctrl vars  
  //==========================================================================

  int nlines=n_est/NI;
  int n_input=nlines*(NC*NI+NS+NI);
  int n_endog=nlines*NI;

  //order is other way round as in main !! necessary sinc rad by columns, worked by lines

  double Cdata[nlines][NC*NI];  // Cdata[n]={C_00n, C_01n,...C_jkn ... C_{NC-1,NI-1,n}
  double Sdata[nlines][NS];  // Sdata[n]={S_0n, ... S_{NS-1,n}
  double ydata[nlines][NI];  // ydata[n]={y_00, ..., y_{NI-1,n}

  for(int n=0; n<nlines; n++){
    for (int j=0; j<NC; j++){
      for (int i=0; i<NI; i++){
	Cdata[n][j*NI+i]=data[n*(NC*NI+NS+NI)+j*NI+i];
      }
    }
    for (int s=0; s<NS; s++){
      Sdata[n][s]=data[n*(NC*NI+NS+NI)+NC*NI+s];
    }
    for (int i=0; i<NI; i++){
      ydata[n][i]=data[n*(NC*NI+NS+NI)+NC*NI+NS+i];
    }

  }

  // extract ctrl data (preserve NI*n_lines elements of data[] for the P_ni)

  bool writeEndogVars=(int)(data[n_input+n_endog]+0.5);


  //########################
  // do the calculation (Probit)
  //########################

  const int NINT=80;
  const double zmin=-4;
  const double zmax=4;
  Math myMath;

  double dz=(zmax-zmin)/(NINT-2);
  double dpNormal[NINT];
  //double pref=1./sqrt(2*PI);

  dpNormal[0]=myMath.norm(zmin);
  double pNormalLast=dpNormal[0];
  for (int iz=1; iz<NINT-1; iz++){
    double z=zmin+iz*dz;
    dpNormal[iz]=myMath.norm(z)-pNormalLast;
    pNormalLast+=dpNormal[iz];
    //cout <<"z="<<z<<" dpNormal[iz]="<<dpNormal[iz]<<" pNormalLast="<<pNormalLast<<endl;
  }
  //cout <<"myMath.norm(zmax)="<<myMath.norm(zmax)<<endl;
  dpNormal[NINT-1]=1-myMath.norm(zmax);
  
 

  double prob[nlines][NI];

  for(int n=0; n<nlines; n++){

    double vdiff[NI][NI];
    for (int i=0; i<NI; i++){
      for (int j=0; j<NI; j++){
	vdiff[i][j]=Vfunc(i,Cdata[n],Sdata[n], beta) -Vfunc(j,Cdata[n],Sdata[n], beta);
	//cout <<"n="<<n<<" i="<<i<<" j="<<j<<" vdiff[i][j]="<<vdiff[i][j]<<endl;
      }
    }


    for (int i=0; i<NI-1; i++){

      // do integral

      double p=0;
      for (int iz=0; iz<NINT; iz++){
	double z=zmin+(iz-0.5)*dz;
        double dp=dpNormal[iz];
	for (int j=0; j<NI; j++){
	  if(j!=i){
	    //cout <<"i="<<i<<" iz="<<iz<<" z="<<z<<" j="<<j<<" vdiff="<<vdiff[i][j]<<endl;
	    dp*=myMath.norm(z+vdiff[i][j]);
	  }
	}
	p+=dp;
	//cout <<"dpNormal[iz]="<<dpNormal[iz]<<" dp="<<dp<<" p="<<p<<endl;
      }
      prob[n][i] = p;
    }
    prob[n][NI-1]=1;
    for (int i=0; i<NI-1; i++){prob[n][NI-1]-=prob[n][i];}

    for(int i=0; i<NI; i++){
      prob[n][i]=max(1e-10, min(1.,prob[n][i]));
      hatz[NI*n+i]=sqrt(-ydata[n][i]*log(prob[n][i]));
    }

    //test
    if(false){
      for (int i=0; i<NI; i++){
	cout <<"n="<<n<<" i="<<i<<" prob[n][i]="<<prob[n][i]<<" hatz[NI*n+i]="<<hatz[NI*n+i]<<endl;
      }
    }


  }

  // write endogeneous vars in data[] container (only after finishing estimation)
  if(writeEndogVars){
    for(int n=0; n<nlines; n++){
      double ysum=0; for(int i=0; i<NI; i++){ysum+=ydata[n][i];}
      for(int i=0; i<NI; i++){
	data[n_input+n*NI+i]=ysum*prob[n][i];
      }
    }
  }

} //calcLSEsummandsProbit







//######################################
// Logit
//######################################

// estimated endogeneous variable hatz= array of sqrt(y_kn*ln(1/P_kn)), 
// so -ln L in "levmar":
// -ln L= sum_{kn} y_kn ln P_kn=sum_kn(hatz-0)^2=min

void calcLSEsummandsLogit(double *beta, 
			   //double (*Vfunc)(int i, const double CdataPers[], 
			  //const double SdataPers[],const double beta[]),
			     double *hatz,
			     int Mparam, int n_est, double *data
			     ){


  // extract exogeneous and endogeneous input data from data container


  // n_est=nlines*I=number of estimated prob =size of hatz
  // nlines=n_est/NI number of data lines/persons/decisions
  // NI=number of alternatives (external const)
  // NC=number of characteristica (external const)
  // NS=number of sociodemographic vars (external const)
  // func: preliminary external function
  //==========================================================================
  // data: data container:
  //      * elements 0.. n_input-1: data exog. and endog. var 
  //        (n_input=nlines*(NC*NI+NS+NI)
  //        - one data line= NC*NI+NS+NI elements:
  //          - NC*NI generic vars C_{jni} (not paramlin. factors!)
  //          - NS socioeconomic or external vars (like weather) S_{jn}
  //          - NI observed endogeneous vars y_{ni}
  //      * elements n_input .. n_input+n_endog-1 (n_endog=nlines*NI):
  //        modelled endogeneous vars haty_{ni}=y_nP_{ni}]
  //      * elements n_input+n_endog,n_input+n_endog+1,...: Ctrl vars  
  //==========================================================================

  int nlines=n_est/NI;
  int n_input=nlines*(NC*NI+NS+NI);
  int n_endog=nlines*NI;

  //order is other way round as in main !! necessary sinc rad by columns, worked by lines

  double Cdata[nlines][NC*NI];  // Cdata[n]={C_00n, C_01n,...C_jkn ... C_{NC-1,NI-1,n}
  double Sdata[nlines][NS];  // Sdata[n]={S_0n, ... S_{NS-1,n}
  double ydata[nlines][NI];  // ydata[n]={y_00, ..., y_{NI-1,n}

  for(int n=0; n<nlines; n++){
    for (int j=0; j<NC; j++){
      for (int i=0; i<NI; i++){
	Cdata[n][j*NI+i]=data[n*(NC*NI+NS+NI)+j*NI+i];
      }
    }
    for (int s=0; s<NS; s++){
      Sdata[n][s]=data[n*(NC*NI+NS+NI)+NC*NI+s];
    }
    for (int i=0; i<NI; i++){
      ydata[n][i]=data[n*(NC*NI+NS+NI)+NC*NI+NS+i];
    }

  }

  // extract ctrl data (preserve NI*n_lines elements of data[] for the P_ni)

  bool writeEndogVars=(int)(data[n_input+n_endog]+0.5);

  // do the calculation

  double prob[nlines][NI];

  //cout <<"in calcLSEsummandsLogit:"<<endl;

  for(int n=0; n<nlines; n++){
    double expV[NI];
    for(int i=0; i<NI; i++){
      expV[i]=exp(Vfunc(i, Cdata[n], Sdata[n], beta));
    }

    double denom=0;
    for(int i=0; i<NI; i++){denom+=expV[i];}
    for(int i=0; i<NI; i++){prob[n][i]=expV[i]/denom;}
    for(int i=0; i<NI; i++){hatz[NI*n+i]=sqrt(-ydata[n][i]*log(prob[n][i]));}
  }

  // write endogeneous vars in data[] container (only after finishing estimation)
  if(writeEndogVars){
    for(int n=0; n<nlines; n++){
      double ysum=0; for(int i=0; i<NI; i++){ysum+=ydata[n][i];}
      for(int i=0; i<NI; i++){
	data[n_input+n*NI+i]=ysum*prob[n][i];
      }
    }
  }

}//calcLSEsummandsLogit



//################################################################################
// Calculate objective function 
// (the obove functions calculate the summands of the LSE minimization
//################################################################################


double objFun(
	      void (*calcLSEsummands)(double *beta, 
				      // double (*V1)(double dT, double *beta), 
				      double *hatz,
				      int Mparam, int n_est, double *data),
	      double *beta, int Mparam, int n_est, double *data
	      ){

  // array of compares=0
  double hatz[NDATA];
  // calculate array of hatz
  calcLSEsummands(beta,hatz, Mparam, n_est, data);
  if(false){ // see also debug output in calcLSEsummands
    cout <<endl;
    for (int i=0; i<n_est; i++){
      cout <<"i="<<i<<" hatz[i]="<<hatz[i]<<endl;} 

  }

  // calculate SSE of (hatz[i]-0) and return it

  double obj=0;
  for(int i=0; i<n_est; i++){
    obj+=hatz[i]*hatz[i];
  }
  return obj;
}


//#############################################################################
//#############################################################################
int main(int argc, char* argv[]){
//#############################################################################
//#############################################################################

  // ####################################################
  // input: File names (for test: without input file)
  // ####################################################

  if ( (argc<3) || (argc>4)){ // 1: 0 args, 2: 1 arg etc
    cerr <<"\nCalling sequence: calibDiscrChoice <projName>  <calcObjLandscape> [<model>]\n"; 
    cerr <<" calcObjLandscape=1 if objective function landscape calulated (runtime!), 0 if not"<<endl;
    cerr<<" Model: 0=Logit, 1=i.i.d. Probit, 2=with inclVals for nLogit\n";
    cerr  <<"Example:\n";
    cerr <<" calibDiscrChoice statedChoiceWS1213 0  [Logit, default]\n";
    cerr <<" calibDiscrChoice statedChoiceWS1213 0 0  [Logit]\n";
    cerr <<" calibDiscrChoice statedChoiceWS1213 0 1  [Probit]\n";
    cerr <<" calibDiscrChoice statedChoiceWS1213 0 2  [outp with inclVals]\n";
    exit (-1);
  }

  char projName[256];
  char datafileName[256];
  char ICfileName[256];
  char outName[256];

  sprintf(projName, "%s", argv[1]);
  bool calcObjLandscape=(atoi(argv[2])==1);
  sprintf(datafileName, "%s.inputData",projName);
  sprintf(ICfileName, "%s.IC",projName);
  sprintf(outName, "%s.outData",projName);

  // ####################################################
  // Input: Type of model (0=Logit, 1=iid Probit, 2=Logit w/ incl for NL)
  // #################################################### 

  int modelType=0;  // Logit as default
  if(argc==4){modelType=atoi(argv[3]);}

  void (*calcLSEsummands)(double *beta, 
			double *hatz,
			int Mparam, int n_est, double *data)
  =(modelType!=1)? calcLSEsummandsLogit : calcLSEsummandsProbit;


  // ####################################################
  // Input: Exogeneous and endogeneous measurement data
  // ####################################################

  InOut inout;
  int nlines;
  double Cdata[NC*NI][NDATA]; // transposed with respect to subroutines!
  double Sdata[NS][NDATA]; // transposed with respect to subroutines!
  double ydata[NI][NDATA]; // transposed with respect to subroutines!
  double ysum[NDATA];

  // characteristica

  for (int i=0; i<NI; i++){
    for (int nc=0; nc<NC; nc++){
      inout.get_col(datafileName,1+nc*NI+i,nlines,Cdata[nc*NI+i]);
    }
  }

  //sociodemographic vars

  for (int ns=0; ns<NS; ns++){
    inout.get_col(datafileName,1+NC*NI+ns,nlines,Sdata[ns]);
  }

  // choice numbers

  for (int i=0; i<NI; i++){
    inout.get_col(datafileName,1+NC*NI+NS+i,nlines,ydata[i]); // y_{ni}
  }

  // get sum of choices made per line

  for (int n=0; n<nlines; n++){
    ysum[n]=0;
    for (int i=0; i<NI; i++){
      ysum[n]+=ydata[i][n];
    }
  }


  // NC=number of characteristica
  // NS=number of socioecon vars
  // NI=number of altern.

  int n_cols=NC*NI+NS+NI;
  int n_input=nlines*n_cols;
  int n_est=nlines*NI;
  int n_endog=n_est;
  int n_ctrl=1;
  int n_data=n_input+n_endog+n_ctrl; 


  for (int n=0; n<nlines; n++){
    cout <<"n="<<n;
    for (int nc=0; nc<NC; nc++){
      for (int i=0; i<NI; i++){
        cout <<" C_{"<<nc<<","<<n<<","<<i<<"}="<<Cdata[nc*NI+i][n]<<"\t";
      }
    }
    cout <<endl;
  }
  cout <<endl;

  for (int n=0; n<nlines; n++){
    cout <<"n="<<n;
    for (int ns=0; ns<NS; ns++){
      cout <<" S_{"<<ns<<","<<n<<"}="<<Sdata[ns][n]<<"\t";
    }
    cout <<endl;
  }
  cout <<endl;

  for (int n=0; n<nlines; n++){
    cout <<"n="<<n;
    for (int i=0; i<NI; i++){
      cout <<" y_{"<<n<<","<<i<<"}="<<ydata[i][n]<<"\t";
    }
    cout <<endl;
  }


  // define and fill data container
  //==========================================================================
  // data: data container:
  //      * elements 0.. n_input-1: data exog. and endog. var 
  //        (n_input=nlines*(NC*NI+NS+NI)
  //        - one data line= NC*NI+NS+NI elements:
  //          - NC*NI generic vars C_{jni} (not paramlin. factors!)
  //          - NS socioeconomic or external vars (like weather) S_{jn}
  //          - NI observed endogeneous vars y_{ni}
  //      * elements n_input .. n_input+n_endog-1 (n_endog=nlines*NI):
  //        modelled endogeneous vars haty_{ni}=y_nP_{ni}]
  //      * elements n_input+n_endog,n_input+n_endog+1,...: output ctrl vars  
  //==========================================================================



  double data[n_data];

  for(int n=0; n<nlines; n++){
    for(int i=0; i<NI; i++){
      for (int nc=0; nc<NC; nc++){
        data[n*(NC*NI+NS+NI)+nc*NI+i]=Cdata[nc*NI+i][n];
      }
    }
  }


  for(int n=0; n<nlines; n++){
    for (int ns=0; ns<NS; ns++){
      data[n*(NC*NI+NS+NI)+NC*NI+ns]=Sdata[ns][n];
    }
  }

  for(int n=0; n<nlines; n++){
    for(int i=0; i<NI; i++){
      data[n*(NC*NI+NS+NI)+NC*NI+NS+i]=ydata[i][n];
    }
  }

  data[n_input+n_endog]=0; // no writing of endog. vars

  //==========================================================================
 
  if(false){
    cout <<endl
	 <<"nlines="<<nlines<<" n_cols="<<n_cols
	 <<" n_est="<<n_est<<" n_data="<<n_data
	 <<endl<<endl;
    for (int n=0; n<nlines; n++){
      for (int icol=0; icol<n_cols; icol++){
        cout <<"n="<<n<<" icol="<<icol<<" data[n*n_cols+icol]="<<data[n*n_cols+icol]<<endl;
      }
      cout <<endl;
    }
  }


  // ####################################################
  // Input: Ctrl data
  // ####################################################


  double opts[LM_OPTS_SZ]; // control parameters for dlevmar_[dif|der] (input)
  double info[LM_INFO_SZ]; // information on the calculation of results (output) 
  const int n_iter_max=1000; // maximum number of iterations

  // passing to levmar NULL instead of opts reverts to defaults 
  opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
  opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used 


  double covar[Mparam*Mparam];
  int ret; // return signal of the dlevmar functions



  // ####################################################
  // Input: Parameter initial conditions
  // ####################################################

  double beta[Mparam]; // parameter vector
  for (int m=0; m<Mparam; m++){
    beta[m]=0; // default
  }
  int MparamInFile;

  if(inout.fileExists(ICfileName)){
    inout.get_col(ICfileName, 1, MparamInFile, beta);
    cout <<"overridden default parameter zero IC by file IC"<<endl;
  }
  else{cout<<"no .IC file exists"<<endl;}


  for (int m=0; m<Mparam; m++){cout <<"beta["<<m<<"]="<<beta[m]<<endl;}


  // ####################################################
  // Test objective function for initial parameter settings
  // ####################################################

  //beta[0]=-10; beta[1]=-0.36; beta[2]=-0.32; beta[3]=-8; //!!!
  cout << "initial estimate of parameter vector:\n";
  for (int m=0; m<Mparam; m++) cout <<"beta["<<m<<"]="<<beta[m]<<"\t";
  cout <<endl;
  cout<<"-logL="<<objFun(calcLSEsummands,beta,Mparam,nlines*NI,data)<<endl;

  //cout <<"development exiting"<<endl; exit(0);

   // ####################################################
   // invoke the optimization function
   // ####################################################

   double  zdata[NDATA];

   ret= dlevmar_dif(calcLSEsummands, beta, zdata, Mparam, nlines*NI, 
  		     n_iter_max, opts, info, NULL, covar, data);

   // check if matrix is not singular 
   // (unfortunately, DOS by output diagnostic info[6], so need to determine
   // indirectly by investigating the correlation matrix


   bool nearSingular=false;

   //!!!
   for (int im=1; im<Mparam; im++){
     for (int im2=0; im2<im; im2++){
       double corr=covar[Mparam*im+im2]/sqrt(covar[Mparam*im+im]*covar[Mparam*im2+im2]);
       if( !(fabs(corr)<1-1e-6)){cerr<<"is singular!!"<<endl; nearSingular=true;}
       cout <<"step original:  corr_"<<im<<im2<<"="<<corr<<endl;
     }
   }
 
  double betaFinal[Mparam];
  for (int i=0; i<Mparam; i++){betaFinal[i]=beta[i];}

  double SSEinit=info[0];
  double SSEfinal=info[1];


  // ####################################################
  // collect calibration result and write to out_string,
  // ####################################################

  cout <<"=============================================\n";
  cout <<"\nResults of fitting:\n";
  cout <<"==================================================\n";
 

  char out_string[5000];

  char reasonstr[256];
  sprintf(reasonstr,"");
  if(info[6]==1){sprintf(reasonstr,"small gradient");}
    else if(info[6]==2){sprintf(reasonstr,"small param change");}
    else if(info[6]==3){sprintf(reasonstr,"itmax reached");}
    else if(info[6]==4){sprintf(reasonstr,"singular matrix");}
    else if(info[6]==7){sprintf(reasonstr,"NaN or Inf");}
    else{sprintf(reasonstr,"other reason: %g",info[6]);}

  sprintf(out_string,"\n##Levenberg-Marquardt returned in %g iter;  reason: %s\n",
     info[5], reasonstr);
  sprintf(out_string,"%s## number of function evaluations: %g\n", out_string,info[7]);
  sprintf(out_string,"%s## number of Jakobian evaluations: %g\n#\n", out_string,info[8]);
  sprintf(out_string,"%s## Max Log-Likelihood: %g\n", out_string,SSEfinal);
  sprintf(out_string,"%s## initial Log-Likelihood:: %g\n#\n", out_string,info[0]);

  char fit_string[4096];  
  sprintf(fit_string,"##Best fit parameters:\n");
  double stddev[Mparam];


  // !!correct parameter covariance matrix since levmar knows nothing of ML
  // (calculates residual errors from SSE)
  // V_LSE=2*\sigeps^2 H_S^{-1} = 2*\sigeps^2 V_ML
  // =>  V_ML=V_LSE/(2*\sigeps^2) approx V_LSE *(nSummands-Mparam-1)/(2*Smin)

  // 0.5*nlines/SSEfinal scheint das Beste zu sein! (reverse engineering levmar)
  //double corrFactor=0.5*(NI*nlines-Mparam-1)/SSEfinal;
  double corrFactor=0.5*nlines/SSEfinal; 
  //double corrFactor=0.18*(NI*nlines-Mparam-1)/SSEfinal;
  //double corrFactor=0.5*(nlines-Mparam-1)/SSEfinal;

  // predictor for covar to get order right (for some reason not exactly consistent)

  for (int im2=0; im2<Mparam*Mparam; im2++){
    covar[im2] *=corrFactor;
  }
  for (int im=0; im<Mparam; im++){
    stddev[im]=sqrt(covar[Mparam*im+im]);
  }
 
  cout <<out_string<<endl;

  if(true){
    cout <<"calibr result:"<<endl;
    for (int im=0; im<Mparam; im++){cout<<"beta["<<im<<"]="<<beta[im]<<endl;}
    cout <<"prelim stddevs:"<<endl;
    for (int im=0; im<Mparam; im++){cout<<"stddev["<<im<<"]="<<stddev[im]<<endl;}
  }
  

  // calculate numerical second derivatives of log-Likelihood

  double Lmax=-SSEfinal;
  double hesseL[Mparam][Mparam];
  double dbrel=0.0001; // difference quotient with dbrel*predictor stddev

  // diagonal elements

  for (int i=0; i<Mparam; i++){beta[i]=betaFinal[i];}

  for (int im=0; im<Mparam; im++){
    beta[im] = betaFinal[im]+dbrel*stddev[im];
    double Lplus=-objFun(calcLSEsummands,beta, Mparam, nlines*NI, data);
    beta[im] = betaFinal[im]-dbrel*stddev[im];
    double Lminus=-objFun(calcLSEsummands,beta, Mparam, nlines*NI, data);
    hesseL[im][im]=(Lplus-2*Lmax+Lminus)/(dbrel*dbrel*stddev[im]*stddev[im]);
    beta[im] = betaFinal[im];
    if(im<=1){
      cout <<"im="<<im<<" Lplus="<<Lplus<<" Lminus="<<Lminus<<" Lmax="<<Lmax
	   <<" hesseL[im][im]="<<hesseL[im][im]<<endl;
    }
  }



  // off-diagonal elements

  for (int i=0; i<Mparam; i++){beta[i]=betaFinal[i];}
  for (int im=0; im<Mparam; im++){
    for (int im2=0; im2<im; im2++){
      beta[im] =betaFinal[im]+dbrel*stddev[im];
      beta[im2]=betaFinal[im2]+dbrel*stddev[im2];
      double Lplusplus=-objFun(calcLSEsummands,beta, Mparam, nlines*NI, data);
      beta[im2]=betaFinal[im2]-dbrel*stddev[im2];
      double Lplusminus=-objFun(calcLSEsummands,beta, Mparam, nlines*NI, data);
      beta[im] =betaFinal[im]-dbrel*stddev[im];
      beta[im2]=betaFinal[im2]+dbrel*stddev[im2];
      double Lminusplus=-objFun(calcLSEsummands,beta, Mparam, nlines*NI, data);
      beta[im2]=betaFinal[im2]-dbrel*stddev[im2];
      double Lminusminus=-objFun(calcLSEsummands,beta, Mparam, nlines*NI, data);
      hesseL[im][im2]=(Lplusplus-Lplusminus-Lminusplus+Lminusminus)/
	(4*dbrel*dbrel*stddev[im]*stddev[im2]);
      hesseL[im2][im]=hesseL[im][im2];
      beta[im] = betaFinal[im];
      beta[im2] = betaFinal[im2];
    }
  }


 

  // calculate covar matrix as matrix inverse of - Hessematrix

  Math myMath;
  double matrix[myMath.NMATRIX][myMath.NMATRIX];
  double inv[myMath.NMATRIX][myMath.NMATRIX];
  for (int im=0; im<Mparam; im++){
    for (int im2=0; im2<Mparam; im2++){
      matrix[im][im2]=-hesseL[im][im2];
    }
  }

  myMath.calcMatrixInverse(matrix, Mparam, inv, false);

  // test for neg. var 
  // (may happen for num second derivatives; then revert to unprecise orig.)

  bool someNegVariances=false;


  for (int im=0; im<Mparam; im++){
    if(inv[im][im]<=0){
     someNegVariances=true;
     cerr<<" Warning: im="<<im<<" Num Variance inv[im][im]="<<inv[im][im]
	 <<"<=0! Using original unprecise covariance matrix!"<<endl;
    }
  }

  // test code

  if(true){
    for (int im=0; im<Mparam; im++){
      cout <<"im="<<im<<" covar[Mparam*im+im]="<<covar[Mparam*im+im]
  	 <<" inv[im][im]="<<inv[im][im]<<endl;
    }
  }

  if(nearSingular){cerr<<" Warning: nearly singular, corr nearer to +/- 1 than 1e-6"
		     <<endl<<" Using original unprecise covariance matrix!"<<endl;
  }
  //exit(-1);

  //if(someNegVariances||nearSingular){
  if(someNegVariances){
    for (int im=0; im<Mparam; im++){
      for (int im2=0; im2<Mparam; im2++){
	inv[im][im2]=covar[Mparam*im+im2];
      }
    }
  }


 

  // test code

  if(false){
   cout <<endl;
   for (int im=0; im<Mparam; im++){
    for (int im2=0; im2<Mparam; im2++){
      cout <<" hesseL["<<im<<"]["<<im2<<"]="<<hesseL[im][im2];
    }
    cout <<endl;
   }

   cout <<endl;
   for (int im=0; im<Mparam; im++){
    for (int im2=0; im2<Mparam; im2++){
      cout <<" inv["<<im<<"]["<<im2<<"]="<<inv[im][im2];
    }
    cout <<endl;
   }

   cout <<"Comparison original:"<<endl;

   for (int im=0; im<Mparam; im++){
    for (int im2=0; im2<Mparam; im2++){
      cout <<" covar["<<im<<"]["<<im2<<"]="<<covar[Mparam*im+im2];
    }
    cout <<endl;
   }
  }

  //exit(-1);


  // overwrite levmar covar matrix with my calculated -H_L^{-1}

  for (int im=0; im<Mparam; im++){
    for (int im2=0; im2<Mparam; im2++){
      covar[Mparam*im+im2]=inv[im][im2];
    }
    stddev[im]=sqrt(covar[Mparam*im+im]);
  }

  if(someNegVariances){
    cout <<"Negative variances encountered!!"<<endl;
    stddev[0]=0.05; //!!!
    stddev[4]=0.5; //!!!
  }

  sprintf(fit_string,"%s#lnLmax=%f\n",fit_string, -SSEfinal);

  for (int im=0; im<Mparam; im++){
    sprintf(fit_string,"%s#beta%i=\t%2.5f; sig%i=%2.5f\n",
	    fit_string,im,betaFinal[im],im, stddev[im]);
  }




  char corr_string[1000];  
  sprintf(corr_string,"#\n##Parameter corr Matrix:\n#");
  for (int im=0; im<Mparam; im++){
    sprintf(corr_string, "%s\n#", corr_string);
    for (int im2=0; im2<Mparam; im2++){
      sprintf(corr_string, "%sr%i%i=%1.3f;\t", 
	     corr_string, im, im2, covar[Mparam*im+im2]/(stddev[im]*stddev[im2]));
    }
   
  } 



  sprintf(out_string,"%s%s%s", out_string, fit_string, corr_string);

  //######################################################
  // get predicted probabilities for calibrated/estimated param vector
  //######################################################

  data[n_input+n_endog]=1; //activate output of endog. vars in data[]

  // calculate haty and fill in data container 
  objFun(calcLSEsummands,beta,Mparam,nlines*NI,data); 

  // define haty and n_array for writing to file

  double haty[NI][nlines];
  double n_array[nlines];
  for (int n=0; n<nlines; n++){
    n_array[n]=n+1;
    for(int i=0; i<NI; i++){
      haty[i][n]=data[n_input+n*NI+i];
    }
  }





  //################################################################
  // output
  //################################################################



  //################################################################
  // !!! transpose of Cdata and Sdata !!! needed as hack since
  // for Vfunc, I need n as first index of Cdata, Sdata!
  //################################################################

  double CdataT[nlines][NC*NI]; 
  for(int n=0; n<nlines; n++){
    for(int ij=0; ij<NC*NI; ij++){
      CdataT[n][ij]=Cdata[ij][n];
    }
  }


  double SdataT[nlines][NS]; 
  for(int n=0; n<nlines; n++){
    for (int js=0; js<NS; js++){
       SdataT[n][js]=Sdata[js][n];
    }
  }



  //############################################
  // Write estimation and fit results to file
  //############################################

  cout <<out_string<<endl;
  sprintf(out_string,"%s\n#\n#", out_string);
  for(int i=0; i<NI; i++){
    sprintf(out_string,"%s%s%i%s\t", out_string,"y",i,"data");
  }
  for(int i=0; i<NI; i++){
    sprintf(out_string,"%s%s%i\t", out_string,"haty",i);
  }
  sprintf(out_string,"%s%s",out_string,"ysum");
  
  cout <<"\nwriting results to "<<outName<<endl;
  if(NI==4){
    inout.write_array(outName, nlines, n_array, ydata[0],ydata[1],ydata[2], ydata[3],
		      haty[0], haty[1], haty[2], haty[3], ysum,
		      out_string);
		      // "y1data","y2data","y3data","y4data",
		      //"haty1","haty2","haty3","haty4","ysum");
  }
  else if(NI==3){
    inout.write_array(outName, nlines, n_array, ydata[0],ydata[1],ydata[2],
		    haty[0], haty[1], haty[2], ysum, 
		    out_string);
  }
  else if(NI==2){
    inout.write_array(outName, nlines, n_array, ydata[0],ydata[1],
		    haty[0], haty[1], ysum, 
		    out_string);
  }


  //############################################
  // Write inclusion values I_n to file .inclVals 
  // if nest of nLogit model
  //############################################


  if(modelType==2){
    char inclfileName[256];
    sprintf(inclfileName, "%s.inclData",projName);
    char title[256];
    sprintf(title, "#inclusion values ln(sum_m exp(V_nlm) for this nest l");
    sprintf(title,"%s\n#n\tinclVal\tn", title);
    double inclVal[nlines];
    cout <<"in write inclusion vals: title="<<title<<endl;

    for(int n=0; n<nlines; n++){
      double expsum=0;
      for(int i=0; i<NI; i++){
	double expV=exp(Vfunc(i, CdataT[n], SdataT[n], beta));
	cout<<"n="<<n<<" i="<<i<<" CdataT[n][i]="<<CdataT[n][i]
	    <<" Vfunc="<<log(expV)<<" expV="<<expV<<endl;
	expsum += expV;
      }
      inclVal[n]=(expsum>0) ? log(expsum) : 0;
    }

    inout.write_array(inclfileName,nlines,n_array,inclVal,n_array,title); 
    cout <<"wrote "<<inclfileName<<endl;
  }

 //############################################
  // Write log-likelihood around the  maximum to files
  //############################################

  double w_stddev=4.; // write +/- 5 standard deviations

  data[n_input+n_endog]=0; // deactivate writing of endogeneous var in mic*Func(..)

  if(calcObjLandscape==1){
    if(Mparam<2){
      cerr<<"Error: need at least two params for plotting param landscape"<<endl;
      exit(-1);
    }
    int nout=61; // number of grid elements in either direction
    char objFunFilename[1024]; // of form outName_betaj_betam
    char titleString[1024];

    double objFunData[inout.NYMAX][inout.NYMAX];
    double dbeta[Mparam];
    double betamin[Mparam];
    double betamax[Mparam];

    for (int j=0; j<Mparam; j++){
      betamin[j]=betaFinal[j]-w_stddev*stddev[j];
      betamax[j]=betaFinal[j]+w_stddev*stddev[j];
      dbeta[j]=(betamax[j]-betamin[j])/(nout-1);
      cout <<"j="<<j<<" betamin[j]="<<betamin[j]<<" betamax[j]="<<betamax[j]<<endl;
    }


    // make (Mparam-1)*(Mparam-2) data sets and files

    for (int ibeta=0; ibeta<Mparam-1; ibeta++){
      //cout <<"ibeta="<<ibeta<<"betamin[ibeta]="<<betamin[ibeta]<<" betamax[ibeta]="<<betamax[ibeta]<<endl;
      for (int jbeta=ibeta+1; jbeta<Mparam; jbeta++){

      // generate the data set for a given beta combination

        for (int k=0; k<Mparam; k++){// revert non-used dimensions
           beta[k]=betaFinal[k];
	}
	for (int i=0; i<nout; i++){
	  for (int j=0; j<nout; j++){
	    beta[ibeta]=betamin[ibeta]+i*dbeta[ibeta];
	    beta[jbeta]=betamin[jbeta]+j*dbeta[jbeta];
	    objFunData[i][j]=-objFun(calcLSEsummands,beta, Mparam, 
				    nlines*NI, data);
	  }
	}


        // write the file for this combination

        sprintf(objFunFilename,"%s_beta%i_beta%i",projName, ibeta, jbeta);
        sprintf(titleString,"#Log-Likelihood for beta%i and beta%i",ibeta, jbeta);
	sprintf(titleString, "%s\n#Base values: beta0=%.2f, beta1=%.2f", 
		titleString, betaFinal[0], betaFinal[1]);
	for(int m=2; m<Mparam; m++){
	  sprintf(titleString, "%s, beta%i=%.2f", titleString, m, betaFinal[m]);
	}

	sprintf(titleString,"%s\n#Min Obj function: %f", titleString, info[1]);
	sprintf(titleString,"%s\n#beta%i\t\tbeta%i\t\tlogL", titleString, ibeta, jbeta);

	inout.write_array2d(objFunFilename, betamin[ibeta],betamax[ibeta],nout,
			    betamin[jbeta],betamax[jbeta],nout,
			    objFunData,titleString);
      }
    }
  }


  //################################################################
  // Write property sums
  //################################################################

  cout <<endl<<"Property sums (Merkmalssummen):"<<endl;
  cout <<"============================================="<<endl<<endl;

    
  for (int im=0; im<Mparam; im++){
    double betaSelect[Mparam];
    for(int jm=0; jm<Mparam; jm++){
      betaSelect[jm]=delta(jm,im);
    }
    double X_data=0;
    double X_MNL0=0; // predicted property sum for beta=0
    double X_MNL=0;  // predicted property sum for calibrated beta vals

    
    for(int n=0; n<nlines; n++){
       for (int i=0; i<NI; i++){
	X_data += Vfunc(i, CdataT[n], SdataT[n], betaSelect)*ydata[i][n];
	X_MNL0 += 1./NI*Vfunc(i, CdataT[n], SdataT[n], betaSelect)*ysum[n];
	X_MNL  += Vfunc(i, CdataT[n], SdataT[n], betaSelect)*haty[i][n];
	if(false){
	  cout <<" n="<<n<<" i="<<i
	       <<" Vfunc="<<Vfunc(i, CdataT[n], SdataT[n], betaSelect)
	       <<" ydata[i][n]="<<ydata[i][n]
	       <<" haty[i][n]="<<haty[i][n]
	       <<" ysum[n]="<<ysum[n]
	       <<endl;
	}
      }
    }

    cout <<"beta"<<im<<": X_data="<<X_data
	 <<" X_MNL0="<<X_MNL0<<" X_MNL="<<X_MNL<<endl;


    
  }


  

  cout <<endl<<"For input in "<<projName<<".gnu"<<endl;
  cout <<"============================================="<<endl<<endl;
  int Nchoice=0; for(int n=0; n<nlines; n++){Nchoice+=ysum[n];}
  cout <<"lnLinit="<<-SSEinit<<endl;
  cout <<"lnLmax="<<-SSEfinal<<endl;
  cout <<"N="<<Nchoice<<"        #(number of single choices)"<<endl;
  cout <<"Mparam="<<Mparam<<endl;
  for (int im=0; im<Mparam; im++){
    cout<<"beta"<<im<<"="<<betaFinal[im]
	<<"; sig"<<im<<"="<<stddev[im]<<endl;
  }
  cout <<endl;

  
  cout <<endl;
  
  if(nearSingular){
    cerr<<"WARNING: nearly singular (at least one corr coeff nearer to +/- 1 than 1e-6)"
	<<endl << "correlations near 1 may be erroneously given in corr matrix!"
	<<endl;
  }
  exit(0); 
}

