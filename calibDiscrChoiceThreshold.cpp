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
const int NC=1; // 1 characteristicum (info needed for unpacking of data[] later on)
const int NS=0; // 0 sociodemographic vars (info needed for unpacking of data[] later on)

// Selector for alternatives
int delta(int k1, int k2){return (k1==k2) ? 1 : 0;}

// some nonlinear functions
double cutofflin(double x, double a){return (fabs(x)>a) ? x/fabs(x)*(fabs(x)-a) : 0;}
double powfun(double x, double p){return (fabs(x)>1e-6) ? x/fabs(x)*pow(fabs(x),p) : 0;}



//######################################
/**
 deterministic utility function V for alternative k for a given person
 - NK=number of alternatives
 - Cdata=characteristica in form C_00, C_01,...C_jk ... C_{NC-1,NK-1}
   e.g., C[k]=Traveltime(k), C[1*NK+k]=Cost(k)
 - Sdata=Socioeconomic vars, 
   e.g. S[0]=gender (0=male,1=female), S[1]=age, ...
 - beta = parameters
*/
//######################################


///*

// func1: tanh with horizontal plateau of arbitrary width

const int NK=2; 
const int Mparam=3;
double Vfunc(int k, int nk, 
	     const double CdataPers[], const double SdataPers[], 
	     const double beta[]){
  double dT=CdataPers[1]-CdataPers[0];
  return delta(k,1)*(beta[0]+beta[1]*(dT-beta[2]*tanh(dT/beta[2])));
}
const double betainit0[]={0,0,1};
//*/



/*

// func1tri: !!!

const int NK=3; 
const int Mparam=3;
double Vfunc(int k, int nk, 
	     const double CdataPers[], const double SdataPers[], 
	     const double beta[]){
  double dT1=CdataPers[1]-CdataPers[0];
  double dT2=CdataPers[2]-CdataPers[0];
  return delta(k,1)*(beta[0]+beta[1]*(dT1-beta[2]*tanh(dT1/beta[2])))
  + delta(k,2)*(beta[0]+beta[1]*(dT2-beta[2]*tanh(dT2/beta[2])));
}
const double betainit0[]={0,0,1};
*/





// func 2: tanh with plateau of arbitrary angle and width

/*
const int NK=2; 
const int Mparam=4;
double Vfunc(int k, int nk, 
	     const double CdataPers[], const double SdataPers[], 
	     const double beta[]){
  double dT=CdataPers[1]-CdataPers[0];
  return delta(k,1)*(beta[0]+beta[1]*(dT-beta[2]*tanh(dT/beta[3])));
}
const double betainit0[]={0,0,0,1};
*/

// func 3: purely linear function

/*
const int NK=2; 
const int Mparam=2;
double Vfunc(int k, int nk, 
	     const double CdataPers[], const double SdataPers[], 
	     const double beta[]){
  double dT=CdataPers[1]-CdataPers[0];
  return delta(k,1)*(beta[0]+beta[1]*dT);
}
const double betainit0[]={0,0.1};
*/


// func 4: stepwise-linear with plateau of arbitrary angle, fixed width

/*
const int NK=2; 
const int Mparam=3;
const double plateauwidth=7;
double Vfunc(int k, int nk, 
	     const double CdataPers[], const double SdataPers[], 
	     const double beta[]){
  double dT=CdataPers[1]-CdataPers[0];
  return delta(k,1)*(beta[0]+beta[1]*dT+beta[2]*cutofflin(dT, plateauwidth));
}
const double betainit0[]={0,0,1};
*/


// func 5: stepwise-linear with plateau of arbitrary angle and arbitrary width

/*
const int NK=2; 
const int Mparam=4;
double Vfunc(int k, int nk, 
	     const double CdataPers[], const double SdataPers[], 
	     const double beta[]){
  double dT=CdataPers[1]-CdataPers[0];
  return delta(k,1)*(beta[0]+beta[1]*dT+beta[2]*cutofflin(dT, beta[3]));
}
const double betainit0[]={0,0,0,1};
*/



// func 6: linear + power-function with variable power

/*
const int NK=2; 
const int Mparam=4;
double Vfunc(int k, int nk, 
	     const double CdataPers[], const double SdataPers[], 
	     const double beta[]){
  double dT=CdataPers[1]-CdataPers[0];
  return delta(k,1)*(beta[0]+beta[1]*dT+beta[2]*powfun(dT,beta[3]));
}
const double betainit0[]={0,0,0,0.5};
*/


// func 7: linear + power-function with const power

/*
const int NK=2; 
const int Mparam=3;
const double power=0.5;
double Vfunc(int k, int nk, 
	     const double CdataPers[], const double SdataPers[], 
	     const double beta[]){
  double dT=CdataPers[1]-CdataPers[0];
  return delta(k,1)*(beta[0]+beta[1]*dT+beta[2]*powfun(dT,power));
}
const double betainit0[]={0,0,0};
*/


// func 8: only power-function with variable power

/*
const int NK=2; 
const int Mparam=3;
double Vfunc(int k, int nk, 
	     const double CdataPers[], const double SdataPers[], 
	     const double beta[]){
  double dT=CdataPers[1]-CdataPers[0];
  return delta(k,1)*(beta[0]+beta[1]*powfun(dT,beta[2]));
}
const double betainit0[]={0,0,0.5};
*/




//######################################
// Binomial logit
//######################################



// estimated endogeneous variable hatz= array of sqrt(y_kn*ln(1/P_kn)), 
// so -ln L in "levmar":
// -ln L= sum_{kn} y_kn ln P_kn=sum_kn(hatz-0)^2=min


void calcLSEsummandsLogit2(double *beta, 
			     //double (*V1)(double dT, double *beta),
			     double *hatz,
			     int Mparam, int n_est, double *data
			     ){

  // extract exogeneous and ctrl. data


  // n_est=nlines*K=number of estimated prob =size of hatz
  // nlines=number of data lines/persons/decisions
  // NK=number of alternatives
  // NC=number of characteristica
  // NS=number of sociodemographic vars

  int nlines=n_est/NK;
  double Cdata[nlines][NC*NK];  // Cdata[n]={C_00n, C_01n,...C_jkn ... C_{NC-1,NK-1,n}
  double Sdata[nlines][NS];  // Sdata[n]={S_0n, ... S_{NS-1,n}
  double ydata[nlines][NK];  // ydata[n]={y_00, ..., y_{NK-1,n}
  for(int n=0; n<nlines; n++){
    for (int j=0; j<NC; j++){
      for (int k=0; k<NK; k++){
	Cdata[n][j*NK+k]=data[n*(NC*NK+NS+NK)+j*NK+k];
      }
    }
    for (int s=0; s<NS; s++){
      Sdata[n][s]=data[n*(NC*NK+NS+NK)+NC*NK+s];
    }
    for (int k=0; k<NK; k++){
      ydata[n][k]=data[n*(NC*NK+NS+NK)+NC*NK+NS+k];
    }
  }


  // do the calculation

  //cout <<"in calcLSEsummandsLogit2:"<<endl;

  for(int n=0; n<nlines; n++){
    double expV[NK];
    for(int k=0; k<NK; k++){
      expV[k]=exp(Vfunc(k, NK, Cdata[n], Sdata[n], beta));
    }

    double denom=0;
    for(int k=0; k<NK; k++){denom+=expV[k];}

    double prob[NK];
    for(int k=0; k<NK; k++){prob[k]=expV[k]/denom;}
   
    for(int k=0; k<NK; k++){hatz[NK*n+k]=sqrt(-ydata[n][k]*log(prob[k]));}
    
    if(false){
      if(n==0){ cout <<endl;}
      cout <<"n="<<n<<" V1="<<log(expV[1])<<" prob1="<<prob[1]
    	 <<" hatz0="<<hatz[NK*n]<<" hatz1="<<hatz[NK*n+1]<<endl;
    }
  }

}



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
  if(true){
    for (int i=0; i<n_est; i++){cout <<"i="<<i<<" hatz[i]="<<hatz[i]<<endl;}
  }

  // calculate SSE of (hatz[i]-0) and return it

  double obj=0;
  for(int i=0; i<n_est; i++){
    obj+=hatz[i]*hatz[i];
  }
  return obj;
}

//(11.10.12) !!! Bis hier!!!!!

//#############################################################################
//#############################################################################
int main(int argc, char* argv[]){
//#############################################################################
//#############################################################################

  // ####################################################
  // input: File names (for test: without input file)
  // ####################################################

  if (argc!=1){ // 1: 0 args, 2: 1 arg etc
    cerr <<"\nCalling sequence: calibDiscrChoice <datafile>  <calcObjLandscape>\n"; 
    cerr <<" calcObjLandscape=1 if objective function landscape calulated (runtime!), 0 if not"<<endl;
    cerr  <<"Example:\n";
    cerr <<" calibDiscrChoice timePlateau.dat 0\n";
    exit (-1);
  }

  char datafileName[256];
  char outName[256];



  // ####################################################
  // Input: Exogeneous and endogeneous measurement data
  // ####################################################

  InOut inout;

  // test: without input file

  // NC=1 number of characteristica
  // NS=0 number of socioecon vars
  // NK=2 number of altern.

  int nlines=18;
  int ncols=NC*NK+NS+NK;
  int ndata=nlines*ncols;
  int n_est=nlines*NK;

  // from gnuplot file

  double T1data[NDATA];
  double T2data[NDATA];
  double y1data[NDATA];
  double y2data[NDATA];

  T1data[0]=1;   y1data[0]=0.50;  // T1data[0]=1 markedly different from =0 for nl funcs!
  T1data[1]=5;   y1data[1]=0.50; // reference: 0.50
  T1data[2]=10;  y1data[2]=0.42; // reference: 0.42
  T1data[3]=15;  y1data[3]=0.10;
  T1data[4]=20;  y1data[4]=0.06;
  T1data[5]=25;  y1data[5]=0.04;
  T1data[6]=30;  y1data[6]=0.03;
  T1data[7]=40;  y1data[7]=0.01;
  T1data[8]=50;  y1data[8]=0.00;
  for (int n=nlines/2; n<nlines; n++){
    T1data[n]=-T1data[n-nlines/2];
    y1data[n]=1-y1data[n-nlines/2];
  }

  if(NK==3){ // trinomial
    double dT2=10;
    double f2f0=1.1; //f2/f0 if T2-T0=dT2
    for (int i=0; i<nlines/2-1; i++){
      T2data[2*i+0]=dT2; y2data[2*i+0]=(1-y1data[2*i+0])/f2f0;
      T2data[2*i+1]=-dT2; y2data[2*i+1]=(1-y1data[2*i+1])*f2f0;
      cout <<"2*i="<<2*i<<" y1data[2*i]="<<y1data[2*i]<<" y2data[2*i]="<<y2data[2*i]<<endl;
    }
    T2data[nlines/2-1]=0; y2data[nlines/2-1]=(1-y1data[nlines/2-1]);
    for (int i=0; i<nlines/2-1; i++){
      double norm=1+y2data[i];
      y1data[i]/=norm;
      y2data[i]/=norm;
    }
    for (int n=nlines/2; n<nlines; n++){
      T1data[n]=-T1data[n-nlines/2];
      T2data[n]=-T2data[n-nlines/2];
      y1data[n]=pow(1-y1data[n-nlines/2]-y2data[n-nlines/2],2)/(y1data[n-nlines/2]+1e-6);
      y2data[n]=pow(1-y1data[n-nlines/2]-y2data[n-nlines/2],2)/y2data[n-nlines/2];
      double norm=1+y1data[n]+y2data[n];
      y1data[n]/=norm;
      y2data[n]/=norm;
    }
  }
  for (int n=0; n<nlines; n++){
    cout <<"n="<<n<<" y1data[n]="<<y1data[n]<<" y2data[n]="<<y2data[n]<<endl;
  }


  double data[ndata];

  for(int n=0; n<nlines; n++){

    data[n*(NC*NK+NS+NK)+0]=0;T1data[n];
    data[n*(NC*NK+NS+NK)+1]=T1data[n];
    data[n*(NC*NK+NS+NK)+NC*NK+NS+0]=1-y1data[n];
    data[n*(NC*NK+NS+NK)+NC*NK+NS+1]=y1data[n];
    if(NK==3){
      data[n*(NC*NK+NS+NK)+2]=T2data[n];
      data[n*(NC*NK+NS+NK)+NC*NK+NS+0]=1-y1data[n]-y2data[n];
      data[n*(NC*NK+NS+NK)+NC*NK+NS+1]=y1data[n];
      data[n*(NC*NK+NS+NK)+NC*NK+NS+2]=y2data[n];
    }

  }



  if(true){
    cout <<endl
	 <<"nlines="<<nlines<<" ncols="<<ncols<<" n_est="<<n_est<<" ndata="<<ndata
	 <<endl<<endl;
    for (int n=0; n<nlines; n++){
      for (int icol=0; icol<ncols; icol++){
        cout <<"n="<<n<<" icol="<<icol<<" data[n*ncols+icol]="<<data[n*ncols+icol]<<endl;
      }
      cout <<endl;
    }
  }
  
  //cout <<"developm exit"<<endl; exit(0);

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
  // initialize parameter vector and define box boundaries
  // ####################################################

  double beta[Mparam]; // parameter vector
  for (int m=0; m<Mparam; m++){
    beta[m]=betainit0[m];
  }

  // ####################################################
  // Test objective function for initial parameter settings
  // ####################################################


  cout << "initial estimate of parameter vector:\n";
  for (int m=0; m<Mparam; m++) cout <<"beta["<<m<<"]="<<beta[m]<<"\t";
  cout<<"-logL="<<objFun(calcLSEsummandsLogit2,beta,Mparam,nlines*NK,data)<<endl;
  //cout <<"Development exiting"<<endl; exit(0);

   // ####################################################
   // invoke the optimization function
   // ####################################################

  double  zdata[NDATA];;

    ret= dlevmar_dif(calcLSEsummandsLogit2, beta, zdata, Mparam, nlines*NK, 
  		     n_iter_max, opts, info, NULL, covar, data);

   // check if matrix is not singular 
   // (unfortunately, DOS by output diagnostic info[6], so need to determine
   // indirectly by investigating the correlation matrix


   bool isSingular=false;


   for (int im=1; im<Mparam; im++){
    for (int im2=0; im2<im; im2++){
      double corr=covar[Mparam*im+im2]/sqrt(covar[Mparam*im+im]*covar[Mparam*im2+im2]);
      cout <<"im="<<im<<" im2="<<im2<<" corr="<<corr<<endl;
      if( !(fabs(corr)<=1)){isSingular=true;}
    }
   }


 
  double betaFinal[Mparam];
  for (int k=0; k<Mparam; k++){betaFinal[k]=beta[k];}

  double SSEinit=info[0];
  double SSEfinal=info[1];


  // ####################################################
  // print result
  // ####################################################

  char out_string[5000];
  sprintf(out_string,"%s#=============================================\n",out_string);
  cout <<"\nResults of fitting:\n";
  cout <<"==================================================\n";
 
  //sprintf(out_string,"%s# info[4]=%f",out_string, info[4]);


  char reasonstr[256];
  sprintf(reasonstr,"");
  if(info[6]==1){sprintf(reasonstr,"small gradient");}
    else if(info[6]==2){sprintf(reasonstr,"small param change");}
    else if(info[6]==3){sprintf(reasonstr,"itmax reached");}
    else if(info[6]==4){sprintf(reasonstr,"singular matrix");}
    else if(info[6]==7){sprintf(reasonstr,"NaN or Inf");}
    else{sprintf(reasonstr,"other reason: %g",info[6]);}

  sprintf(out_string,"%s\n#Levenberg-Marquardt returned in %g iter;  reason: %s\n",
     out_string, info[5], reasonstr);
  sprintf(out_string,"%s# number of function evaluations: %g\n", out_string,info[7]);
  sprintf(out_string,"%s# number of Jakobian evaluations: %g\n\n", out_string,info[8]);
  sprintf(out_string,"%s# resulting SSE: %g\n", out_string,SSEfinal);
  sprintf(out_string,"%s# initial  SSE: %g\n", out_string,info[0]);

  char fit_string[4096];  
  sprintf(fit_string,"\n##Best fit parameters (v0=2*max(vdata) if choice_model=1):\n");
  double stddev[Mparam];



  for (int im=0; im<Mparam; im++){
    stddev[im]=sqrt(covar[Mparam*im+im]);
    sprintf(fit_string,"%s#beta[%i]=\t%2.5f; sig_%i=%2.5f\n",
	    fit_string,im,betaFinal[im],im, stddev[im]);
  }


  char corr_string[1000];  
  sprintf(corr_string,"\n##Parameter corr Matrix:\n");
  for (int im=0; im<Mparam; im++){
    sprintf(corr_string, "%s\n#", corr_string);
    for (int im2=0; im2<Mparam; im2++){
      sprintf(corr_string, "%sr%i%i=%1.5f;\t", 
	     corr_string, im, im2, covar[Mparam*im+im2]/(stddev[im]*stddev[im2]));
    }
   
  } 

  sprintf(out_string,"%s%s%s", out_string, fit_string, corr_string);



  //############################################
  // Write result to file
  //############################################

  cout <<out_string<<endl;

  cout <<"\nwriting results to "<<outName<<endl;

  /*
  inout.write_array(outName, ndata, tdata, sdata, vdata, vldata, adata,
		      hataFinal, out_string,
		      "\n\n#time\t\tsdata\t","vdata\t","vldata\t","adata\t", "hata");
  */
  exit(0);
}

