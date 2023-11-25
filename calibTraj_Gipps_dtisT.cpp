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

const int Ndata=10000; // n <= Ndata measurements 


// IDM  (v0,T,s0,a,b)

double accIDM(double s, double v, double dv, const double beta[]){
  double v0=beta[0];
  double  T=maxsmooth(0,beta[1], 0.1);
  double s0=maxsmooth(0,beta[2], 0.1); 
  double  a=maxsmooth(0,beta[3], 0.1); 
  double  b=maxsmooth(0,beta[4], 0.1); 
  //double  b=beta[4];

  double sstar=s0+max(0.,v*T+0.5*v*dv/sqrt(a*b));

  double accIDM=a*(1-pow(v/v0,4) - pow( sstar/s, 2));
  if(isnan(accIDM)){
    cerr<<"s="<<s<<" v="<<v<<" dv="<<dv<<endl;
    cerr<<"v0="<<v0<<endl;
    cerr<<"T="<<T<<endl;
    cerr<<"s0="<<s0<<endl;
    cerr<<"a="<<a<<endl;
    cerr<<"b="<<b<<endl;
    cerr<<"accIDM="<<accIDM<<endl;
    exit(-1);
  }

  return accIDM;
}

// HDM   (v0,T,s0,a,b,T_react) (presently just IDM with s,v,...=old quantities)

double accHDM(double s, double v, double dv, double aOld, const double beta[]){
  double v0=beta[0];
  double  T=maxsmooth(0,beta[1], 0.1);

  double s0=maxsmooth(0,beta[2], 0.1); 
  double  a=maxsmooth(0,beta[3], 0.1); 
  double  b=maxsmooth(0,beta[4], 0.1); 
  double Treact=beta[5];// not yet used here !! -> see micGlobalFunc
  //double Treact=maxsmooth(0,beta[5], 0.001); // not yet used !!

  //double  b=beta[4];

  double sstar=s0+max(0.,v*T+0.5*v*dv/sqrt(a*b));

  double accHDM=a*(1-pow(v/v0,4) - pow( sstar/s, 2));
  if(isnan(accHDM)){
    cerr<<"s="<<s<<" v="<<v<<" dv="<<dv<<endl;
    cerr<<"v0="<<v0<<endl;
    cerr<<"T="<<T<<endl;
    cerr<<"s0="<<s0<<endl;
    cerr<<"a="<<a<<endl;
    cerr<<"b="<<b<<endl;
    cerr<<"accHDM="<<accHDM<<endl;
    exit(-1);
  }

  return accHDM;
}


// IDM with v0 fixed (T,s0,a,b)

double accIDM_v0fixed(double s, double v, double dv, double v0, const double beta[]){
  double  T=maxsmooth(0,beta[0], 0.1);

  double s0=maxsmooth(0,beta[1], 0.1); 
  double  a=maxsmooth(0,beta[2], 0.1); 
  double  b=maxsmooth(0,beta[3], 0.1); 
  double sstar=s0+max(0.,v*T+0.5*v*dv/sqrt(a*b));
  return a*(1-pow(v/v0,4) - pow( sstar/s, 2));
}


// IIDM-ACC with a_lead unknown (=0) (v0,T,s0,a,b,cool)

double accACC(double s, double v, double dv, double a_lead, const double beta[]){
  double delta=4.;
  //double delta=beta[6];
  double cool=0.5*(1.+tanh(beta[5]-1.)); // maps to values between 0 and 1
  double v0=(beta[0]<100) ? beta[0] : beta[0]+log(beta[0]/100);  // NaN for pow((v/v0),delta) and extreme v0
  double  T=maxsmooth(0,beta[1], 0.1);

  double s0=maxsmooth(0,beta[2], 0.1); 
  double  a=maxsmooth(0,beta[3], 0.1); 
  double  b=maxsmooth(0,beta[4], 0.1); 

  double sstar=s0+max(0.,v*T+0.5*v*dv/sqrt(a*b));

  double z=sstar/max(s,0.01);
  double accEmpty=(v<=v0) ? a*(1- pow((v/v0),delta))
    : -b*(1- pow((v0/v),a*delta/b));
  double accPos=accEmpty*(1.- pow(z, min(2*a/accEmpty, 100.))  );
  double accInt=a*(1-z*z);

  double accIIDM=(v<v0) 
    ?  (z<1) ? accPos : accInt 
    :  (z<1) ? accEmpty : accInt+ accEmpty;

  double accIDM=a*(1-pow(v/v0,4) - pow( sstar/s, 2));//test with IDM

  //double a_lead_loc=0;
  double a_lead_loc=0; //!!!
  //double a_lead_loc=a_lead;
  //double a_lead_loc=max(a_lead,0.);
  //double a_lead_loc=min(a_lead,0.);

  double dvp=max(dv, 0.0);
  double v_lead = v-dvp;
  double denomCAH        =  v_lead*v_lead - 2 * s * a_lead_loc;

  double accCAH   = ( (v_lead*dvp  < - 2 * s * a_lead_loc) &&(denomCAH!=0))
    ? v*v*a_lead_loc/denomCAH
    : a_lead_loc - 0.5*dvp*dvp/max(s, 0.0001);
    
  // mix IIDM-CAH


  accIIDM=accIDM; //!!! test
  double accACC_IIDM=(accIIDM>accCAH)
    ? accIIDM
    : (1-cool)*accIIDM + cool*( accCAH+b*tanh((accIIDM-accCAH)/b));

  if(isnan(accACC_IIDM)){
    cerr<<"s="<<s<<" v="<<v<<" dv="<<dv<<endl;
    cerr<<"v0="<<v0<<endl;
    cerr<<"T="<<T<<endl;
    cerr<<"s0="<<s0<<endl;
    cerr<<"a="<<a<<endl;
    cerr<<"b="<<b<<endl;
    cerr<<"cool="<<cool<<endl;

    cerr<<"accACC_IIDM="<<accACC_IIDM<<endl;
    exit(-1);
  }
  return accACC_IIDM;
  //return accIIDM; //test
  //return accIDM; // test
 }


// VDIFF with triangular fundamental diagram (v0,T,s0,a,gamma)

double accVDIFF(double s, double v, double dv, const double beta[]){
  double v0=(beta[0]<100) ? beta[0] : beta[0]+log(beta[0]/100);
  double  T=beta[1];
  double s0=maxsmooth(0,beta[2], 0.1); 
  double  a=maxsmooth(0,beta[3], 0.1);
  double gamma=beta[4]; // sensitivity to speed differences

  double vopt=max( min((s-s0)/T, v0), 0.);  // optimal velocity
  return a*(vopt-v)/v0 - gamma*dv;
  //return a*(vopt-v)/v0;

}


// OVM with triangular fundamental diagram (v0,T,s0,a)

double accOVM(double s, double v, double dv, const double beta[]){
  double betaOVM[4];
  for (int k=0; k<4; k++){betaOVM[k]=beta[k];} betaOVM[4]=0;
  return accVDIFF(s,v,dv,betaOVM);
}


// Gipps (v0,T,s0,a,b)

double accGipps(double s, double v, double dv, const double beta[]){
  double v0=beta[0];
  double T=maxsmooth(0,beta[1], 0.1);
  double s0=maxsmooth(0,beta[2], 0.1); 
  double  a=maxsmooth(0,beta[3], 0.1);
  double  b=maxsmooth(0,beta[4], 0.1); 

  double vp=v-dv;
  //T=1.1;

  double vsafe=-b*T+sqrt(b*b*T*T+vp*vp+2*b*max(s-2*s0,0.)); // safe velocity!!!2*s0 statt s0
  double vnew=min(vsafe, min(v+a*T, v0));
  return (vnew-v)/T;
}


// ADAS v5 (v0,T,s0,a,c1)

double accADAS(double s, double v, double dv, const double beta[]){
  double v0=(beta[0]<100) ? beta[0] : beta[0]+log(beta[0]/100);  // NaN for pow((v/v0),delta) and extreme v0
  double   T=maxsmooth(0,beta[1], 0.1);
  double s0=maxsmooth(0,beta[2], 0.1); 
  double   a=maxsmooth(0,beta[3], 0.1);
  double  c1=beta[4];
  double eta=beta[5];

  double sloc=max(s,0.001);
  double vd=max(0., (s-s0)/T);
  double vopt=min(v0, vd);
  double aOVM=a/v0*(vopt-v);
  return min(aOVM, a/v0*(vd-v) -2*c1*dv*exp(s0/sloc)*(1+s0*dv/(eta*sloc*sloc)));
}





// ####################################################
/* Function for local calibration 
   hata are n IDM accelerations for the exogeneous variables  
   s[i]=data[i],
   v[i]=data[Ndata+i],
   vL[i]=data[2*Ndata+i],
   aL[i]=data[3*Ndata+i],
   i=0 ... n-1
   will be fitted to ydata = n accelerations
*/
// ####################################################

// I have modified void* data and void* adata to double* data and double* adata, 
// respectively, everywhere since I do not know how to use the pointer-to-void structure

//######################################


// estimated endogeneous variable = acceleration hatacc
void micLocalFunc(double *beta, double *hatacc, int Mparam, int ndata, double *data){
  register int i;
  double sdata[ndata];
  double vdata[ndata];
  double vldata[ndata];
  double aldata[ndata];
  int choice_model=0; // only default

  for(i=0; i<ndata; ++i){
    sdata[i]=data[i];
    vdata[i]=data[ndata+i];
    vldata[i]=data[2*ndata+i];
    aldata[i]=data[3*ndata+i];
    choice_model=(int)(data[4*ndata]);
  }
  Statistics stat;
  double v0=2*stat.getmax(vdata,ndata); //!! for restrained calibration (model 1)

  if( choice_model==0)
    for(i=0; i<ndata; ++i) hatacc[i]=accIDM(sdata[i],vdata[i],vdata[i]-vldata[i], beta); 
  else if( choice_model==1)
    for(i=0; i<ndata; ++i) hatacc[i]=accIDM_v0fixed(sdata[i],vdata[i],vdata[i]-vldata[i], v0, beta);
  else if( choice_model==2)
    for(i=0; i<ndata; ++i) hatacc[i]=accACC(sdata[i],vdata[i],vdata[i]-vldata[i], aldata[i], beta); 
  else if( choice_model==3)
    for(i=0; i<ndata; ++i) hatacc[i]=accGipps(sdata[i],vdata[i],vdata[i]-vldata[i], beta); 
  else if( choice_model==4)
    for(i=0; i<ndata; ++i) hatacc[i]=accOVM(sdata[i],vdata[i],vdata[i]-vldata[i], beta); 
  else if( choice_model==5)
    for(i=0; i<ndata; ++i) hatacc[i]=accVDIFF(sdata[i],vdata[i],vdata[i]-vldata[i], beta); 
  else if( choice_model==6)
    for(i=0; i<ndata; ++i) hatacc[i]=accADAS(sdata[i],vdata[i],vdata[i]-vldata[i], beta); 
  else if( choice_model==7)
    for(i=0; i<ndata; ++i) hatacc[i]=accHDM(sdata[i],vdata[i],vdata[i]-vldata[i],aldata[i], beta); 
  else{ cerr<<"Error: choice_model="<<choice_model<<" mus be <= 7"<<endl; exit(-1);}
}


// estimated endogeneous variable = gap hats
void micGlobalFunc(double *beta, double *hats_or_lns, int Mparam, int ndata, double *data){

  //cout <<"beta[5]="<<beta[5]<<endl;

  // extract exogeneous and ctrl. data

  double sdata[ndata];
  double vdata[ndata];
  double vldata[ndata];
  double aldata[ndata];

  int choice_model=0; // only default
  double dtData=0; // only default
  int calType=1; // only default

  for(int i=0; i<ndata; ++i){
    sdata[i]=data[i];
    vdata[i]=data[ndata+i];
    vldata[i]=data[2*ndata+i];
    aldata[i]=data[3*ndata+i]; //!!!
  }
  choice_model=(int)(data[4*ndata]);
  dtData=data[4*ndata+1];
  calType=(int)(data[4*ndata+2]);


  // make vl and a consistent => redefine vldata and adata by basic sdata and vdata

  if(true){
    vldata[0]=vdata[0]+(sdata[1]-sdata[0])/dtData;
    vldata[ndata-1]=vdata[ndata-1]+(sdata[ndata-1]-sdata[ndata-2])/dtData;
    for (int i=1; i<ndata-1; i++){
      vldata[i]=vdata[i]+0.5*(sdata[i+1]-sdata[i-1])/dtData;
      //vldata[i]=vdata[i]+(sdata[i+1]-sdata[i])/dtData;
    }

    aldata[0]=(vldata[1]-vldata[0])/dtData;
    aldata[ndata-1]=(vldata[ndata-1]-vldata[ndata-2])/dtData;
    for (int i=1; i<ndata-1; i++){
      aldata[i]=0.5*(vdata[i+1]-vdata[i-1])/dtData
	+(sdata[i+1]-2*sdata[i]+sdata[i-1])/(dtData*dtData);
    }
    

  }



  // define time step used for simulation (may be different from timestep dtData of data)
  // Gipps: dt=T=Tr=tau_relax; OVM: smaller dt

  const double dtSim=(choice_model==4) ? 0.1 : (choice_model==3) ? beta[1] : 0.2; //!!
  //const double dtSim=(choice_model==4) ? 0.1 : (choice_model==3) ? 0.1 : 0.2; 
  //const double dtSim=(choice_model==4) ? 0.1 : 0.2; //!!

  // initialize

  const int nSim=(int)((ndata-1)*dtData/dtSim);
  double tmax=(ndata-1)*dtData;
  double sSim[nSim]; // sSim: nSim points; hats_or_lns (=endog.argument): ndata points
  double vSim[nSim];

  sSim[0]=sdata[0];
  vSim[0]=vdata[0];
  double accmax=10; //!! for tests for new targets (accmax=100)
  Statistics stat;
  double v0=stat.getmax(vdata,ndata); // for restrained minimum

  double sSimTr=sSim[0]; // initialization; only used for finite reaction time
  double vSimTr=vSim[0]; // initialization; only used for finite reaction time
  double aSimTr=0;
  double vlDataTr=vldata[0];

  // micro-simulation 

  for (int i=1; i<nSim; i++){

    double tSimLast=(i-1)*dtSim;
    double tSim=i*dtSim;
    if(tSim>tmax){
      cerr<<"ndata="<<ndata<<" i="<<i<<" nSim="<<nSim<<": tSim="<<tSim<<">tmax="<<tmax<<endl;
      exit(-1);
    }
    if((sSim[i-1]<0) || (vSim[i-1]<0)){
      cerr<<"i="<<i<<"sSim[i-1]="<<sSim[i-1]<<" vSim[i-1]="<<vSim[i-1]
	  <<" vSim[i-1]-vldata[i-1]="<<vSim[i-1]-vldata[i-1] <<endl;
      for (int k=0; k<Mparam; k++){
	cerr<<"beta["<<k<<"]="<<beta[k]<<endl;
      }
      //exit(-1);
    }

    double vintpOld=intp(vdata,ndata,tSimLast,0,tmax);
    double vlintpOld=intp(vldata,ndata,tSimLast,0,tmax);
    double alintpOld=intp(aldata,ndata,tSimLast,0,tmax);
    double sintpOld=intp(sdata,ndata,tSimLast,0,tmax);

    double sintp=intp(sdata,ndata,tSim,0,tmax);
    double vintp=intp(vdata,ndata,tSim,0,tmax);
    double vlintp=intp(vldata,ndata,tSim,0,tmax);
    if(choice_model==7){
      //double beta5Loc=maxsmooth(0,beta[5], 0.001); //!!
      double beta5Loc=beta[5]; //!!
      double tDelay=min(tmax, max(0., tSimLast-beta5Loc));
      //sSimTr=(i<3) ? sSim[i-1] : intp(sSim,i,tDelay,0,tSimLast);
      //vSimTr=(i<3) ? vSim[i-1] : intp(vSim,i,tDelay,0,tSimLast);
      sSimTr=sSim[i-1]; vSimTr=vSim[i-1];
      vlDataTr=(i<3) ? vlintp : intp(vldata,ndata,tDelay,0,tmax);
      //aSimTr=(i<3) ? 0 : (vSim[i-1]-vSim[i-2])/dtSim;
      if(fabs(accIDM(sSimTr,vSimTr,vSimTr-vlDataTr,beta)-accIDM(sSim[i-1],vSim[i-1],vSim[i-1]-vlintpOld,beta))>5){
	cerr<<"iSimLast="<<i-1<<" tSimLast="<<tSimLast
	    <<" beta[5]="<<beta[5]<<" beta5Loc="<<beta5Loc
	    <<" tDelay="<<tDelay<<endl;
        cerr<<" sSimTr="<<sSimTr<<" sSim[i-1]="<<sSim[i-1]<<endl;
        cerr<<" vSimTr="<<vSimTr<<" vSim[i-1]="<<vSim[i-1]<<endl;
        cerr<<" vlDataTr="<<vlDataTr<<" vlintpOld="<<vlintpOld<<endl;
	cerr<<" accIDM(sSimTr,..)="<<accIDM(sSimTr,vSimTr,vSimTr-vlDataTr,beta)
	    <<" accIDM(sSim[i-1],..)="
	    <<accIDM(sSim[i-1],vSim[i-1],vSim[i-1]-vlintpOld,beta)<<endl;
      }

    }


    double acc=
     (choice_model==0)    ?   accIDM(sSim[i-1],vSim[i-1],vSim[i-1]-vlintpOld,beta)
      : (choice_model==1) ? accIDM_v0fixed(sSim[i-1],vSim[i-1],vSim[i-1]-vlintpOld,15,beta)
      : (choice_model==2) ?   accACC(sSim[i-1],vSim[i-1],vSim[i-1]-vlintpOld,alintpOld,beta)
      : (choice_model==3) ? accGipps(sSim[i-1],vSim[i-1],vSim[i-1]-vlintpOld,beta)
      : (choice_model==4) ?   accOVM(sSim[i-1],vSim[i-1],vSim[i-1]-vlintpOld,beta)
      : (choice_model==5) ? accVDIFF(sSim[i-1],vSim[i-1],vSim[i-1]-vlintpOld,beta)
      : (choice_model==6) ?  accADAS(sSim[i-1],vSim[i-1],vSim[i-1]-vlintpOld,beta)
      : (choice_model==7) ?  accHDM(sSimTr,vSimTr,vSimTr-vlDataTr, aSimTr,beta)
      : 0;
    vSim[i]=max(0., vSim[i-1]+acc*dtSim);
    // check in which sense the observed s is platoon consistent!!!

    sSim[i]=max(0.1, sSim[i-1] + 0.5*dtSim*(vlintp+vlintpOld)
			- 0.5*dtSim*(vSim[i]+vSim[i-1]));

    if((choice_model==3) && (sSim[i]<0.1*beta[2])){sSim[i]=0.1*beta[2];vSim[i]=0;}//!!!
    // test if new target detected
    double sdiff=sintp-sintpOld -( 0.5*dtSim*(vlintp+vlintpOld-vintp-vintpOld));
    if(fabs(sdiff)>0.5*accmax*dtSim*dtSim){
      if(false){
        cout <<"i="<<i<<" tSim="<<tSim<<" sdiff="<<sdiff<<">0.5*accmax*dtSim*dtSim="<<0.5*accmax*dtSim*dtSim
	   <<":New target detected!"<<endl;
        cout <<"  sintp-sintpOld="<<sintp-sintpOld <<endl;
        cout <<"  0.5*dtSim*(vldata-vdata)="<<0.5*dtSim*(vlintp-vintp)<<endl;
        cout <<"  0.5*dtSim*(vldataOld-vdataOld)="<<0.5*dtSim*(vlintpOld-vintpOld)<<endl;
      }
      sSim[i]=sintp;
      vSim[i]=vintp;
    }
   
  }

  // define test variable hats by interpolation of sSim at data time instants

  if(calType==1){
    for (int i=0; i<ndata; i++){
      hats_or_lns[i]=intp(sSim, nSim, i*dtData, 0, tmax); // hats
    }
  }
  else{// calType=2
    for (int i=0; i<ndata; i++){
      hats_or_lns[i]=log(intp(sSim, nSim, i*dtData, 0, tmax)); // ln(hats)
    }
  }



} // end micGlobalFunc

// Calculate objective function (the obove functions calculate hats_or_lns[i]
// now only for global calibration, i.e., y=gap s or log(gap s)
//(use ../localCalibrationBook for local plots)


double objFun(
	      void (*func)(double *beta, double *haty, int Mparam, int ndata, double *data),
	      double *beta, int Mparam, int ndata, double *data
	      ){
  double haty[ndata];
  double ydata[ndata]; // local calibr: adata->dvdata/dt; global calibr: sdata

  int calType=(int)(data[4*ndata+2]);
  double dtData=data[4*ndata+1];

  if(calType==0){ //ydata=adata=>vdata; vdata= second array of data[]
    ydata[0] =(data[ndata+1]-data[ndata])/dtData;
    ydata[ndata-1] =(data[2*ndata-1]-data[2*ndata-2])/dtData;
    for(int i=1; i<ndata-1; i++){ydata[i] =0.5*(data[ndata+i+1]-data[ndata+i-1])/dtData;}
  }

  else if(calType==1){ //ydata=ydata, first array of data[]
    for(int i=0; i<ndata; ++i){ydata[i]=data[i];}
  }

  else if(calType==2){//ydata=ydata, first array of data[]
    for(int i=0; i<ndata; ++i){ydata[i]=log(data[i]);}
  }


  func(beta,haty,Mparam,ndata,data);


  double obj=0;
  for(int i=0; i<ndata; ++i){
    obj+=(haty[i]-ydata[i])*(haty[i]-ydata[i]);
  }
  return obj;
}




//#############################################

int main(int argc, char* argv[]){

  // ####################################################
  // input: File names
  // ####################################################

  if (argc!=5){
    cout <<"\nCalling sequence: calibTraj <calType> <modelnumber> <EFC datafile> <calcObjLandscape>\n"; 
    //cout <<" projectName contains a .ctrl file for input\n";
    //cout <<" and a .cal* file indicating the calibr results\n";
    cout <<" calType={0 (local), 1 (global, obj=s), 2 (global, obj=ln s)}\n";
    cout <<" modelnumber: 0=IDM, 1=IDM_v0fixed, 2=ACC, 3=Gipps, 4=OVM, 5=VDIFF, 6=ADAS, 7=IDMdelay\n";
    cout <<" EFC data in Bosch format:\n";
    cout <<" # t(s)    vL(m)   v(m/s)   a(m^2/s)      s(m)    dv(m/s) \n";
    cout <<" calcObjLandscape=1 if objective function landscape calulated (runtime!), 0 if not"<<endl;
    cout  <<"Exampe if locally calibrating the IDM to the EFC data file Boschdata3 inside project localProjGen:\n";
    cout <<" calibTraj 0 0 Boschdata3 0\n";
    exit (-1);
  }

  char EFCname[256];
  //char projName[256];
  //char ctrlName[256];
  char outName[256];

  //sprintf(projName,"%s",argv[1]);
  int calType=atoi(argv[1]);
  int choice_model=atoi(argv[2]);
  sprintf(EFCname,"%s",argv[3]);
  //sprintf(ctrlName,"%s.ctrl",projName);
  int calcObjLandscape=atoi(argv[4]);

  const char* model0="IDM";
  const char* model1="IDM-v0fix";
  const char* model2="IIDM-ACC";
  const char* model3="Gipps";
  const char* model4="OVM";
  const char* model5="VDIFF";
  const char* model6="ADAS";
  const char* model7="IDM-delay";
  const char* modstring[]={model0, model1, model2, model3, model4, model5, model6, model7};

  sprintf(outName,"cal%s_%s_%s",
	  ((calType==0) ? "Loc" : ((calType==1 ? "Glob_s" : "Glob_lns"))),
	  modstring[choice_model],EFCname);


  // ####################################################
  // input: Chose model
  // ####################################################

  //{0=IDM(5 params), 1=IDM_v0fixed(4), 2=ACC(6), 3=Gipps(5), 4=OVM(4),
  // 5=VDIFF(5), 6=ADAS(6)}
 
  const int Mparam=
    (choice_model==0) ? 5 :
    (choice_model==1) ? 4 :
    (choice_model==2) ? 6 :
    (choice_model==3) ? 5 :
    (choice_model==4) ? 4 :
    (choice_model==5) ? 5 :
    (choice_model==6) ? 6 : 6;

  const char* pstring0=(choice_model==1) ? "T" : "v0";
  const char* pstring1=(choice_model==1) ? "s0" : "T";
  const char* pstring2=(choice_model==1) ? "a" : "s0";
  const char* pstring3=(choice_model==1) ? "b" : "a";
  const char* pstring4=(choice_model==5) ? "gamma" : (choice_model==6) ? "c1" : "b";
  const char* pstring5=(choice_model==2) ? "cool" :
    (choice_model==6) ? "eta" :
    (choice_model==7) ? "Treact" :
    "";

  const char* pstring[]={pstring0, pstring1, pstring2, pstring3, pstring4, pstring5};
 
  double beta[Mparam]; // parameter vector

  // ####################################################
  // input: define control variables
  // ####################################################

  InOut inout;

  // later possibly ctrl variables

  //double ctrlVals[256];
  //int nCtrl=0;
  //inout.get_array(ctrlName,nCtrl,ctrlVals);



  // following may be put into .ctrl file of project

  double data[4*Ndata+2]; // exogeneous data vector+ some ctrl data
  double opts[LM_OPTS_SZ]; // control parameters for dlevmar_[dif|der] (input)
  double info[LM_INFO_SZ]; // information on the calculation of results (output) 
  const int n_iter_max=1000; // maximum number of iterations



 
  // ####################################################
  // exogeneous and endogeneous measurement data
  // ####################################################


  double tdata[Ndata];
  double sdata[Ndata]; // endogeneous for global calibration
  double lnsdata[Ndata]; // endogeneous for global calibration
  double vdata[Ndata];
  double vldata[Ndata];
  double adata[Ndata];  // endogeneous for local calibration
  double aldata[Ndata];
  int ndata;

  cout <<"Check/verify that n_array<=Ndata="<<Ndata<<"!\n";
  
  inout.get_col (EFCname, 1,  ndata, tdata);
  if(ndata>Ndata){
     cerr <<"Error when reading EFC data: ndata="<<ndata<<" greater than Ndata="<<Ndata<<"!\n";
     exit(-1);
  }
  inout.get_col (EFCname, 3,  ndata, vdata);
  inout.get_col (EFCname, 5,  ndata, sdata);
  inout.get_col (EFCname, 2,  ndata, vldata);

  // prevent netto gaps=0; set to minimum value smin
  double smin=0.9;  // 0.5 or 0.888
  for (int i=0; i<ndata; i++){if (sdata[i]<smin) sdata[i]=smin;}
  // prevent negative speeds
  for (int i=0; i<ndata; i++){vdata[i]=max(0., vdata[i]);}
  //define lnsdata
  for (int i=0; i<ndata; i++){lnsdata[i]=log(sdata[i]);}


 
  // calculate accelerations from speed differences
  // (for objective fun of local calibr and accACC model)

  double dt=tdata[1]-tdata[0];
  if( !( (dt>1e-6)||(dt<5))){
      cerr<<"after reading data: in lnL: error: dt not in range between 1e-6 and 5\n"
	<<"Something went wrong!"<<endl;
      exit(-1);
  }

  adata[0] =(vdata[1]-vdata[0])/dt;
  aldata[0]=(vldata[1]-vldata[0])/dt;
  adata[ndata-1] =(vdata[ndata-1]-vdata[ndata-2])/dt;
  aldata[ndata-1]=(vldata[ndata-1]-vldata[ndata-2])/dt;
  for (int i=1; i<ndata-1; i++){
    adata[i]=0.5*(vdata[i+1]-vdata[i-1])/dt;
    aldata[i]=0.5*(vldata[i+1]-vldata[i-1])/dt;
  }

  // wrap up in data container
 
  for(int i=0; i<ndata; i++){
    data[i]=sdata[i];
    data[ndata+i]=vdata[i];
    data[2*ndata+i]=vldata[i];
    data[3*ndata+i]=aldata[i]; //!!! aldata (correct) or vldata (better)
  }

  // add ctrl data to the exogeneous data

  data[4*ndata]=choice_model;
  data[4*ndata+1]=dt;
  data[4*ndata+2]=calType;






  // ####################################################
  // initialize parameter vector
  // ####################################################

  Statistics stat;

  double Tdata[Ndata];
  double vthr=5;
  for (int i=0; i<ndata; i++){Tdata[i]=sdata[i]/max(vdata[i],vthr);}

  beta[0]=stat.getmax(vdata,ndata);  //v0
  beta[1]=0.7*stat.arithmeticMeans(Tdata,ndata);  //T
  beta[2]=stat.getmin(sdata,ndata);  //s0
  beta[3]=0.5*stat.getmax(adata,ndata); // a;  adata=accdata
  beta[4]=-0.5*stat.getmin(adata,ndata);  //b
  if(choice_model==5){beta[4]=0.2;} //VDIFF: speed diff sensitivity gamma
  if(choice_model==6){beta[4]=0.1;} //ADAS v5: c1
  beta[5]=(choice_model==2) ? 0.5 : 0.2; // cool or eta (model 6) or Treact(model7)

  if(choice_model==1){// IDM w/o v0
    for(int k=1; k<5; k++){beta[k-1]=beta[k];}
  }

  cout << "initial estimate of parameter vector:\n";
  for (int ip=0; ip<Mparam; ip++){
    cout <<"beta["<<ip<<"]="<<beta[ip]<<" "<<pstring[ip]<<endl;
  }

  // ####################################################
  // optimization control parameters
  // ####################################################

  // passing to levmar NULL instead of opts reverts to defaults 
  opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
  opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used 



  // ####################################################
  // invoke the optimization function
  // ####################################################



  double covar[Mparam*Mparam];
  int ret; // return signal of the dlevmar functions

  if(calType==0)
    ret= dlevmar_dif(micLocalFunc, beta, adata, Mparam, ndata, 
		     n_iter_max, opts, info, NULL, covar, data);
  else if(calType==1)
    ret= dlevmar_dif(micGlobalFunc, beta, sdata, Mparam, ndata, 
		     n_iter_max, opts, info, NULL, covar, data);
  else if(calType==2)
    ret= dlevmar_dif(micGlobalFunc, beta, lnsdata, Mparam, ndata, 
		     n_iter_max, opts, info, NULL, covar, data);

  double betaFinal[Mparam];
  for (int k=0; k<Mparam; k++){betaFinal[k]=beta[k];}

  double SSEinit=info[0];
  double SSEfinal=info[1];

  // ####################################################
  // print result
  // ####################################################

  char out_string[5000];
  if(calType==0)sprintf(out_string,"#Results of fitting %s locally to %s\n",
     modstring[choice_model],EFCname);
  else if(calType==1)sprintf(out_string,
     "#Results of fitting %s globally to %s with respect to gaps\n",
     modstring[choice_model],EFCname);
  else if(calType==2)sprintf(out_string,
     "#Results of fitting %s globally to %s with respect to ln(gaps)\n",
     modstring[choice_model],EFCname);
  sprintf(out_string,"%s#=============================================\n",out_string);
  cout <<"\nResults of fitting "<<modstring[choice_model]<<" to "<<EFCname<<":\n";
  cout <<"==================================================\n";
 
  char reasonstr[256];
  sprintf(reasonstr,"");
  if(info[6]==1){sprintf(reasonstr,"small gradient");}
    else if(info[6]==2){sprintf(reasonstr,"small param change");}
    else if(info[6]==3){sprintf(reasonstr,"itmax reached");}
    else if(info[6]==4){sprintf(reasonstr,"singular matrix");}
    else if(info[6]==7){sprintf(reasonstr,"NaN or Inf");}
    else{sprintf(reasonstr,"other reason: %g",info[6]);}
  char relStringInit[256];
  char relStringFinal[256];
  if(calType==0){
    sprintf(relStringInit,"(avg error %1.2f m/s^2)", sqrt(info[0]/ndata));
    sprintf(relStringFinal,"(avg error %1.2f m/s^2)", sqrt(info[1]/ndata));
  }
  else if(calType==1){
    sprintf(relStringInit,"(avg error %2.1f m)", sqrt(info[0]/ndata)); 
    sprintf(relStringFinal,"(avg error %2.1f m)", sqrt(info[1]/ndata));
  }
  else if(calType==2){
    sprintf(relStringInit,"(avg error %1.2f \% )", 100*sqrt(info[0]/ndata));
    sprintf(relStringFinal,"(avg error %1.2f \% )", 100*sqrt(info[1]/ndata));
  }

  sprintf(out_string,"%s\n#Levenberg-Marquardt returned in %g iter;  reason: %s\n",
     out_string, info[5], reasonstr);
  sprintf(out_string,"%s# number of function evaluations: %g\n", out_string,info[7]);
  sprintf(out_string,"%s# number of Jakobian evaluations: %g\n\n", out_string,info[8]);
  sprintf(out_string,"%s# resulting SSE: %g %s\n", out_string,info[1],relStringFinal);
  sprintf(out_string,"%s# initial  SSE: %g %s\n", out_string,info[0],relStringInit);

  char fit_string[4096];  
  sprintf(fit_string,"\n##Best fit parameters:\n");
  double stddev[Mparam];

  // de-scale some bounded quantities

  //covar[Mparam*2+2]*=pow(0.5+0.25*beta[2]/sqrt(0.25*beta[2]*beta[2]+0.1),2);
  betaFinal[1]=maxsmooth(0,beta[1], 0.1); // covar before redef beta!
  betaFinal[2]=maxsmooth(0,beta[2], 0.1); // covar before redef beta!
  betaFinal[3]=maxsmooth(0,beta[3], 0.1); // covar before redef beta!

  if(choice_model<4){
    //covar[Mparam*4+4]*=pow(0.5+0.25*beta[4]/sqrt(0.25*beta[4]*beta[4]+0.1),2);
    betaFinal[4]=maxsmooth(0,beta[4], 0.1); // covar before redef beta!
  }

  if(choice_model==2){
    //covar[Mparam*5+5]/=pow(cosh(beta[5]-1),2);
    betaFinal[5]=0.5*(1.+tanh(beta[5]-1.));
  }

  //if(choice_model==3){betaFinal[1]=1;} //!!!

  if(choice_model==7){
    //covar[Mparam*5+5]...
    //betaFinal[5]=maxsmooth(0,beta[5], 0.001);
    betaFinal[5]=beta[5];
  }


  for (int im=0; im<Mparam; im++){
    stddev[im]=sqrt(covar[Mparam*im+im]);
    sprintf(fit_string,"%s#%s=\t%2.5f; sig_%s=%2.5f\n",
	    fit_string,pstring[im],betaFinal[im],pstring[im], stddev[im]);
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


  // write result in file

  // get local accelerations/global trajectories for estimated parameters

  double hataFinal[Ndata];
  micLocalFunc(betaFinal, hataFinal, Mparam, ndata, data);
  double hats_or_lnsFinal[Ndata];
  micGlobalFunc(betaFinal, hats_or_lnsFinal, Mparam, ndata, data);

  sprintf(out_string,"%s%s%s", out_string, fit_string, corr_string);
  if(calType==0){
    inout.write_array(outName, ndata, tdata, sdata, vdata, vldata, adata, hataFinal, out_string, 
		      "\n\n#time\t\tsdata\t","vdata\t","vldata\t","adata\t", "hata");
  }
  else if(calType==1){
    inout.write_array(outName, ndata, tdata, sdata, vdata, vldata, hats_or_lnsFinal, out_string, 
		    "\n\n#time\t\tsdata\t","vdata\t","vldata\t", "hats");
  }
  else if(calType==2){
    inout.write_array(outName, ndata, tdata, sdata, vdata, vldata, lnsdata, hats_or_lnsFinal, out_string, 
		    "\n\n#time\t\tsdata\t","vdata\t","vldata\t","lnsdata\t", "hatlns");
  }

  cout <<out_string<<endl;
  cout <<"\nwriting results to "<<outName<<endl;


  //############################################
  // Write log-likelihood around the  maximum to files
  //############################################

  if(calcObjLandscape==1){


  double w_stddev=30; // half-width of scanning range in stddev
  int nout=41; // number of grid elements in either direction
  double v0min=8;
  double v0max=25;
  char objFunFilename[1024];
  char titleString[1024];
  double objFunData[inout.NYMAX][inout.NYMAX];

  double dbeta[Mparam];
  double betamin[Mparam];
  double betamax[Mparam];


  for (int j=0; j<Mparam; j++){
    betamin[j]=betaFinal[j]-w_stddev*stddev[j];
    betamax[j]=betaFinal[j]+w_stddev*stddev[j];
    if(betamin[j]<stddev[j]){betamin[j]=stddev[j];}
  }
 
  // special for beta[0]=v0 
  betamin[0]=max(betamin[0],v0min);
  betamax[0]=min(betamax[0],v0max);

  // or explicit

  //if(calType>0){
  if(true){
    betamin[0]=10;  betamax[0]=25;
    betamin[1]=0.7; betamax[1]=2.0;
    betamin[2]=0.1; betamax[2]=4;
    betamin[3]=0.2; betamax[3]=3.5;
    betamin[4]=0.2; betamax[4]=3.5;
    betamin[5]=0.0; betamax[5]=1;
  }

  for (int j=0; j<Mparam; j++){
    dbeta[j]=(betamax[j]-betamin[j])/(nout-1);
  }


  // make (Mparam-1)*(Mparam-2) data sets and files

  for (int ibeta=0; ibeta<Mparam-1; ibeta++){
    for (int jbeta=ibeta+1; jbeta<Mparam; jbeta++){

      // generate the data set for a given beta combination

      for (int k=0; k<Mparam; k++){beta[k]=betaFinal[k];}// revert non-used dimensions
      for (int i=0; i<nout; i++){
	for (int j=0; j<nout; j++){
	  beta[ibeta]=betamin[ibeta]+i*dbeta[ibeta];
	  beta[jbeta]=betamin[jbeta]+j*dbeta[jbeta];
	  objFunData[i][j]=(calType==0)
	    ? objFun(micLocalFunc,beta, Mparam, ndata, data)
	    : objFun(micGlobalFunc,beta, Mparam, ndata, data);
	}
      }


      // write the file for this combination

      sprintf(objFunFilename,"%s_beta%i_beta%i",outName, ibeta, jbeta);
      sprintf(titleString,"#Objective function for beta%i and beta%i",ibeta, jbeta);
      sprintf(titleString,"%s\n#Base values: beta0=%.1f, beta1=%.2f, beta2=%.2f, beta3=%.2f, beta4=%.2f",
	      titleString,betaFinal[0],betaFinal[1],betaFinal[2],betaFinal[3],betaFinal[4]);

      inout.write_array2d(objFunFilename, betamin[ibeta],betamax[ibeta],nout,
			  betamin[jbeta],betamax[jbeta],nout,
			  objFunData,titleString);
    }
  }
  }
  exit(0);
}

