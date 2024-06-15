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

//  c++ 
#include <iostream>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

#include "levmar-2.5/levmar.h"

#include "InOut.h"
#include "Math.h"
#include "Statistics.h"
//#include "RandomUtils.h"
#include "general.h"



//######################################
// user defined global settings/functions
//######################################


const int MPARAM_MAX=10; // max number of params
const double DT=0.2; // sim timestep (newTarget detect needs dtData%dtSim=0)
                      // !!! adapt it such dtData=multiple
const double DT_OVM=0.05;
const double IDM_B=2;
const double GAP_MIN=0.2;  // need minimum gap for GOF=SSE(ln s)

const bool FVDM_FULL_LINEAR=false;

const double BMAX=9; // maximum braking deceleration (all models except OVM)

// use dt*0.5*(vlLast+vlNew) to update gaps.
// Otherwise, only follower's speed considered second order
// but setting it to false is generally more robust
// (no speed wiggling in combination with other special rules if gap<0)

const bool SECOND_ORDER_VLEAD=false; 

// use new data entry to set a signal for a new target
// (otherwise, last data entry used, see INTERNAL CALCULATION OF NEW TARGETS)
const bool USE_NEW_DATA_TO_DETERMINIE_TARGETS=true;
const bool DEBUG_TARGETS=false;

// Not using vlSpeed coumn and always calculate vLead internally
// makes d8_0900_0930_road2_veh145_IDM_s better,
// otherwise often slightly worse
// generally, set to false (=state < 2022-11-20)
// if the data vLead column is not platoon consistent
// with respect to gap and follower speeds (these are always used)
// see INTERNAL CALCULATION OF NEW TARGETS AND LEADING SPEED
const bool USE_VLINFO_IF_AVAILABLE=true; 


double get_v0fixed(const double vdata[], int ndata){
  Statistics stat;
  return 2.0*stat.getmax(vdata,ndata);
}

double get_dtSim(int choice_model){
  return (choice_model==4) ? DT_OVM 
    : DT;
}


//###############################################################
// in case of external restrictions, the global params used
// for the unconstrained optimization in dlevmar_dif(.) 
// differ from the local parms used in accMic and the underlying
// CF models
// betaLoc can have one element more: def double betaLoc[Mparam+1]
//###############################################################

void relocateParameters(int choice_model, int Mparam, 
			  double v0fix, double bfix, 
			  const double betaGlob[], double betaLoc[]){

   // default for unconstrained optimization: betaLoc=betaGlob

  for(int p=0; p<Mparam; p++){
    betaLoc[p]=betaGlob[p];
  }


  // the different restrained IDM variants MparamLoc>paramGlob=Mparam
  // (notice that IDM_delay is the reverse, a superset of params)

  if(choice_model==1){ // v0fixed
    betaLoc[0]=v0fix; 
    for (int p=0; p<Mparam; p++){betaLoc[p+1]=betaGlob[p]; }
  }

  else if(choice_model==8){ // v0 and b is fixed
    betaLoc[0]=v0fix;
    betaLoc[4]=bfix;
    for (int p=0; p<Mparam; p++){betaLoc[p+1]=betaGlob[p]; }
  }

  else if(choice_model==9){ // v0 and b is fixed
    betaLoc[4]=bfix;
    for (int p=0; p<Mparam; p++){betaLoc[p]=betaGlob[p]; }
  }

  if(false){
    cout <<"micGlobalFunc: choice_model="<<choice_model<<endl;
    cout <<" global function for levmar, local for calcul. obj fun"<<endl;
    for(int p=0; p<Mparam; p++){
      cout <<"Global params: p="<<p<<" betaGlob[p]="<<betaGlob[p]<<endl;
    }
    for(int p=0; p<6; p++){
      cout <<"Local params: p="<<p<<" betaLoc[p]="<<betaLoc[p]<<endl;
    }
  }
}



//###############################################################
// IDM  (v0,T,s0,a,b); formulated with vl=v_leader; dv=v-vl
// handles:
// choice_model=0: IDM
// choice_model=1: IDM_v0fix
// choice_model=8: IDM_v0bfix
// choice_model=9: IDM_bfix
// choice_model=7: IDM_delay (IDM gets old state, see micGlobalFunc)
//###############################################################

double accIDM(double s, double v, double vl, const double beta[]){
  // s already guaranteed to be >=GAP_MIN in the data processing stage 
  double v0=beta[0];
  double  T=maxsmooth(0,beta[1], 0.01);
  double s0=maxsmooth(0,beta[2], 0.01);
  double  a=maxsmooth(0,beta[3], 0.01); //Math.h -> maxsmooth to prevent
  double  b=maxsmooth(0,beta[4], 0.01); //neg params. CHECK!
  double sstar=max(GAP_MIN,s0+v*T+0.5*v*(v-vl)/sqrt(a*b)); //!!!
  double accIDM=max(-BMAX, a*(1-pow(v/v0,4) - pow( sstar/s, 2)));

  //if((v>5.661)&&(v<5.663)){ // without exit
  if(isnan(accIDM)){  //!! with exit
    cerr<<"s="<<s<<" v="<<v<<" vl="<<vl<<endl;
    cerr<<"v0="<<v0<<endl;
    cerr<<"T="<<T<<endl;
    cerr<<"s0="<<s0<<endl;
    cerr<<"a="<<a<<endl;
    cerr<<"b="<<b<<" beta[4]="<<beta[4]
	<<" maxsmooth(0.01,beta[4], 0.01)="
	<<maxsmooth(0.01,beta[4], 0.01)<<endl;
    cerr<<"sstar="<<sstar<<" accIDM="<<accIDM<<endl;
    //exit(-1);
  }

  return accIDM;
}



//###############################################################
// IIDM-ACC model with a_lead unknown (=0) (v0,T,s0,a,b,cool)
// notice that v0 needs to be restricted because of pow((v/v0),delta)
//###############################################################

double accACC(double s, double v, double vl, double a_lead,
	      const double beta[]){
  double delta=4.;
  double cool=0.5*(1.+tanh(beta[5]-1.)); // maps to values between 0 and 1
  double v0=(beta[0]<100) ? beta[0] : beta[0]+log(beta[0]/100);
  double  T=maxsmooth(0,beta[1], 0.01);

  double s0=maxsmooth(0,beta[2], 0.01); 
  double  a=maxsmooth(0,beta[3], 0.01); 
  double  b=maxsmooth(0,beta[4], 0.01); 

  double sstar=max(s0, s0+v*T+0.5*v*(v-vl)/sqrt(a*b));

  double z=sstar/max(s,0.01);
  double accEmpty=(v<=v0) ? a*(1- pow((v/v0),delta))
    : -b*(1- pow((v0/v),a*delta/b));
  double accPos=accEmpty*(1.- pow(z, min(2*a/accEmpty, 100.))  );
  double accInt=a*(1-z*z);

  double accIIDM=(v<v0) 
    ?  (z<1) ? accPos : accInt 
    :  (z<1) ? accEmpty : accInt+ accEmpty;

  double a_lead_loc=0; //!!

  double dvp=max(v-vl, 0.0);
  double v_lead = v-dvp;
  double denomCAH        =  v_lead*v_lead - 2 * s * a_lead_loc;

  double accCAH   = ( (v_lead*dvp  < - 2 * s * a_lead_loc) &&(denomCAH!=0))
    ? v*v*a_lead_loc/denomCAH
    : a_lead_loc - 0.5*dvp*dvp/max(s, 0.0001);
    
  // mix IIDM-CAH 

  double accACC_IIDM=(accIIDM>accCAH)
    ? accIIDM
    : (1-cool)*accIIDM + cool*( accCAH+b*tanh((accIIDM-accCAH)/b));


  if(isnan(accACC_IIDM)){
    cerr<<"s="<<s<<" v="<<v<<" vl="<<vl<<endl;
    cerr<<"v0="<<v0<<endl;
    cerr<<"T="<<T<<endl;
    cerr<<"s0="<<s0<<endl;
    cerr<<"a="<<a<<endl;
    cerr<<"b="<<b<<endl;
    cerr<<"cool="<<cool<<endl;

    cerr<<"accACC_IIDM="<<accACC_IIDM<<endl;
    exit(-1);
  }
  return max(accACC_IIDM, -BMAX);
 }


//###############################################################
// OVM/FVDM/LCM with triangular fundamental diagram (v0,T,s0,a,gamma)
//###############################################################

double accFVDM(double s, double v, double vl, const double beta[]){
  double v0=(beta[0]<100) ? beta[0] : beta[0]+log(beta[0]/100);
  double  T=maxsmooth(0,beta[1], 0.1);
  double s0=maxsmooth(0,beta[2], 0.1); 
  double  a=maxsmooth(0,beta[3], 0.1); //a=beta[3]; => big errors
  //a=minsmooth(40, a, 10);
  //a=min(40., a);
  double  gamma=beta[4]; // sensitivity to speed differences
  double vopt=max( min((s-s0)/T, v0), 0.);  // optimal velocity

  //!!NEW: no param restrictions as demo test. SSE should be as model10!

  if(FVDM_FULL_LINEAR){
    v0=beta[0];
    T=beta[1];
    s0=beta[2];
    a=beta[3];
    //vopt=(s-s0)/T; // full linear
    vopt=min(v0, max(0.,(s-s0)/T)); // with desired speed
  }

  return max(-BMAX, a*(vopt-v)/v0 - gamma*(v-vl));
}

double accOVM(double s, double v, double vl, const double beta[]){
  double betaFVDM[5];
  for (int k=0; k<4; k++){betaFVDM[k]=beta[k];} betaFVDM[4]=0;
  return accFVDM(s,v,vl,betaFVDM); // no BMAX restriction
}

// linear control model (LCM)
// beta[0]=-a*s0/(v0*T), beta{1]=a/(v0*T), beta[2]=-a/v0, beta[3]=gamma

double accLCM(double s, double v, double vl, const double beta[]){
  double v0=20; //!!! fixed
  double accLin=beta[0]+beta[1]*s+beta[2]*v+beta[3]*(vl-v);
  double acc=(v<v0) ? accLin : min(0.,accLin);
  return min(BMAX, max(-BMAX, acc));
}


//###############################################################
// Gipps (v0,T,s0,a,b)
//###############################################################

double accGipps(double s, double v, double vl, const double beta[]){
  
  double v0=beta[0];
  double T=maxsmooth(0,beta[1], 0.01);
  double s0=maxsmooth(0,beta[2], 0.01); 
  double  a=maxsmooth(0,beta[3], 0.01);
  double  b=maxsmooth(0,beta[4], 0.01);



  // simple: my simplified Gipps

  // double aFree=a;
  // double vSafe=-b*T+sqrt(b*b*T*T+vl*vl+2*b*max(s-s0,0.)); // simple!!!
  // double vGipps=min(min(v+a*T,v0), vSafe);

  // full Gipps: assumptions theta=0.5 and bl=b or bl!=b (see mic->Gipps.cpp)

  double bdivbl=1.0; // b/bl for full Gipps !!!
  double vFree=v+2.5*T*a*(1-v/v0)*sqrt(0.025+v/v0);
  double vSafe=-b*T+sqrt(b*b*T*T+vl*vl*bdivbl+2*b*max(s-s0,0.)-b*v*T);
  double vGipps=min(vFree, vSafe);


  //double vnew=v+DT/T*(vGipps-v);  // vNew is only part of vGipps if DT/T<1
  //return max(-BMAX, (vnew-v)/DT); //
  return max(-BMAX, (vGipps-v)/T); // this is equivalent to above 2 lines
}


//###############################################################
// ADAS v5 (v0,T,s0,a,c1)
// need to restrict vo since NaN for pow((v/v0),delta) and extreme v0
//###############################################################

double accADAS(double s, double v, double vl, const double beta[]){
  double v0=(beta[0]<100) ? beta[0] : beta[0]+log(beta[0]/100); 
  double   T=maxsmooth(0,beta[1], 0.01);
  double s0=maxsmooth(0,beta[2], 0.01); 
  double   a=maxsmooth(0,beta[3], 0.01);
  double  c1=beta[4];
  double eta=beta[5];

  double sloc=max(s,0.001);
  double vd=max(0., (s-s0)/T);
  double vopt=min(v0, vd);
  double aOVM=a/v0*(vopt-v);
  double dv=v-vl;
  return max(-BMAX,
	     min(aOVM,
		 a/v0*(vd-v) -2*c1*dv*exp(s0/sloc)*(1+s0*dv/(eta*sloc*sloc))));
}



// ####################################################
// generic local acceleration function of the car-following model
// ####################################################

double accMic(int choice_model, double s, double v, double vl, 
	   double al, const double beta[]){

  if(choice_model>10){ 
    cerr<<"Error: choice_model="<<choice_model
	<<" must be <= 10"<<endl; exit(-1);
  }

  return
    (choice_model<2)      ?   accIDM(s,v,vl,beta)
      : (choice_model==8) ?   accIDM(s,v,vl,beta)
      : (choice_model==9) ?   accIDM(s,v,vl,beta)
      : (choice_model==7) ?   accIDM(s,v,vl,beta)
      : (choice_model==2) ?   accACC(s,v,vl,al,beta)
      : (choice_model==3) ? accGipps(s,v,vl,beta)
      : (choice_model==4) ?   accOVM(s,v,vl,beta)
      : (choice_model==5) ? accFVDM(s,v,vl,beta)
      : (choice_model==6) ?  accADAS(s,v,vl,beta)
      : (choice_model==10) ?  accLCM(s,v,vl,beta)
      : 0;
}


/*###################################################################
INTERNAL CALCULATION OF NEW TARGETS AND LEADING SPEED
For robustness and to tackle all sorts of data, 
new leading vehicles (new targets) can and all leading speeds will
be claculated internally from the gaps and the own speed alone. 4 steps:


(1) calculate all time indices (multiples of data sampling period) 
were new leaders appear (active/passive lane changes, no leader at all
=> array newTarget
heuristically if deviation between constant-speed heuristics 2s[i-2]-s[i-1]
and actual s[i] is greater than a minimum deviation sdevmin AND a ballistic
deviation using accmax

control parameters: sdevmin, accmax


(2) calculate leading speed from s and v (multiples of data sampling period)
N[i]   N[i+1]    Vl[i]
n      n         v+0.5*(s[i+1]-s[i-1])/dtData
n      y         v+(s[i]-s[i-1])/dtData
y      n         v+(s[i+1]-s[i])/dtData
y      y         v

(3) At actual simulation (dtSim <=dtData!!)

(4) final reset @ data level (no feedback to sim)

Exanple (3),(4): dtSim=2/3*dtData
newTarget (nT) assumed to be {0,1,0,0,1,1,0}

iSim  iData  iDataBefore iDtata2Before nT(iDataBefore) data@SSE  result
1     0      0           0             0               -         -
2     1      0           0             0               1         reset
3     2      1           0             1               0         sim(dtSim)
4     2      2           1             0               0         -
5     3      2           2             0               0         sim(4*dtSim)
6     4      3           2             0               0         sim(5*dtSim)
7     4      4           3             0               0         -
8     5      4           4             0               1         reset
9     6      5           4             1               1         reset
10    6      6           5             1               1         -
11    7      6           6             1 or 0          0         sim(dtSim) or 2*dtSim
12    8      7           6             0               0         sim(2*dtsim) or 3*dtSim  => simReset if nT(iDataBefore && (iData2Before !=iDataBefore)
  ##################################################################*/


// ####################################################
// determine new targets based on the time series of the gaps:
// new target[i] in data point i true if
// (0) explicit change of pair
// (1) ballistic heuristics implies |accel|>accmax
// (2) or constant-speed heuristics leads to deviation > sdevmin
// ####################################################

void getNewTargetTimes(const double *sdata, const int *pairdata, int ndata, 
		       double dtData, bool *newTarget){
  double accmax=10; 
  double sdevmin=5; 
  newTarget[0]=true; // start with a new target

  for (int i=1; i<ndata; i++){

    // discriminate explicitely if pairinfo exists

    if(pairdata[i]>-999){ // pairinfo exists; pairdata may be <0 for TL etc
      newTarget[i]=(pairdata[i]!=pairdata[i-1]);
      if(DEBUG_TARGETS&&newTarget[i]){
	cout<<"getNewTargetTimes: use explicit target info: "
	    <<" new target at idata=i="<<i<<" t-tmin="<<i*dtData<<endl;
      }
    }

    // otherwise discriminate according to ballistic considerations
    // newTarget[0]=true=>access [i-2] for i>=2 ok
    
    else{ 
      double sConstSpeed
        =(newTarget[i-1]==false) ? 2*sdata[i-1]-sdata[i-2] : sdata[i-1];
      double sdeviation=fabs(sdata[i]-sConstSpeed);
      //newTarget[i]= (!newTarget[i-1])
      newTarget[i]= (sdeviation>sdevmin)
        && (sdeviation>0.5*accmax*dtData*dtData);

      if(DEBUG_TARGETS&&newTarget[i]){
        cout <<"getNewTargetTimes: new explicit target info=>use ballistic:"
	     <<"\n i="<<i<<" t-tmin="<<i*dtData
	     <<" !newTarget[i-1]="<<(!newTarget[i-1])
	     <<" sdata[i]="<<sdata[i]
	     <<" sdata[i-1]="<<sdata[i-1]
	     <<" sdata[i-2]="<<sdata[max(0,i-2)]
	     <<" sConstSpeed="<<sConstSpeed
	     <<" sdeviation="<<sdeviation
	     <<" sdevmin="<<sdevmin
	     <<" 0.5*accmax*dtData*dtData="<<0.5*accmax*dtData*dtData
	     <<" newTarget[i]="<<newTarget[i]
	     <<endl;
      }

    }
  }

  
}



// ####################################################
/* Function for local calibration 
   hata are n IDM accelerations for the exogeneous variables  
   s[i]=data[i],
   v[i]=data[Ndata+i],
   vL[i]=data[2*Ndata+i],
   aL[i]=data[3*Ndata+i],
   pair[i]=data[4*Ndata+i],
   i=0 ... n-1
   will be fitted to ydata = n accelerations
*/
// ####################################################

// I have modified void* data and void* adata to double* data and double* adata, 
// respectively, everywhere since I do not know how to use the pointer-to-void structure

//######################################


// estimated endogeneous variable = array of accelerations hatacc
void micLocalFunc(double *beta, double *hatacc, int Mparam, int ndata, double *data){

  // extract exogeneous and ctrl. data
  // data[5*ndata+i] and data[6*ndata+i] defined by simulation
  register int i;
  double sdata[ndata];
  double vdata[ndata];
  double vldata[ndata];
  double aldata[ndata];
  int pairdata[ndata];

  // sdata already guaranteed to be >=GAP_MIN in the data processing stage   

  for(i=0; i<ndata; i++){
    sdata[i]=data[i];
    vdata[i]=data[ndata+i];
    vldata[i]=data[2*ndata+i];
    aldata[i]=data[3*ndata+i];
    pairdata[i]=(int)(data[4*ndata+i]);
  }
  int choice_model=(int)(data[7*ndata]);


  // implement possible restraints

  double betaLoc[Mparam+1];
  double v0fix=get_v0fixed(vdata, ndata);
  relocateParameters(choice_model, Mparam, v0fix, IDM_B, beta, betaLoc);

  // do the calculation

  for(i=0; i<ndata; i++){
    hatacc[i]=accMic(choice_model,sdata[i],vdata[i],
		     vldata[i],aldata[i],betaLoc);
  } 


}//micLocalFunc






//#####################################################################
// Global calibration by simulation
// estimated endogeneous variable calType=1: array of gaps
//                                calType=2: array of ln(gaps)
// calType contained in data
// hats_or_lns is the output array of simulated s or ln s
//####################################################################

void micGlobalFunc(double *beta, double *hats_or_lns, 
		   int Mparam, int ndata, double *data){

  // (1) extract exogeneous and ctrl. data from container
  // data[5*ndata+i] and data[6*ndata+i] defined by simulation
  
  // MT 2022-11 dynamic memory allocation
  double* sdata=new double[ndata];
  double* vdata=new double[ndata];
  double* vldata=new double[ndata];
  double* aldata=new double[ndata];
  int* pairdata=new int[ndata];
  bool* newTarget=new bool[ndata];

  for(int i=0; i<ndata; i++){
    sdata[i]=data[i];
    vdata[i]=data[ndata+i];  
    vldata[i]=data[2*ndata+i];
    aldata[i]=data[3*ndata+i];
    pairdata[i]=(int)(data[4*ndata+i]);
  }
  //cout<<" in micGlobalFunc: sdata[0]="<<sdata[0]
  //   <<" sdata[ndata-1]="<<sdata[ndata-1]<<endl;

  int choice_model    =(int)(data[7*ndata]);
  double dtData       =data[7*ndata+1];
  int calType         =(int)(data[7*ndata+2]);

  //!! set outside 0=false,1=true
  bool writeEndogVars =(int)(data[7*ndata+3]+0.5);



  // (2) calculate instances of new targets: bool newTarget[ndata]

  getNewTargetTimes(sdata, pairdata, ndata, dtData, newTarget);

  // (3) implement possible restraints
  // (notice: for IDMdelay, beta[5]=delay only used in micGlobalFunc
  // but not in underlying accMic

  double betaLoc[Mparam+1];
  double v0fix=get_v0fixed(vdata, ndata);
  relocateParameters(choice_model, Mparam, v0fix, IDM_B, beta, betaLoc);


  // (4) define time step used for simulation 
  // (may be different from timestep dtData of data)

  const double dtSim=get_dtSim(choice_model);
    //2 instances of def dtSim. Here for global calibration


  // (5) initialize simulation
  //!! int(nSimReal)=723 even if nSimReal=724
  
  const double nSimReal=ndata*dtData/dtSim;
  const int nSim=int(nSimReal+1e-6);
  //cout<<"init! ndata="<<ndata
  //    <<" nSimReal="<<nSimReal<<" nSim="<<nSim<<endl;
  double tmax=(ndata-1)*dtData;
  double sSim[nSim]; // sSim: nSim points; hats_or_lns: ndata points
  double vSim[nSim];
  double aSim[nSim];

  sSim[0]=sdata[0];
  vSim[0]=vdata[0];

  double sSimTr=sSim[0]; // following 3 lines only for finite reaction time
  double vSimTr=vSim[0];
  double vlDataTr=vldata[0];


  // ========================================================
  // (6) micro-simulation !!
  // i=time after calc of acc and state update. Input @ step i-1
  // dtSim<=dtData
  // => INTERNAL CALCULATION OF NEW TARGETS ...
  // ======================================================

  
  for (int i=1; i<nSim; i++){
  
    // time control and consistency check

    double tSim=i*dtSim;
    double tSimLast=tSim-dtSim; // before sim step (def since often used)

    //!! watch out! Sometimes int(4./2.)=1 happens!
    // finite precision of doubles can here be a problem!
    // (no problem in intp/intpextp since if i=4./2.=1, [i=1] taken anyway)
    int iData2Before=max(0, (int)((tSim-2*dtSim)/dtData+1e-6));  // test
    int iDataBefore=(int)(tSimLast/dtData+1e-6);  // before sim step
    int iData=(int)(tSim/dtData+1e-6);            // tData=iData*dtSim<=tSim

    if((sSim[i-1]<0) || (vSim[i-1]<0)){
      cerr<<"micGlobalFunc: gap or v negative! i="<<i
	  <<"sSim[i-1]="<<sSim[i-1]<<" vSim[i-1]="<<vSim[i-1]
	  <<endl;
      vSim[i-1]=max(vSim[i-1], 0.);

      for (int k=0; k<Mparam; k++){
	cerr<<"betaLoc["<<k<<"]="<<betaLoc[k]<<endl;
      }
      //exit(-1);
    }
			
    // ========================================================
    // (6a) in micro-simulation time loop:
    // calculate state vars for local acceleration micromodel
    // (t is new time, state var at old time)
    // ========================================================

 
    double sLast=sSim[i-1];
    double vLast=vSim[i-1];
   
    double vlLast=intp(vldata,ndata,tSimLast,0,tmax);
    double alLast=intp(aldata,ndata,tSimLast,0,tmax);


    
    // determine past state if reaction delay
  
    if(choice_model==7){
      double dtDelay=maxsmooth(0,beta[5],0.1);
      double diDelayDouble=dtDelay/dtSim; // from beta[5]
      double tOld=min(tSimLast, max(0.,tSimLast-diDelayDouble*dtSim));
 
      // get delayed state variables 

      int choice=1;  // 0=with intpextp, 1=explicit (more robust here!)

      if(choice==0){  
 	int di=(int)(diDelayDouble)+1; //delayed models
        sLast=(i<di) ? sSim[0] : intpextp(sSim,i,tOld,0,tSimLast);
        vLast=(i<di) ? vSim[0] : intpextp(vSim,i,tOld,0,tSimLast);
        vlLast=intpextp(vldata,ndata,tOld,0,tmax);
      }
      
      else if(choice==1){ // delayed/memory
	int diDelay=(int)(diDelayDouble);
	int di=diDelay+1;
	int i1=min(i-1, max(0,i-1-di));
	int i2=min(i-1, max(0,i-di));
	double r=diDelayDouble-diDelay;
	sLast=r*sSim[i1]+(1-r)*sSim[i2];
	vLast=r*vSim[i1]+(1-r)*vSim[i2]; 
	vlLast=intpextp(vldata,ndata,tOld,0,tmax);

	// this branch only for choice_model==7 (with memory)
	if(false&&(i<8)){cout<<"sim model 7: choice=0, i="<<i
		    <<" dtDelay="<<dtDelay
		    <<" i1="<<i1<<" i2="<<i2<<" r="<<r
		    <<" sSim[i1]="<<sSim[i1]
		    <<" sSim[i2]="<<sSim[i2]
		    <<" tSimLas="<<tSimLast
		    <<" tOld="<<tOld
		    <<" sLast="<<sLast
		    <<" vlLast="<<vlLast
		    <<endl;
	}
      }
    } // choice_model==7 branch for model(s) with memory


    
    // ========================================================
    // (6b) in micro-simulation time loop: !! calculate acceleration
    // ========================================================

    aSim[i-1]=accMic(choice_model,sLast,vLast,vlLast,alLast,betaLoc);
 

 
    // ========================================================
    // (6c) in micro-simulation time loop:
    // update speeds and positions with ballistic integration
    // !! first-order update for leader's contribution
    // to gap sometimes slightly better than second-order (0.5*vlLast+vlNew)
    // update (new targets => first order anyway)
    // ========================================================
    
    vSim[i]=max(0., vSim[i-1]+aSim[i-1]*dtSim);

    if(SECOND_ORDER_VLEAD){ // tSim=new time; tSimLast=old
      double vlNew=intp(vldata,ndata,tSim,0,tmax); 
      sSim[i]=max(0.1, sSim[i-1]
		  + 0.5*dtSim*(vlLast+vlNew)
		- 0.5*dtSim*(vSim[i-1]+vSim[i]));
    }
    else{
      sSim[i]=max(0.1, sSim[i-1]  // 1st order
		  + 1.0*dtSim*vlLast
		  - 0.5*dtSim*(vSim[i-1]+vSim[i]));
    }

    if(false){
      //if((iData>1425)&&(iData<1455)){//!!!
      cout<<"iData="<<iData<<" i="<<i
	  <<" sSim[i]="<<sSim[i]
	  <<" sSim[i-1]="<<sSim[i-1]
	  <<" vSim[i]="<<vSim[i]
	  <<endl;
    }


    //!!
    if(false){
    //if(i<10){
      cout<<"after accMic: it="<<i
	  <<" tSim="<<tSim//<<" nSim="<<nSim
	  <<" sLast="<<sLast
	  <<" vLast="<<vLast
	  <<" vlLast="<<vlLast
	  <<" aSim[i-1]="<<aSim[i-1]
	//  <<"\n    vSim[i-1]="<<vSim[i-1]
	// <<" sSim[i-1]="<<sSim[i-1]
	//  <<" vSim[i-1]="<<vSim[i-1]
	//  <<" aSim[i-1]="<<aSim[i-1]
	  <<endl;
    }
    
    
    // i=iSim
    
    bool debug=false; 

    // ========================================================
    // (6d) in micro-simulation time loop:
    // reset gap s to data if new target is detected
    // (array newTarget is defined before the simulation)
    // before sim is OK because, later,
    // checked on data timesteps (here sim) against actual change
    // (but no feedback to sim=>simply "if(false){.. " won't work)
    // => INTERNAL CALCULATION OF NEW TARGETS ...
    // ========================================================

    //first-order better! (newTarget=^ first order by principle)
    bool newTargetDetected=false;
    for(int i=iDataBefore; i<=iData; i++){//!!!
      if(newTarget[i]){newTargetDetected=true;}
    }
    //bool newTargetDetected=(USE_NEW_DATA_TO_DETERMINIE_TARGETS)
    //  ? newTarget[iData] : newTarget[iDataBefore]; //!!!

   
    // reset follower gap if new leader
    // (!! in data, gaps restricted to sdata>=GAP_MIN)
    if(newTargetDetected){
    //if(newTarget[iDataBefore]&&(iDataBefore!=iData2Before)){ // 2nd order
      
      sSim[i]=sdata[iData];  // v not reset!!! (=> (6e)
                  // v not reset!!! (=> (6e)
      //vSim[i]=vdata[iData];  // v not reset!!! (=> (6e)
      //if(true){
      if(DEBUG_TARGETS){ //!!!
	cout<<"new target! iData="<<iData
	    <<" tSim (always start with 0)="<<tSim
	    <<": reset sSim[i]=sdata[iData]="<<sSim[i]
	    <<", unchanged vSim[i]="<<vSim[i]<<endl;
      }
	
    }

    // ========================================================
    // (6e) in micro-simulation time loop:

    // reset follower speed to the data if new follower changes
    // can only happen for concatenated trajectory pair for
    // pooled calibration. Always automatically detected by
    // violation of local consistency acc !=(vNew-vOld)/dt

    if((iData==0)||(fabs(vdata[iData]-vdata[max(0,iData-1)])>BMAX*dtData)){
      vSim[i]=vdata[iData];
    }
    
    // ========================================================


    // !! debug output (dtSim != dtData; vldata calculated previously:
    // => "calculate dependent quantities vldata and aldata"

    //if(debug){
    if(false&&(betaLoc[0]>6.8)&&(betaLoc[0]<7.0)){

      cout //<<"After calc acc and possible new target manipulations"
	   //<<" (calc. with input from step i-1, now all at step i)"<<endl
	   <<" i="<<i<<" iDataBefore="<<iDataBefore
	   <<" iData="<<iData
	   <<" newTarget[iDataBefore]="<<newTarget[iDataBefore]
	//<<" newTarget[iData]="<<newTarget[iData]
	//<<endl
	   <<" sLast="<<sLast
	   <<" vlLast="<<vlLast
	   <<" sSim[i]="<<sSim[i]
	   <<" sdata[iData]="<<sdata[iData]
	//<<" vLast="<<vLast
	//<<" vlLast="<<vlLast<<" acc="<<aSim[i-1]
	   <<endl;
    }

  }
  // end simulation loop 
  
   
  // micGlobalFunc: define output arrays to be used in dlevmar or objFun 
  // at data time instants (intp!)

  if(calType==1){ // SSE(s)
    for (int i=0; i<ndata; i++){
      int isim=int(i*dtData/dtSim+1e-6); //!! 1e-6 crucial
      //hats_or_lns[i]=(newTarget[i]) //!!!
      //	? sdata[i] : intp(sSim, nSim, i*dtData, 0, tmax);
      hats_or_lns[i]=intp(sSim, nSim, i*dtData, 0, tmax);//!!!

      if(false){
	//if((i>1425)&&(i<1455)){//!!!
	cout<<"on outpu for dlevmar: i="<<i
	    <<" hats_or_lns[i]="<<hats_or_lns[i]
	    <<" i*dtData="<<i*dtData
	    <<" nSim="<<nSim<<" tmax="<<tmax
	    <<" isim="<<isim
	    <<" sSim[isim]="<<sSim[isim]
	    <<endl;
      }

    }
  }


  else if(calType==2){// SSE(ln s)
    for (int i=0; i<ndata; i++){
      hats_or_lns[i]=(newTarget[i])
	? log(sdata[i])
	: log(intp(sSim, nSim, i*dtData, 0, tmax));
    }
  }

  // !!! micGlobalFunc: restrict unphysical parameters/parameter excursions
  // by soft box-constraints
  // via gap mismatch increase
  // since output is just array of simulated ln s values (hatlns),
  // not an objective function, just increase distance from hats from data
  // Do this only during calibration, not for final output for .out file
  // (writeEndogVars=true), otherwise huge gap differences in OVM since there
  // constraints are really effective (Levmar just runs away otherwise)
  
  if((choice_model<=5)&&(!writeEndogVars)){
    //if(((choice_model<=5)||(choice_model==10))&&(!writeEndogVars)){
    double v0extreme=500; // at this value, dist to data value doubled
    double Tmin=0.2;
    double ds0=1;
    double amin=0.3;// 0.5
    double aextreme=10000; // at this value, dist to data value doubled
    double bmin=0.3; // 0.5
    double bextreme=100; // at this value, dist to data value doubled
    int iv0=0;
    int iT=(choice_model==1)  ? 0 : 1;
    int is0=(choice_model==1) ? 1 : 2;
    int ia=(choice_model==1)  ? 2 : 3;
    int ib=(choice_model==1)  ? 3 : 4;

    for (int i=0; i<ndata; i++){
      double hats_or_lnsdata=(calType==1) ? sdata[i] : log(sdata[i]);


      // v0: only from above (in model 1 v0 is fixed)
      // soft restraints "from above" are more robust if applied for
      // *all* values. However, for physical values, the restraint is
      // not relevant (e.g. v0=25: distance increased by factor (1+1/400);
      // SSE contribution increased by factor (1+1/200)

      if(choice_model!=1){ 
	  hats_or_lns[i]
	    +=SQR(beta[iv0]/v0extreme)*(hats_or_lns[i]-hats_or_lnsdata);
	}

      // T: only from below

      if(beta[iT]<Tmin){
	hats_or_lns[i]
	  +=SQR((Tmin-beta[iT])/Tmin)*(hats_or_lns[i]-hats_or_lnsdata);
      }
      
      // s0: only for <0

      if(beta[is0]<0){ // soft-restrict IDM-like s0
	hats_or_lns[i] +=SQR(-beta[is0]/ds0)*(hats_or_lns[i]-hats_or_lnsdata);
      }

      // a: soft restraints both from below and above
      
      if(beta[ia]<amin){ // soft-restrict IDM-like a
	hats_or_lns[i]
	  +=(amin/max(beta[ia],0.001)-1)*(hats_or_lns[i]-hats_or_lnsdata);
      }

      hats_or_lns[i]+=SQR(beta[ia]/aextreme)*(hats_or_lns[i]-hats_or_lnsdata);

      
      if(beta[ib]<bmin){ // soft-restrict IDM-like b
	hats_or_lns[i]
	  +=(bmin/max(beta[ib],0.001)-1)*(hats_or_lns[i]-hats_or_lnsdata);
      }
      hats_or_lns[i]+=SQR(beta[ib]/bextreme)*(hats_or_lns[i]-hats_or_lnsdata);
    }
  }


  // micGlobalFunc: define simulated output to be used in output file 
  // for plotting and write them to the data container data
  // e.g., speed if the GOF is gap related (needed, mar 18)
  // writeEndogVars=(data[7*ndata+3]==1) set in main if objLandscape=true=1;
  // or in final output .out
  // also increase of gap mismatch deactivated which is
  // needed for actual calib

  if(writeEndogVars){ // micGlobalFunc
    for (int i=0; i<ndata; i++){
      data[5*ndata+i]=intp(vSim, nSim, i*dtData, 0, tmax); 
      data[6*ndata+i]=intp(aSim, nSim, i*dtData, 0, tmax); 
    }
  }

  delete[] sdata;
  delete[] vdata;
  delete[] vldata;
  delete[] aldata;
  delete[] pairdata;
  delete[] newTarget;
  
  if(false){
    cout <<"micGlobalFunc: choice_model="<<choice_model<<endl;
    cout <<" global function for levmar, local for calcul. obj fun"<<endl;
    for(int p=0; p<Mparam; p++){
      cout <<"Global params: p="<<p<<" beta[p]="<<beta[p]<<endl;
    }
    for(int p=0; p<5; p++){
      cout <<"Local params: p="<<p<<" betaLoc[p]="<<betaLoc[p]<<endl;
    }
    exit(0);
  }

} // end micGlobalFunc


//########################################################################
// Platoon calibration by simultaneous simulation of a platoon of followers
//
// number ntraj of trajectories extracted from data[7*ndata+4]
// number of data points  must be multiple of ntraj: ndata=ntraj*nTimes
// (notice: local var nSim=number of *simulated* time steps)
//
// each trajectory must begin and end at the same time, respectively
// trajectory itraj follows immediaely itraj-1; 
// no active or passive lane changes are allowed
//
// calType contained in data
// hats_or_lns is the output array of simulated s or ln s
//########################################################################

//TODO: (1) Data preparation and check 
//          (no target changes, equal tmin, tmax) in main program
//      (2) implement in main and pass ntraj=data[7*ndata+4] as 
//          command-line param

void micPlatoonFunc(double *beta, double *hats_or_lns, int Mparam, 
int ndata, double *data){


  // (1) extract exogeneous and ctrl. data from container
  // data[5*ndata+i] and data[6*ndata+i] defined by simulation


  double sdata[ndata]; // gaps of all ntraj trajectories
  double vdata[ndata];
  double vldata[ndata];
  double aldata[ndata];
  int pairdata[ndata];

  for(int i=0; i<ndata; i++){
    sdata[i]=data[i];
    vdata[i]=data[ndata+i];  
    vldata[i]=data[2*ndata+i];
    aldata[i]=data[3*ndata+i];
    pairdata[i]=(int)(data[4*ndata+i]);
  }
  int choice_model =(int)(data[7*ndata]);


  double dtData    =data[7*ndata+1];
  int calType      =(int)(data[7*ndata+2]);

  //!! micPlatoonFunc: set outside 0=false,1=true
  bool writeEndogVars=(int)(data[7*ndata+3]+0.5);

  int ntraj=(int)(data[7*ndata+4]+0.5); //first veh: traj 0; last veh: traj ntraj-1
  if(ndata%ntraj!=0){cerr<<" calibtraj, micPlatoonFunc: ndata="<<ndata
			  <<" is not a multiple of ntraj="<<ntraj<<endl;
    exit(-1);
  }
  int nTimes=ndata/ntraj; 

  // (2) no Target changes since not allowed (must be catched in main program)

  // (3) calculate possible parameter permutations and externally fixed values
  // (as in micGlobalFunc)

  double betaLoc[Mparam+1];
  double v0fix=get_v0fixed(vdata, ndata);
  relocateParameters(choice_model, Mparam, v0fix, IDM_B, beta, betaLoc);


  // (4) define time step used for simulation (may be different from timestep dtData of data)
  // (as in micGlobalFunc)


  const double dtSim=get_dtSim(choice_model);
  // 2 instances of def dtSim. Here for Platoon calibration

  

  // (5) initialize (no finite Tr)
  // sSim: nSim*ntraj points; 
  // hats_or_lns (=endog.argument): ndata=nTimes*ntraj points

  const int nSim=(int)((nTimes-1)*dtData/dtSim);
  double tmax=(nTimes-1)*dtData;
  double sSim[ntraj][nSim]; 
  double vSim[ntraj][nSim];
  double aSim[ntraj][nSim];

  for (int itraj=0; itraj<ntraj; itraj++){
    sSim[itraj][0]=sdata[itraj*nTimes];
    vSim[itraj][0]=vdata[itraj*nTimes];
  }

  // ==============
  // (6) micro-simulation platoon calibration
  // ==============

  for (int itraj=0; itraj<ntraj; itraj++){
    for (int i=1; i<nSim; i++){

      // times and consistency check

      double tSimLast=(i-1)*dtSim;
      double tSim=i*dtSim;

      if((sSim[itraj][i-1]<0) || (vSim[itraj][i-1]<0)){
        cerr<<"itraj="<<itraj<<" i="<<i
          <<" sSim[itraj][i-1]="<<sSim[itraj][i-1]
	  <<" or vSim[itraj][i-1]="<<vSim[itraj][i-1]
	  <<" negative" <<endl;
        for (int k=0; k<Mparam; k++){
	  cerr<<"betaLoc["<<k<<"]="<<betaLoc[k]<<endl;
	}
        //exit(-1);
      }

      // calculate exogeneous vars of leading veh

      double vl=(itraj==0) ? intp(vldata,ndata,tSim,0,tmax)
	: vSim[itraj-1][i];
      double vlold=(itraj==0) ? intp(vldata,ndata,tSimLast,0,tmax)
	: vSim[itraj-1][i-1];
      double alold=(itraj==0) ? intp(aldata,ndata,tSimLast,0,tmax)
	: aSim[itraj-1][i-1];

      // calculate acceleration

      aSim[itraj][i-1]=accMic(choice_model,sSim[itraj][i-1],
			      vSim[itraj][i-1],vlold,alold,betaLoc);

      // update speeds and positions with ballistic integration


      vSim[itraj][i]=max(0., vSim[itraj][i-1]+aSim[itraj][i-1]*dtSim);
      sSim[itraj][i]=max(0.1, sSim[itraj][i-1]
			 + 0.5*dtSim*(vl+vlold)
			 - 0.5*dtSim*(vSim[itraj][i]+vSim[itraj][i-1]));
    }
  } // end of platoon simulation



  // define output arrays to be used in dlevmar or objFun 
  // at data time instants (intp!)

  if(calType==1){ // hats
    for (int itraj=0; itraj<ntraj; itraj++){
      for (int i=0; i<nTimes; i++){
        hats_or_lns[nTimes*itraj+i]= intp(sSim[itraj], nSim, i*dtData, 0, tmax); 
      }
    }
  }
  else if(calType==2){ // ln(hats)
    for (int itraj=0; itraj<ntraj; itraj++){
      for (int i=0; i<nTimes; i++){
        hats_or_lns[nTimes*itraj+i]= log(intp(sSim[itraj], nSim, i*dtData, 0, tmax)); 
      }
    }
  }


  // define output arrays to be used in output file for plotting
  //  and write them to the data container data

  if(writeEndogVars){  // micPlatoonFunc 
    for (int itraj=0; itraj<ntraj; itraj++){
      for (int i=0; i<nTimes; i++){
        data[5*ndata+nTimes*itraj+i]=intp(vSim[itraj], nSim, i*dtData, 0, tmax);
        data[6*ndata+nTimes*itraj+i]=intp(aSim[itraj], nSim, i*dtData, 0, tmax);
      }
    }
  }


}// end micPlatoonFunc(




//#######################################################################
// Calculate objective function SSE (haty[i]-ydata[i])
// (the above functions calculate the arrays of haty, 
// i.e.,predicted accel, gaps, ln gaps)
// Notice: the levmar optimizer dlevmar does not need objFun but only
// the function func=> micLocalFunc, micGlobalFunc calculating the array of haty
//#######################################################################

double objFun(void (*func)(double *beta, double *haty, 
			   int Mparam, int ndata, double *data),
	      double *beta, 
	      int Mparam, int ndata, 
	      double *data){

  double haty[ndata]; 
  double ydata[ndata]; 

  double dtData=data[7*ndata+1];
  int calType=(int)(data[7*ndata+2]);

  // calculate array of ydata

  if(calType==0){ //ydata=adata=finite diff of vdata;vdata[i]=data[ndata+i]
    ydata[0] =(data[ndata+1]-data[ndata])/dtData;
    ydata[ndata-1] =(data[2*ndata-1]-data[2*ndata-2])/dtData;
    for(int i=1; i<ndata-1; i++){
      ydata[i] =0.5*(data[ndata+i+1]-data[ndata+i-1])/dtData;
    }
  }

  else if(calType==1){ //ydata=sdata, sdata[i]=data[i]
    for(int i=0; i<ndata; i++){ydata[i]=data[i];}
  }

  else if(calType==2){//ydata=ln(sdata)
    for(int i=0; i<ndata; i++){ydata[i]=log(data[i]);}
  }


  // calculate array of haty

  func(beta,haty,Mparam,ndata,data);


  // calculate SSE of (haty[i]-ydata[i]) and return it

  double obj=0;
  for(int i=0; i<ndata; i++){
    obj+=(haty[i]-ydata[i])*(haty[i]-ydata[i]);
    //if((haty[i]-ydata[i])*(haty[i]-ydata[i])>10){cout<<"i="<<i<<" ydata="<<ydata[i]<<" haty="<<haty[i]<<endl;}
  }
  return obj;
}




//######################################################################
int main(int argc, char* argv[]){
//######################################################################


  // ####################################################
  // input: File names
  // ####################################################

  if ( !( ((argc>=5) && (argc<=7)) || (argc>=9)) ){

    cerr <<"\n\nCalling sequence: calibTraj GOF modelnumber pairdata calcLandscape [dt_smooth] [typeID]\n"
	 <<"                  calibTraj GOF modelnumber pairdata <5 or more fixed params>\n"
         <<"                  (no calcLandscape; superfluous params ignored but too few params not caught,\n"
	 <<"                  for OVM use FVDM with last param=0)\n\n"
    <<"  GOF={0: SSE(a), 1: SSE(s, 2: SSE(ln s)}\n";
    cerr <<"  modelnumber={IDM,IDM_v0fix,ACC,GIP,OVM,FVDM,ADAS,IDMdelay,IDM_v0bfix,IDM_bfix,LCM}\n";
    cerr <<"  pairdata.FCdata including cols with headers t[s]  v[m/s]  s[m] and optionally headers pairs and/or id\n";
    cerr<<"       - all others kinematic vars (vl, al) derived internally\n";
    cerr<<"       - if no \"pairs\" or \"ID\" column, new leaders transitions are estimated from the data\n";
    cerr <<"  calcLandscape=1 if objective function landscape calulated\n\n";

    cerr<<"Examples: \n\n";
    cerr <<"  Local IDM calibration w/o landscape calc:"<<endl
	 <<"    calibTraj 0 0 Boschdata3 0"<<endl
	 <<"  Global ACC calibration lns with landscape calc and 1.5 s smoothing:"<<endl
	 <<"  Global IDM calibration s w/o landscape, 2 s smoothing, vehtype 2:"<<endl
	 <<"    calibTraj 1 0 ChennaiPairs 0 2 2"<<endl
	 <<"    [do not now now what typeselection does (2022-12-01)]"<<endl
         <<"  Just running the simulation with ACC with fixed parameters:"<<endl
	 <<"    calibTraj 1 2 NGSIM_veh50 15 1.2 3 2 2 1"<<endl
	 <<"  Running FVDM with fixed parameters:"
	 <<" calibTraj 1 5 NGSIM_veh50 15 1.2 3 2 0.2"<<endl
	 <<"  Running OVDM with fixed parameters:"
	 <<" calibTraj 1 5 NGSIM_veh50 15 1.2 3 2 0"<<endl;
   exit (-1);
  }


    
  char EFCname[256];
  char projName[256];
  char outName[256];

  int calType=atoi(argv[1]);
  int choice_model=atoi(argv[2]);
  sprintf(projName,"%s",argv[3]);
  sprintf(EFCname,"%s.FCdata",argv[3]);
  sprintf(outName,"%s.out",argv[3]);   // "projName" instead of argv[3]
                                       // compiler warning

  bool calcObjLandscape=(atoi(argv[4])==1); // overwritten for fixed params
  double dtSmooth=(argc==6)||(argc==7) ? atof(argv[5]) : 0;
  int typeSelected=(argc==7) ? atoi(argv[6]) : -1;
  bool simFixedBeta=(argc>=8);

  double betaFixed[argc-4]; // max of argv is argv[argc-1] as usual

  if(simFixedBeta){
    calcObjLandscape=false;
    cout <<"Calculating only trajectories with fixed values"<<endl;
    for (int m=0; m<argc-4; m++){betaFixed[m]=atof(argv[m+4]);}

    if(calType==0){
      cerr<<" Calculating of trajectories with fixed values only useful for"
	  <<" calType=1 (s) or calType=2 (ln s)"<<endl;
      exit(-1);
    }
  }


  const char* model0="IDM";
  const char* model1="IDM-v0fix";
  const char* model8="IDM-v0bfix";
  const char* model9="IDM-bfix"; // choice_model==9
  const char* model2="IIDM-ACC";
  const char* model3="Gipps";
  const char* model4="OVM";
  const char* model5="FVDM";
  const char* model10="LCM";
  const char* model6="ADAS";
  const char* model7="IDM-delay";
  const char* modstring[]={model0, model1, model2, model3, model4, 
			   model5, model6, model7, model8, model9, model10};


  // ####################################################
  // input: Chose model
  // ####################################################

  //{0=IDM(5 params), 1=IDM_v0fixed(4), 2=ACC(6), 3=Gipps(5), 4=OVM(4),
  // 5=FVDM(5), 6=ADAS(6), 7=IDMdelay, 8=IDM_v0bfixed, 9=IDM_bfixed}
 
  const int Mparam=
    (choice_model==0) ? 5 :
    (choice_model==1) ? 4 :
    (choice_model==8) ? 3 :
    (choice_model==9) ? 4 :
    (choice_model==2) ? 6 :
    (choice_model==3) ? 5 :
    (choice_model==4) ? 4 :
    (choice_model==5) ? 5 :
    (choice_model==10) ? 4 :
    (choice_model==6) ? 6 : 6;

  bool IDMv0constr=(choice_model==1)||(choice_model==8);
  bool lin=(choice_model==10);
  const char* pstring0=IDMv0constr ? "T" : (lin) ? "beta0" : "v0";
  const char* pstring1=IDMv0constr ? "s0": (lin) ? "beta1" :  "T";
  const char* pstring2=IDMv0constr ? "a" : (lin) ? "beta2" :  "s0";
  const char* pstring3=IDMv0constr ? "b" : (lin) ? "beta3" :  "a";
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

   // later possibly ctrl variables

  //double ctrlVals[256];
  //int nCtrl=0;
  //inout.get_array(ctrlName,nCtrl,ctrlVals);


  double opts[LM_OPTS_SZ]; // control parameters for dlevmar_[dif|der] (input)
  double info[LM_INFO_SZ]; // information on the calculation of results (output) 
  const int n_iter_max=1000; // maximum number of iterations



  // ####################################################
  // Input: exogeneous and endogeneous measurement data
  // ####################################################

  InOut inout;
  Statistics stat;

  // parse header to get data columns for t,s,v 
  // (the rest, e.g., vl, is calculated by 
  // internal and platoon consistency)

  string str_t="#t[s]";
  string str_t2="#t(s)";
  string str_t3="#time[s]";
  string str_s="s[m]";
  string str_s2="s(m)";
  string str_s3="gap[m]";
  string str_s4="dist(m)";
  string str_v="v[m/s]";
  string str_v2="v(m/s)";
  string str_v3="speed[m/s]";
  string str_v4="vx[m/s]";
  string str_vl="lead_vx";
  string str_vl2="vl[m/s]";
  string str_pair="pair";
  string str_pair2="pairs";
  string str_pair3="leadID";
  string str_type="type";
  string str_type2="type2";

  int col_t=inout.get_colnumber(EFCname,str_t);
  if(col_t<1){col_t=inout.get_colnumber(EFCname,str_t2);}
  if(col_t<1){col_t=inout.get_colnumber(EFCname,str_t3);}
  if(col_t<1){cerr<<"err: time header of "<<EFCname<<" neither "<<str_t
		  <<" nor "<<str_t2<<" nor "<<str_t3<<endl
		  <<" notice: no space between # and "
		  <<str_t<<" etc"; exit(-1);}

  int col_s=inout.get_colnumber(EFCname,str_s);
  if(col_s<1){col_s=inout.get_colnumber(EFCname,str_s2);}
  if(col_s<1){col_s=inout.get_colnumber(EFCname,str_s3);}
  if(col_s<1){col_s=inout.get_colnumber(EFCname,str_s4);}
  if(col_s<1){cerr<<"err: gap header of "<<EFCname<<" neither "<<str_s
		  <<" nor "<<str_s2<<" nor "<<str_s3<<endl; exit(-1);}


  int col_v=inout.get_colnumber(EFCname,str_v);
  if(col_v<1){col_v=inout.get_colnumber(EFCname,str_v2);}
  if(col_v<1){col_v=inout.get_colnumber(EFCname,str_v3);}
  if(col_v<1){col_v=inout.get_colnumber(EFCname,str_v4);}
  if(col_v<1){cerr<<"err: gap header of "<<EFCname<<" neither "<<str_v
		  <<" nor "<<str_v2<<" nor "<<str_v3<<endl; exit(-1);}

  int col_vl=inout.get_colnumber(EFCname,str_vl);
  if(col_vl<1){col_vl=inout.get_colnumber(EFCname,str_vl2);}
  bool vlInfoExists=(col_vl>0);
  if(vlInfoExists){
    cout <<"found column "<<col_vl<<" for the leading speed"<<endl;
  }
  else{
    cout <<"found no column for the leading speed;"
	 <<" will calculate it internally from speed and gap"<<endl;
  }

  if(!USE_VLINFO_IF_AVAILABLE){vlInfoExists=false;}
  
  int col_pair=inout.get_colnumber(EFCname,str_pair);
  if(col_pair<1){col_pair=inout.get_colnumber(EFCname,str_pair2);}
  if(col_pair<1){col_pair=inout.get_colnumber(EFCname,str_pair3);}
  bool pairInfoExists=(col_pair>0);

  bool typeInfoExists=false; int col_type=-1;
  if(typeSelected>=0){
    col_type=inout.get_colnumber(EFCname,str_type);
    if(col_type<1){col_type=inout.get_colnumber(EFCname,str_type2);}
    typeInfoExists=(col_type>=1);
    
    if(!typeInfoExists){
      cerr<<"err: type selection chosen but no suitable header "
	  <<str_type<<" nor "<<str_type2
	  <<" continuing w/o selection"<<endl;
    }
  }

  if(false){
    cout <<"col_t="<<col_t
	 <<" col_s="<<col_s
	 <<" col_v="<<col_v
	 <<" col_vl="<<col_vl
	 <<" col_pair="<<col_pair
	 <<" col_type="<<col_type
	 <<" vlInfoExists="<<vlInfoExists
	 <<" pairInfoExists="<<pairInfoExists
	 <<" typeInfoExists="<<typeInfoExists
	 <<endl;
    exit(0);
  }



  // Measure of Performance (MoP) is kind of used variable to calculated GoF
  // Goodness of Fit (GoF) a.k.a. objective function,
  // e.g, SSE, Theil's coefficient

  // MT 2020-02 dynamic memory allocation
 
  int ndata=inout.getNumberOfDataLines(EFCname);
  if(ndata<10){
    cerr<<"Error: the extended FC data contain too few data lines,"
	<<" ndata="<<ndata<<endl;
    exit(-1);
  }
  int ndatamax=2*ndata;             // take care of missing lines; worst case
  double* tdata=new double[ndatamax];
  double* sdata=new double[ndatamax];  // MoP1 global calibr
  double* vdata=new double[ndatamax];  // MoP2 global calibr
  int* pairdata=new int[ndatamax];    
  int* typedata=new int[ndatamax];    

  double* lnsdata=new double[ndatamax];// hybrid MoP1-GoF global calibr
  double* vldata=new double[ndatamax]; // leader's speed data (input)
  double* adata=new double[ndatamax];  // MoP local calibr
  double* aldata=new double[ndatamax]; // leader's accel for ACC model (input)
  double* Tdata=new double[ndatamax];  // timegap; needed to define
                                    // box boundaries of optimization

  
  inout.get_col (EFCname, col_v,  ndata, vdata);
  inout.get_col (EFCname, col_s,  ndata, sdata);
  inout.get_col (EFCname, col_t,  ndata, tdata);
  if(vlInfoExists){inout.get_col (EFCname, col_vl,  ndata, vldata);}
  else for(int i=0; i<ndata; i++){vldata[i]=-999;} // error value, overwritten
  if(pairInfoExists) inout.get_col (EFCname, col_pair,  ndata, pairdata);
  else for(int i=0; i<ndata; i++){pairdata[i]=-999;}
  if(typeInfoExists) inout.get_col (EFCname, col_type,  ndata, typedata);
  else for(int i=0; i<ndata; i++){typedata[i]=-1;}


  //#######################################################
  // extract vehicle type typeSelected if typeInfoExists
  // (in-array replacement)
  //#######################################################

  if(typeInfoExists){
    int isel=0;
    for(int i=0; i<ndata; i++){
      if(typedata[i]==typeSelected){
	tdata[isel]=tdata[i];
	sdata[isel]=sdata[i];
	vdata[isel]=vdata[i];
	pairdata[isel]=pairdata[i];
	typedata[isel]=typedata[i];
	isel++;
      }
    }
    ndata=isel;

    if(ndata<10){
      cerr<<"Error: only ndata="<<ndata
	  <<" < 10 data lines for selected vehicle type "
	  <<typeSelected<<". Select another type"<<endl;
      exit(-1);
    }

  }

  if(false){
    cout <<"after possible in-array-replacement: ndata="<<ndata<<endl;
    exit(0);
  }

  //#######################################################
  // do in-program smoothing and resampling if dtSmooth>0
  // (only useful if pairInfoExists or fixed leader-follower pair)
  //#######################################################
  
  cout <<"dtSmooth="<<dtSmooth<<endl;
  if(dtSmooth>1e-6){
    int nOut=0;              // true value to be determined
    double tdataOut[ndatamax];
    int pairOut[ndatamax];
    double vdataSmooth[ndatamax];
    double sdataSmooth[ndatamax];

    // determine natural sampling w/o missing lines

    double dt1=tdata[1]-tdata[0];
    double dt2=tdata[2]-tdata[1];
    double dtOut=dt1;
    if(fabs(dt1-dt2)>1e-6){
      dtOut=0.5;
      cerr<<"Warning: cannot determine natural sampling interval"<<endl
	  <<"  => resampling to dtOut="<<dtOut<<endl;
    }
    cout <<"dtOut="<< dtOut<<endl;

    // case w/o pair info: Smooth w/o interruptions

    if(!pairInfoExists){
      cerr<<"Warning: smoothing dtSmooth="<<dtSmooth
	  <<">0 but no vehicle pair info exists"<<endl
	  <<"  => makes only sense if no changing leaders or followers"<<endl
	  <<"  and increasing timestaps"
	  <<endl;

      double tmin=stat.getmin(tdata, ndata);
      double tmax=stat.getmax(tdata, ndata);
      int nOut=(int)((tmax-tmin)/dtOut)+1;
      stat.calcSymmetricEMA(tdata, vdata, ndata, dtSmooth,
			    dtOut, tdataOut, vdataSmooth);
      stat.calcSymmetricEMA(sdata, vdata, ndata, dtSmooth,
			    dtOut, tdataOut, vdataSmooth);
      ndata=nOut;
      for(int i=0; i<nOut; i++){
	tdata[i]=tdataOut[i];
	sdata[i]=sdataSmooth[i];
	vdata[i]=vdataSmooth[i];
	pairdata[i]=0;
      }

    }

    // case with pair info: Smooth with interruptions at pair changes
    // separate data arrays in parts with unchanged leader-follower pair

    else{  // pairInfoExists

      // (1) determine the episodes with same pair (leader-follower id)

      int npairs=0;  // to be determined
      int npair[ndata];  // pair ipair=0..npairs-1 has npair[ipair] lines 
      int ibeginpair[ndata];  // first original data line of pair ipair
      int oldpair=pairdata[0];
      int it=1; 
      int ipair=0;

      for(int i=1; i<ndata; i++){

	// extract the data lines for a fixed pair

	if((pairdata[i]==oldpair) &&(i<ndata-1)){
	  it++;
	}

	// pair finished=>store line info

	else{

	  if(i==ndata-1){// add last data point to last pair; otherw dropped
	    it++; i++;
	  }

	  npair[ipair]=it;
	  ibeginpair[ipair]=i-npair[ipair];
	  it=1;
	  ipair++;
	  npairs=ipair;
	  oldpair=pairdata[i];
	}
      }


      if(true){
	cout <<"main: pairInfoExists, after separating vehile pairs:"<<endl;
	for(int ipair=0; ipair<npairs; ipair++){
	  cout <<"ipair="<<ipair<<" npair[ipair]="<<npair[ipair]
	       <<" ibeginpair[ipair]="<<ibeginpair[ipair]
	       <<endl;
	}
	cout<<"ndata="<<ndata<<endl;
      }




      //(2) Do separate smoothing/resampling for each pair

      double tdataPair[ndatamax];
      double sdataPair[ndatamax];
      double vdataPair[ndatamax];

      double tdataTmp[ndatamax];
      double sdataTmp[ndatamax];
      double vdataTmp[ndatamax];
      int pairdataTmp[ndatamax];

      int iOutStart=0;
      for(int ipair=0; ipair<npairs; ipair++){
	int istart=ibeginpair[ipair];
        double tmin=tdata[istart];
        double tmax=tdata[istart+npair[ipair]-1];
        for(int it=0; it<npair[ipair]; it++){
	  tdataPair[it]=tdata[istart+it];
	  sdataPair[it]=sdata[istart+it];
	  vdataPair[it]=vdata[istart+it];
	}
	stat.calcSymmetricEMA(tdataPair, sdataPair, npair[ipair], dtSmooth,
			      dtOut, tdataOut, sdataSmooth);
	stat.calcSymmetricEMA(tdataPair, vdataPair, npair[ipair], dtSmooth,
			      dtOut, tdataOut, vdataSmooth);
 
	int nOut=(int)((tmax-tmin)/dtOut)+1;
        for(int it=0; it<nOut; it++){
	  tdataTmp[iOutStart+it]=tdataOut[it];
	  sdataTmp[iOutStart+it]=sdataSmooth[it];
	  vdataTmp[iOutStart+it]=vdataSmooth[it];
	  pairdataTmp[iOutStart+it]=pairdata[istart];
	}
	iOutStart += nOut;
      }
 
      //(3) finally copy back into data for calibration

      ndata=iOutStart;
      for(int i=0; i<ndata; i++){
	tdata[i]=tdataTmp[i];
	sdata[i]=sdataTmp[i];
	vdata[i]=vdataTmp[i];
	pairdata[i]=pairdataTmp[i];
      }


      if(true){
	cout<<"main: pairInfoExists, resampled+smoothed data:"<<endl;
	for(int i=0; i<ndata; i++){
	  if( (i<60)||(i>ndata-50)){
	    cout <<"i="<<i<<" tdata[i]="<<tdata[i]
		 <<"\t sdata[i]="<<sdata[i]
		 <<"\t  vdata[i]="<<vdata[i]
		 <<"\t  pairdata[i]="<<pairdata[i]
		 <<endl;
	  }
	}
	cout<<"resampled ndata="<<ndata<<endl;
      }

      //exit(0); 

    

    }

  }

  //#######################################################
  // start with actual calibration and/orr GoF function
  //#######################################################

  double smin=GAP_MIN; // needed for GOF=SSE(ln s)
  for (int i=0; i<ndata; i++){if (sdata[i]<smin) sdata[i]=smin;}
  for (int i=0; i<ndata; i++){vdata[i]=max(0., vdata[i]);} // v>=0
  for (int i=0; i<ndata; i++){lnsdata[i]=log(sdata[i]);}


 
  // calculate all relevant accelerations from primary quantities s and v
  // (for objective fun of local calibr and accACC model)

  double dtData=tdata[1]-tdata[0];
  if( !( (dtData>1e-6)||(dtData<5))){
      cerr<<"after reading data: in lnL: error: dtData not in range between 1e-6 and 5\n"
	<<"Something went wrong!"<<endl;
      exit(-1);
  }



  // determine target changes due to active or passive lane changes

  bool newTarget[ndata];  // {no new target, passive or active lane change}
  getNewTargetTimes(sdata,pairdata,ndata,dtData,newTarget);

 
  // calculate dependent quantities vldata and aldata by 
  // the conditions of local and platoon consistency
  // (we have at least 10 data points, so newTarget[1] etc defined)

  // => INTERNAL CALCULATION ... data level

  if(!vlInfoExists){
    vldata[0]=(newTarget[1]) ? vdata[0] : vdata[0]+(sdata[1]-sdata[0])/dtData;
    vldata[ndata-1]=(newTarget[ndata-1])
      ? vdata[ndata-1]
      : vdata[ndata-1]+ (sdata[ndata-1]-sdata[ndata-2])/dtData;
  }

  adata[0] =(vdata[1]-vdata[0])/dtData;
  adata[ndata-1] =(vdata[ndata-1]-vdata[ndata-2])/dtData;
  aldata[0]=0;
  aldata[ndata-1]=0;

  for (int i=1; i<ndata-1; i++){

    // interally calculate vLead if not available or not to be used
    
    if(!vlInfoExists){
      vldata[i]=(newTarget[i])
        ? vdata[i]+(sdata[i+1]-sdata[i])/dtData: (newTarget[i+1])
        ? vdata[i]+(sdata[i]-sdata[i-1])/dtData 
        : vdata[i]+0.5*(sdata[i+1]-sdata[i-1])/dtData;

      // take care of pathological spec cases => INTERNAL CALCULATION ... data

      if( newTarget[i] && (newTarget[i+1]||newTarget[i-1] )){
        vldata[i]=vdata[i];
      }
      vldata[i]=max(0.,vldata[i]);
    }
  
    adata[i]=0.5*(vdata[i+1]-vdata[i-1])/dtData;
    aldata[i]=((newTarget[i])||(newTarget[i+1]))
      ? 0 : 0.5*(vldata[i+1]-vldata[i-1])/dtData; // simplified
  }


 


  //##########################################################
  // wrap up everything in data container
  // exogeneous and endogeneous data vectors+ some ctrl data
  // {sdata[], vdata[], vldata[], aldata[], hatv[], hata[], 
  // choice_model,dtData,calType, writeEndogVars, ntrajIfPlatoonCalib}
  //##########################################################

  double* data = new double[7*ndata+4]; 

  for(int i=0; i<ndata; i++){
    data[i]=sdata[i];
    data[ndata+i]=vdata[i];
    data[2*ndata+i]=vldata[i];
    data[3*ndata+i]=aldata[i];
    data[4*ndata+i]=pairdata[i];
    // data[5*ndata+i]=hatv, data[6*ndata+i]=hata defined in simulation
  }

  // add ctrl data to the positions data[7*ndata+i] of the container

  data[7*ndata]=choice_model;
  data[7*ndata+1]=dtData;
  data[7*ndata+2]=calType;
  data[7*ndata+3]=false; // sets writeEndogVars and deactivates param control
  // during calibration, no writeEndogVars and gap mismatch increase
  // activated to control unphysical parameters driven by levmar




  // ####################################################
  // initialize parameter vector and define box boundaries
  // ####################################################

  double v0fixed=get_v0fixed(vdata,ndata);  // extra small function
  double bfixed=IDM_B;  // extra small function

  double vthr=5;
  for (int i=0; i<ndata; i++){Tdata[i]=sdata[i]/max(vdata[i],vthr);}


  beta[0]=stat.getmax(vdata,ndata);  //v0
  beta[1]=min(2., 0.7*stat.arithmeticMeans(Tdata,ndata));  //T
  beta[2]=stat.getmin(sdata,ndata);  //s0
  beta[3]=min(5.,0.5*stat.getmax(adata,ndata)); // a;  adata=accdata
  beta[4]=min(5.,max(1., -0.5*stat.getmin(adata,ndata)));  //b

  if(choice_model==4){beta[3]=15;} //OVDM: high acceleration
  if(choice_model==5){beta[4]=0.2;} //FVDM: speed diff sensitivity gamma
  if(choice_model==6){beta[4]=0.1;} //ADAS v5: c1
  beta[5]=(choice_model==2) ? 0.5 :
    (choice_model==6) ? 0.2 : 0.25; // cool (model2), eta (6), Treact(7)


  if(choice_model==10){
    beta[0]=-0.2; beta[1]=0.2; beta[2]=-0.2; beta[3]=0.2;
  }

  // box boundaries

  double betamin[Mparam];
  double betamax[Mparam];

  betamin[0]=0.2*beta[0]; betamax[0]=100;  // v0
  betamin[1]=0.2; betamax[1]=5;  // T
  betamin[2]=0.2; betamax[2]=10;  // s0
  betamin[3]=0.1; betamax[3]=100;//a
  betamin[4]=0.1; betamax[4]=100;//b
  betamin[5]=0; betamax[5]=3;//various


  if((choice_model==1)||(choice_model==8)){// IDM with v0fixed
    for(int k=1; k<5; k++){
      beta[k-1]=beta[k];
      betamin[k-1]=betamin[k];
      betamax[k-1]=betamax[k];
    }
  }

  if((choice_model==10) || (FVDM_FULL_LINEAR&&(choice_model==5))){
    for(int p=0; p<Mparam; p++){
      betamin[p]=-100;
      betamax[p]=100;
    }
  }

  if(simFixedBeta){ // no calibration, only simulation with prescribed beta
    for(int k=0; k<Mparam; k++){
      beta[k]=betaFixed[k];
    }
  }


  cout << "initial estimate of parameter vector:\n";
  for (int ip=0; ip<Mparam; ip++){
    cout <<"beta["<<ip<<"]="<<beta[ip]<<" "<<pstring[ip]<<endl;
  }

  // ####################################################
  // tests: just return single value(s) of SSE (objFun) for testing
  // ####################################################

  if(false){

    double v0= 10.91091; beta[0]=v0;
    double T=	0.69880; beta[1]=T;
    double s0=	3.10180; beta[2]=s0;
    double a=	2.71146; beta[3]=a;
    double b=	8.60354; beta[4]=b;
    double Tr=0.3;       beta[5]=Tr;// not yet used, test later

    double SSE=(calType==0)
      ? objFun(micLocalFunc,beta, Mparam, ndata, data)
      :objFun(micGlobalFunc,beta, Mparam, ndata, data);

    cout<<"\nIn Test: calType="<<calType<<" choice_model="<<choice_model<<endl;
    cout<<"beta=("<<beta[0]<<","<<beta[1]<<","<<beta[2]
	<<","<<beta[3]<<","<<beta[4]<<","<<beta[5]<<")"
	<<", SSE="<<SSE<<endl;
    exit(0);
  }

  // ####################################################
  // define optimization control and output parameters
  // ####################################################

  // passing to levmar NULL instead of opts reverts to defaults 
  opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
  opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used 


  double covar[Mparam*Mparam];
  int ret; // return signal of the dlevmar functions



  // ################################################################
  // The main optimization function dlevmar_dif

  // the box-constraint variant dlevmar_bc_dif is slower 
  // by a factor of 7 or 8 
  // => try unconstrained first since faster by about a factor of 5-10
  // => possibly use "soft" constraints by tanh functions in obj itself!
   // ################################################################

  /* ################################################################
  hint on structure: dlevmar_dif(func,beta,x,nbeta,nx,ctrlargs,addtlData)
  - the 1st arg is a function pointer to a function 
    func(beta,hatx,nbeta,nx,addtlData)
    defining on output an array *hatx of estim. values 
    out of the param array *beta
    addtlData are all the additional data needed to calculate func
  - the 2nd arg repeats beta but here used to control it =>hatbeta on output
  - the 3rd arg is the measured data x 
     with sum_i(x_i-hatx_i(beta))^2 to be minimized w/respect to beta,
  - arg 4,5,11 repeat args 3,4,5 of func
  - args 6-10 control the optimization 
  // ################################################################
  */

  if (!simFixedBeta){ 
    if(calType==0)
       ret= dlevmar_dif(micLocalFunc, beta, adata, Mparam, ndata, 
  		     n_iter_max, opts, info, NULL, covar, data);
    else if(calType==1)
      ret= dlevmar_dif(micGlobalFunc, beta, sdata, Mparam, ndata, 
  		     n_iter_max, opts, info, NULL, covar, data);
    else if(calType==2)
      ret= dlevmar_dif(micGlobalFunc, beta, lnsdata, Mparam, ndata, 
		     n_iter_max, opts, info, NULL, covar, data);


   // check if matrix is not singular 
   // (unfortunately, DOS by output diagnostic info[6], so need to determine
   // indirectly by investigating the correlation matrix

   bool isSingular=false;

   for (int im=1; im<Mparam; im++){
     for (int im2=0; im2<im; im2++){
       double corr=covar[Mparam*im+im2]
	 /sqrt(covar[Mparam*im+im]*covar[Mparam*im2+im2]);
       //cout <<"im="<<im<<" im2="<<im2<<" corr="<<corr<<endl;
       if( !(fabs(corr)<=1)){isSingular=true;}
     }
   }


   // optimization with box constraints dlevmar_bc_dif
   // if unconstrained optimization resulted in singular matrix

   if(isSingular){
     cerr<<"Warning: Singular matrix for trajectory pair "<<EFCname<<endl;
     cout <<"Singular matrix at unconstrained optimization,"
	  <<" arrived at following point:"<<endl;
     for (int m=0; m<Mparam; m++){
        cout <<"beta["<<m<<"]="<<beta[m]<<endl;
     }
     cout <<"with SSE="
	  <<objFun(micGlobalFunc,beta, Mparam, ndata, data)<<endl;
     cout <<" retry with my box-constraints!"<<endl;

     for (int m=0; m<Mparam; m++){
       beta[m]=min( max(beta[m], betamin[m]), betamax[m]);
       cout <<"starting with beta["<<m<<"]="<<beta[m]<<endl;
     }

     // only if isSingular

     if(calType==0)
       ret= dlevmar_bc_dif(micLocalFunc, beta, adata, Mparam, ndata, 
			  betamin, betamax,n_iter_max, opts,
			  info, NULL, covar, data);
     else if(calType==1)
       ret= dlevmar_bc_dif(micGlobalFunc, beta, sdata, Mparam, ndata, 
			  betamin, betamax,n_iter_max, opts, 
			  info, NULL, covar, data);
     else if(calType==2)
       ret= dlevmar_bc_dif(micGlobalFunc, beta, lnsdata, Mparam, ndata, 
			  betamin, betamax,n_iter_max, opts,
			  info, NULL, covar, data);
   }// isSingular
  }// !simFixedBeta


  // activates writeEndogVars, deactivates param excursion control
  data[7*ndata+3]=1; 
  double betaFinal[Mparam];
  for (int k=0; k<Mparam; k++){betaFinal[k]=beta[k];}

  double SSEinit=info[0];
  double SSEfinal=(!simFixedBeta)
    ? info[1]
    : objFun(micGlobalFunc,betaFinal, Mparam, ndata, data);


  // ####################################################
  // print result
  // ####################################################

  cout <<"\nResults of fitting "<<modstring[choice_model]<<" to "<<EFCname<<":\n";
  cout <<"==================================================\n\n";


  
  stringstream ss_header;
  ss_header<<"# produced by the call"<<argv[0];
  for(int icmd=1; icmd<argc; icmd++){
    ss_header<<" "<<argv[icmd];
  }
  
  ss_header<<"\n# Results of fitting "<< modstring[choice_model]
	   <<( (calType==0) ? " locally" : " globally")
	   <<" to "<<EFCname<<" with respect to"
	   <<( (calType==0) ? " accel" : (calType==1) ? " gaps": " ln gaps")
	   <<endl;


 
  string str_reason=(info[6]==1) ? "small gradient" :
    (info[6]==2) ? "small param change" :
    (info[6]==3) ? "itmax reached" :
    (info[6]==4) ? "singular matrix" :
    (info[6]==7) ? "NaN or Inf" :  "other reason";

  ss_header<<setprecision(5);  
  ss_header<<"\n# Levenberg-Marquardt returned in "
	   <<info[5]<<" iterations; reason: "<<str_reason
	   <<"\n# number of function evaluations: "<<info[7]
	   <<"\n# number of Jakobian evaluations: "<<info[8]
	   <<"\n\n# resulting_SSE:\t"<<SSEfinal;

  ss_header<<setprecision(2);  
 
  ss_header<<"\tavg_error="
	   <<((calType<2) ? sqrt(info[1]/ndata) : 100*sqrt(info[1]/ndata))
	   <<((calType==0) ? "m/s^2" : (calType==1) ? "m" : "%");

  ss_header<<setprecision(5);  
    
  ss_header<<"\n# initial_SSE\t"<<info[0]
	   <<endl;


  ss_header<<"\n##Best fit parameters:";
  if(IDMv0constr){
    ss_header<<"\nv0fixed=, e.g., 2*stat.getmax(vdata,ndata)="<<v0fixed;
  }
  
  if((choice_model==8)||(choice_model==9)){
    ss_header<<"\n#bfixed=, e.g., -stat.getmin(adata,ndata)="<<bfixed;
  }


  // header: fitted parameters and stdev

  ss_header<<setprecision(3)<<fixed;

  double stddev[Mparam];
  for (int im=0; im<Mparam; im++){
    stddev[im]=sqrt(covar[Mparam*im+im]);
    ss_header<<"\n#\t"<<pstring[im]<<"="<<betaFinal[im]
	     <<"\tsig_"<<pstring[im]<<"="<<stddev[im];
  }


  // correlation matrix

  if(!simFixedBeta){
    ss_header<<setprecision(3)<<fixed;
    ss_header<<"\n\n##Parameter corr Matrix:";
    for (int im=0; im<Mparam; im++){
     ss_header<<"\n# ";
     for (int im2=0; im2<Mparam; im2++){
       ss_header<<"r"<<im<<im2<<"="
		<<covar[Mparam*im+im2]/(stddev[im]*stddev[im2])<<"\t";
     }
   }
  }


  // get local accelerations/global trajectories for estimated parameters

  double hataFinal[ndata];
  double hats_or_lnsFinal[ndata];
  double hatvFinal[ndata];

  // data[7*ndata+3]=writeEndogVars=1=true: output final results
  // data[7*ndata+3]=writeEndogVars=0=false:
  // no writing and simultaneously increase of gap mismatch activated
  // to prevent param excursions
  
  data[7*ndata+3]=1; // 0=false, 1=true

  if(calType==0){
    micLocalFunc(betaFinal, hataFinal, Mparam, ndata, data);
  }
  else{
    micGlobalFunc(betaFinal, hats_or_lnsFinal, Mparam, ndata, data);
    for (int i=0; i<ndata; i++){ // micLocalFunc
      hatvFinal[i]=data[5*ndata+i];
      hataFinal[i]=data[6*ndata+i];
    }
  }


  //############################################
  // Write result to file
  //############################################

  cout <<ss_header.str()<<endl;

  cout <<"\nwriting results to "<<outName<<endl;

  char header[4096];
  ss_header<<"\n\n";
  
  if(calType==0){
    ss_header<<"\n#time\t\tsdata\tvdata\tvldata\tadata\thata";
    //const string str_header=ss_header.str(); // NOT directly .str.c_str()
    sprintf(header,"%s",ss_header.str().c_str()); // test
    inout.write_array(outName, ndata, tdata, sdata, vdata, vldata, 
		      adata, hataFinal, header);
  }

  else if(calType==1){
    ss_header<<"\n#time\t\tsdata\tvdata\tvlderiv\thats\thatv\thata";
    sprintf(header,"%s",ss_header.str().c_str());
    inout.write_array(outName, ndata, tdata, sdata, vdata, vldata,
		      hats_or_lnsFinal, hatvFinal, hataFinal,
		      header);
  }
  else if(calType==2){
    ss_header<<"\n#time\t\tsdata\tvdata\tvldata\tlnsdata\thatlns\thatv\thata";
    sprintf(header,"%s",ss_header.str().c_str());
    inout.write_array(outName, ndata, tdata, sdata, vdata, vldata, lnsdata,
		      hats_or_lnsFinal, hatvFinal, hataFinal,
		      header);
  }



  //############################################
  // Write objective fun landscape=log-likelihood around the  maximum
  // write min(SSE, SSEmax) because of plotting comfort
  //############################################

  double SSEmax=10*SSEfinal; 
  data[7*ndata+3]=0; // =0 deactivates writing of endogeneous var in objFun
                     // and activates gap mismatch for param excursion
                     // control. Still hardly mismatch at plotted vals
                     // and faster -> OK

  if(calcObjLandscape){

  double dtCorr=1; // assume calculated variances erroneous for dt<dtCorr
  double w_stddev=4*max(1.,sqrt(dtCorr/dtData)); // half-width of scanning range in (true) stddev
  int nout=31; // number of grid elements in either direction
  double v0minLimit=0.5*stat.getmax(vdata,ndata);
  double v0maxLimit=50;
  double TminLimit=-0.5;
  double TmaxLimit=3.0;
  double s0minLimit=-1.0;
  double s0maxLimit=8;
  double aminLimit=0.05; // not 0.5*stat.getmax(adata,ndata) bec of fluct!;
  double amaxLimit=(choice_model==4) ? 60 : 10;
  double bminLimit=0.05;
  double bmaxLimit=10;
  char objFunFilename[1024];
  char titleString[1024];
  double objFunData[inout.NYMAX][inout.NYMAX];

  double dbeta[Mparam];
  double betamin[Mparam];
  double betamax[Mparam];

  // check for undefined or too small stddev; then revert to standard stddev
  // for the IDM and Gipps variants
  // model={IDM,IDM_v0fix,ACC,GIP,OVM,FVDM,ADAS,
  //        IDMdelay,IDM_v0bfix,IDM_bfix,LCM}
  if((choice_model<4) || (choice_model>6)){
    if(!(stddev[0]>0.15)){stddev[0]=0.15;}
    if(!(stddev[1]>0.1)){stddev[1]=0.1;}
    if(!(stddev[2]>0.1)){stddev[2]=0.1;}
    if(!(stddev[3]>0.1)){stddev[3]=0.1;}
    if(!(stddev[4]>0.2)){stddev[4]=0.2;}
    if(!(stddev[5]>0.2)){stddev[5]=0.2;} // not relevant for Gipps
  }


  
  for (int j=0; j<Mparam; j++){
    betamin[j]=betaFinal[j]-w_stddev*stddev[j];
    betamax[j]=betaFinal[j]+w_stddev*stddev[j];
  }


  // !! IDM special (treatment of too small or undef stddev later)
  // (if s0<0 because of data, often var s0 approx 0, therefoe max(min()) )

  if(choice_model<10){
  betamin[0]=max(betamin[0],v0minLimit);
  betamax[0]=min(betamax[0],v0maxLimit);
  betamin[1]=max(betamin[1],TminLimit);
  betamax[1]=min(betamax[1],TmaxLimit);
  betamin[2]=max(betamin[2],s0minLimit);
  betamax[2]=min(betamax[2],s0maxLimit);
  betamin[3]=max(betamin[3],aminLimit);
  betamax[3]=min(betamax[3],amaxLimit);
  betamin[4]=max(betamin[4],bminLimit);
  betamax[4]=min(betamax[4],bmaxLimit);
  }
  

  if(choice_model==4){betamin[3]=4; betamax[3]=40;}
  if(choice_model==5){betamin[3]=1; betamax[3]=12;}
  if(choice_model==7){betamin[5]=0; betamax[5]=1.1;}

  // linear control model (LCM) +free restr dv/dt<=0 if v>=fixed v0=20
  // dv/dt=a/v0*((s-s0)/T-v)+gamma*(vl-v)
  //      =beta0*beta1*s+beta2*v+beta3*(vl-v))
  // beta[0]=-a*s0/(v0*T), beta[1]=a/(v0*T), beta[2]=-a/v0, beta[3]=gamma


  cout <<"LCM: betamin[0]="<<betamin[0]<<" beta[0]="<<beta[0]<<endl;
  if(false){
    //if(choice_model==10){
    betamin[0]=-1.5; betamax[0]=0.;
    betamin[1]=0.; betamax[1]=0.5;
    betamin[2]=-0.8; betamax[2]=0;
    betamin[3]=0; betamax[3]=1.5;
  }
  
  
  for (int j=0; j<Mparam; j++){
    dbeta[j]=(betamax[j]-betamin[j])/(nout-1);
    cout <<"j="<<j<<" beta[j]="<<beta[j]
	 <<" betamin[j]="<<betamin[j]<<" betamax[j]="<<betamax[j]<<endl;
  }



  // make (Mparam-1)*(Mparam-2) data sets and files

  for (int ibeta=0; ibeta<Mparam-1; ibeta++){
    for (int jbeta=ibeta+1; jbeta<Mparam; jbeta++){

      // generate the data set for a given beta combination

      // revert non-used dimensions
      for (int k=0; k<Mparam; k++){beta[k]=betaFinal[k];}
      for (int i=0; i<nout; i++){
	for (int j=0; j<nout; j++){

	  beta[ibeta]=betamin[ibeta]+i*dbeta[ibeta];
	  beta[jbeta]=betamin[jbeta]+j*dbeta[jbeta];
	  if(false){
	    cout <<"ibeta="<<ibeta<<" jbeta="<<jbeta<<" nout="<<nout
	         <<" SSEmax="<<SSEmax<<" ndata="<<ndata<<endl;
	    for(int k=0; k<Mparam; k++){
	      cout<<"  beta["<<k<<"]="<<beta[k]<<endl;
	    }
	  }
	  
	  objFunData[i][j]=(calType==0)
	    ? min(objFun(micLocalFunc,beta, Mparam, ndata, data), SSEmax)
	    : min(objFun(micGlobalFunc,beta, Mparam, ndata, data), SSEmax);
	}
      }


      // write the file for this combination

      sprintf(objFunFilename,"%s.beta%i_beta%i",projName, ibeta, jbeta);

      stringstream ss_header1;
      ss_header1<<"#Objective function for beta"<<ibeta
		<<" and beta"<<jbeta
		<<"\n#Base values: ";
      for(int k=0; k<Mparam; k++){
	ss_header1<<"beta"<<k<<"="<<betaFinal[k];
      }
      ss_header1<<"\n#Min Obj function: "<<info[1]
		<<"\n#beta"<<ibeta<<"\t\tbeta"<<jbeta<<"\t\tSSE";

      sprintf(header,"%s",ss_header1.str().c_str());

      inout.write_array2d(objFunFilename, betamin[ibeta],betamax[ibeta],nout,
			  betamin[jbeta],betamax[jbeta],nout,
			  objFunData,header);
    }
  }
  }

  // MT 2020-02 dynamic memory allocation

  delete[] data;
  delete[] tdata;
  delete[] sdata;   
  delete[] lnsdata; 
  delete[] vdata;   
  delete[] vldata;  
  delete[] adata;   
  delete[] aldata;  
  delete[] Tdata;  
  delete[] pairdata;  
  delete[] typedata;  

  cout <<" pairInfoExists="<<pairInfoExists<<" ndata="<<ndata<<endl;
  exit(0);
}

