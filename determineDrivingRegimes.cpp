

/*
(2023-10-31, answer to a question of Venkat)

Hi Venkat,

from my program identifying the regimes, I obtain
double vc=0.8*v0data;
double Tc=1.5*Tdata;
double ac=0.2*adata;
double bc=0.2*bdata;
double sc=1.5*s0data;
Here, v0data is the maximum speed in the data, adata and bdata are the maximum observed accelerations and decelerations, s0data the minimum gap, and Tdata is either a plausible value (e.g., 1.5 s) or the arithmetic average of sdata/vdata over all points with (v<vc) AND (a in [-bc,ac])
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
// c++ 
#include <iostream>
using namespace std;


#include "InOut.h"
#include "Math.h"
#include "Statistics.h"
#include "RandomUtils.h"


//######################################
// user defined global settings
//######################################

const int Ndata=10000; // n <= Ndata measurements 


//#############################################

int main(int argc, char* argv[]){

  // ####################################################
  // input
  // ####################################################


  if (argc!=6){
    cout <<"\nCalling sequence: determineDrivingRegimes <FCname> <v0min> <Tmax> <a> <b>\n"; 
    cout <<"Input: \n";
    cout <<" FC data file(s) in format as the .Bosch<i> files:\n";
    cout <<" t(s)   vL(m/s)    v(m/s)  a(SI)   s(m)    dv(m/s) ...\n";
    cout <<" v0min=minimal desired speed (to take care of purely follow data)\n";
    cout <<" Tmax=max. desired time headway (to take care of completely free data)\n";
    cout <<" a=desired free-flow starting acceleration\n";
    cout <<" b=desired deceleration\n";
    cout <<"Output: fraction of regimes (self-explaining) on standard output \n";
    cout <<" and the regime numbers are written (together with tdata) in file FCname.regimes\n";
    cout <<"\nExample:\n";
    cout <<"determineDrivingRegimes  Boschdata3 15 2 1 2\n";
    exit (-1);
  }

  char FCname[256];
  char outname[256];
  sprintf(FCname,"%s",argv[1]);
  sprintf(outname,"%s.regimes",argv[1]);

  double v0min=atof(argv[2]);
  double Tmax=atof(argv[3]);
  double a=atof(argv[4]);
  double b=atof(argv[5]);

  InOut inout;
  int ndata;
  double tdata[Ndata];
  double sdata[Ndata];
  double vdata[Ndata];
  double adata[Ndata];
  double dvdata[Ndata];

  inout.get_col(FCname, 1, ndata, tdata);
  inout.get_col(FCname, 3, ndata, vdata);
  inout.get_col(FCname, 4, ndata, adata);
  inout.get_col(FCname, 5, ndata, sdata);
  inout.get_col(FCname, 6, ndata, dvdata);


  // ####################################################
  // analysis
  // ####################################################

  double vthr=0.1;
  double s0max=3;
  double Tdata[Ndata];
  double TdataFollow[Ndata];

  int iFollow=0;
  for (int i=0; i<ndata; i++){
    Tdata[i]=sdata[i]/max(vdata[i],vthr);
    if(Tdata[i]<=Tmax) {TdataFollow[iFollow]=Tdata[i]; iFollow++;}
  }
  double ndataFollow=iFollow;
  

  Statistics stat;
  double v0=max(v0min,stat.getmax(vdata,ndata));
  double T=(ndataFollow>0) ? stat.arithmeticMeans(TdataFollow,ndataFollow):Tmax;
  double s0=min(s0max, stat.getmin(sdata,ndata)); // a and b already from cmd-line
  //a=max(a, stat.getmax(adata,ndata)); // override input
  //a=stat.getmax(adata,ndata); // override input
  //b=-stat.getmin(adata,ndata);

  double regimeIndices[Ndata];


  double vc=0.8*v0;
  double Tc=1.5*T;
  double ac=0.2*a;
  double bc=0.2*b;
  double sc=1.5*s0;
  double rc=0.1; // inverse time-of collision (collision rate)

  for (int i=0; i<ndata; i++){
    double bkin=(dvdata[i]<=0) ? 0 : 0.5*dvdata[i]*dvdata[i]/sdata[i];
    double r=dvdata[i]/sdata[i];


    if(sdata[i]<sc){
      regimeIndices[i]=(vdata[i]<vc) ? 2 : 7;
    }

    else{

      if((vdata[i]>vc)&&(Tdata[i]>Tc)){
	if(adata[i]>ac) regimeIndices[i]=3;
	else if(adata[i]>=-bc) regimeIndices[i]=0;//nCruise++;
	else regimeIndices[i]=(bkin<bc) ? 6 : 4;
      }

      else if((vdata[i]>vc)&&(Tdata[i]<=Tc)){
 	if(adata[i]>ac) regimeIndices[i]=7; //nAggressive++;
	else if(adata[i]>=-bc) regimeIndices[i]=1; //nSteadyFollow++;
	else regimeIndices[i]=(bkin<bc) ? 5 : 4;
      }
 
      else if((vdata[i]<=vc)&&(Tdata[i]>Tc)){
 	if(adata[i]>ac) regimeIndices[i]=3;//nAccFree++;
	else if(adata[i]>=-bc) regimeIndices[i]=5; //nOsc++;
	else regimeIndices[i]=(bkin<bc) ? 6 : 4;
      }

      else{
 	if(adata[i]>ac) regimeIndices[i]=5; //nOsc++;
	else if(adata[i]>=-bc) regimeIndices[i]=1; //nSteadyFollow++;
	else regimeIndices[i]=(bkin<bc) ? 5 : 4;
      }
    }

    // override by oscillat regime if too high inverse time of collisions
    if( (fabs(r)>rc) &&(adata[i]<=ac) && (adata[i]>=-bc)){
      if(false){
	cout <<"t="<<tdata[i]<<"dv/dt small and r="<<r
	     <<" has abs value greater than "<<rc<<endl;
      }
      regimeIndices[i]=5;
    }
  }

  int nCruise=0;   // v0
  int nSteadyFollow=0; // T
  int nStanding=0; // s0
  int nAccFree=0; //a
  int nAppr=0; //b
  int nOsc=0; //following in non-steady-state (a,b)
  int nDefensive=0; // inconsistent-defensive
  int nAggressive=0; // inconsistent-aggressive

  for (int i=0; i<ndata; i++){
    if(regimeIndices[i]==0) nCruise++;
    else if(regimeIndices[i]==1) nSteadyFollow++;
    else if(regimeIndices[i]==2) nStanding++;
    else if(regimeIndices[i]==3) nAccFree++;
    else if(regimeIndices[i]==4) nAppr++;
    else if(regimeIndices[i]==5) nOsc++;
    else if(regimeIndices[i]==6) nDefensive++;
    else nAggressive++;
  }

  // ####################################################
  // output
  // ####################################################

  cout <<"\nAnalysis for "<<FCname<<":\n\n";

  cout <<"Estimated driving properties:\n";
  cout <<" v0="<<v0<<" vc="<<vc<<"\n T="<<T<<" Tc="<<Tc<<endl;
  cout <<" s0="<<s0<<" sc="<<sc<<endl;
  cout <<" a="<<a<<" ac="<<ac<<endl;
  cout <<" b="<<b<<" bc="<<bc<<endl;

  cout <<"\nFraction of driving regimes:\n\n";
  cout <<"Cruising:                   n = "<<nCruise<<"\t("<<100*nCruise/ndata<<"%)"<<endl;
  cout <<"Steady-state following:     n = "<<nSteadyFollow<<"\t("<<100*nSteadyFollow/ndata<<"%)"<<endl;
  cout <<"Standing:                   n = "<<nStanding<<"\t("<<100*nStanding/ndata<<"%)"<<endl;
  cout <<"Free acceleration:          n = "<<nAccFree<<"\t("<<100*nAccFree/ndata<<"%)"<<endl;
  cout <<"Approaching:                n = "<<nAppr<<"\t("<<100*nAppr/ndata<<"%)"<<endl;
  cout <<"Non-steady-state following: n = "<<nOsc<<"\t("<<100*nOsc/ndata<<"%)"<<endl;
  cout <<"Inconsistent-aggressive:    n = "<<nAggressive<<"\t("<<100*nAggressive/ndata<<"%)"<<endl;
  cout <<"Inconsistent-defensive:     n = "<<nDefensive<<"\t("<<100*nDefensive/ndata<<"%)"<<endl;

  char timeString[1024];
  char sString[256];
  char vString[256];
  char regimeString[256];

  sprintf(timeString,"produced by determineDrivingRegimes %s %s %s %s %s\n#regimeIndex: 0=cruising,1=steady following,2=standing,3=free acc,4=approaching,5=osc following,6=inconsistent-aggress,7=inconsistent defens)\n#time[s]",FCname, argv[2], argv[3], argv[4], argv[5]);
  
  sprintf(sString, "gap[m]\t");
  sprintf(vString, "v[m/s]\t");
  sprintf(regimeString, "regimeIndex");
  cout <<"\nwriting regimes into file "<<outname<<" ...\n";
  inout.write_array(outname, ndata, tdata, sdata, vdata,regimeIndices,timeString,sString,vString,regimeString);

  exit(0);
}
