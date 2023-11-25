
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


  if (argc!=3){
    cout <<"\nCalling sequence: makeEFCfromFCdata <FCname>  <EFCname>\n"; 
    cout <<"Input: \n";
    cout <<" FC data file(s) in format as the .car<i> fies of mic:\n";
    cout <<" t(s)   x(m)    v(m/s)  a(SI)   s(m)    dv(m/s) ...\n";
    cout <<"Output: \n";
    cout <<" EFC data in format in Bosch formata:\n";
    cout <<" t(s)   vL(m/s) v(m/s)  a(SI)   s(m)    dv\n";
    cout <<"Example:\n";
    cout <<"makeEFCfromFCdata calibrProfile1.car2 calibrProfile1.EFC2_1\n";
    exit (-1);
  }

  char FCname[256];
  // char FCnameLeader[256];
  char EFCname[256];
  char EFCtitleString[1024];

  sprintf(FCname,"%s",argv[1]);
  cout <<"FCname="<<FCname<<endl;
  sprintf(EFCname,"%s",argv[2]);
  cout <<"EFCname="<<EFCname<<endl;
  sprintf(EFCtitleString,"Produced by makeEFCfromFCdata %s %s\n#this output file is third argument\n#time(s)",
	  FCname, EFCname);
  cout <<"EFCtitleString="<<EFCtitleString<<endl;

  InOut inout;
  int ndata;
  double tdata[Ndata];
  double sdata[Ndata];
  double vdata[Ndata];
  double adata[Ndata];
  double dvdata[Ndata];
  double vldata[Ndata];
  inout.get_col(FCname, 1, ndata, tdata);
  inout.get_col(FCname, 3, ndata, vdata);
  inout.get_col(FCname, 4, ndata, adata);
  inout.get_col(FCname, 5, ndata, sdata);
  inout.get_col(FCname, 6, ndata, dvdata);

  // ####################################################
  // make transform and output
  // ####################################################


  for (int i=0; i<ndata; i++){
    vldata[i]=vdata[i]-dvdata[i];
  }
  inout.write_array(EFCname, ndata, tdata, vldata, vdata, adata, sdata, dvdata,
		    EFCtitleString, "vldata", "vdata", "adata", "sdata", "dvdata");

  cout <<"\nwrote results to "<<EFCname<<endl;
  exit(0);
}
