#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// c++ 
#include <vector>
#include <iostream>
#include <sstream>
//#include <fstream>
//#include <istream>

using namespace std;


// own

#define SQR(x) ((x)*(x))
#include "Statistics.h"

//##################################################
// (feb17) Densities, distribution functions and inverse (quantile) functions
// of standard distributions =>
// StatisticalFunctions.cpp, .h
//##################################################


////////////////////////////////////////////////////////////////////////
/// general-purpose statistical operation
////////////////////////////////////////////////////////////////////////

double Statistics::sum(const double data[], int ndata) const{
  double sum=0;
  for (int i=0; i<ndata; i++){sum+=data[i];}
  return(sum);
}

int Statistics::sum(const int data[], int ndata) const{
  int sum=0;
  for (int i=0; i<ndata; i++){sum+=data[i];}
  return(sum);
}

////////////////////////////////////////////////////////////////////////
/// general-purpose statistical operation
////////////////////////////////////////////////////////////////////////

double Statistics::getmax(const double data[], int ndata) const
{
  double x=data[0];
  for(int i=1; i<ndata; i++){
    if(data[i]>x){x=data[i];}
  }
  return x;
}

////////////////////////////////////////////////////////////////////////
/// general-purpose statistical operation
////////////////////////////////////////////////////////////////////////

double Statistics::getmin(const double data[], int ndata) const
{
  double x=data[0];
  for(int i=1; i<ndata; i++){
    if(data[i]<x){x=data[i];}
  }
  return x;
}

////////////////////////////////////////////////////////////////////////
/// general-purpose statistical operation
////////////////////////////////////////////////////////////////////////

int Statistics::getmax(const int data[], int ndata) const
{
  int x=data[0];
  for(int i=1; i<ndata; i++){
    if(data[i]>x){x=data[i];}
  }
  return x;
}

////////////////////////////////////////////////////////////////////////
/// general-purpose statistical operation
////////////////////////////////////////////////////////////////////////

int Statistics::getmin(const int data[], int ndata) const
{
  int x=data[0];
  for(int i=1; i<ndata; i++){
    if(data[i]<x){x=data[i];}
  }
  return x;
}


////////////////////////////////////////////////////////////////////////
/// Arithmetic Means
/// take data[0] ... data[ndata-1]
////////////////////////////////////////////////////////////////////////

double Statistics::arithmeticMeans(const double data[], int ndata) const{
  if(ndata<1) {
      cerr<<"Statistics.arithmeticMeans: error: trying to calculate"
	  <<" average from less than 1 data point!\n";
      exit(-1);
  }

  double sum=0;
  for (int i=0; i<ndata; i++){sum+=data[i];}
  //cout <<"Statistics.arithmeticMeans:  ndata="<<ndata<<" sum="<<sum<<endl;
  return(sum/ndata);
}

double Statistics::arithmeticMeans(const int data[], int ndata) const{
  if(ndata<1) {
      cerr<<"Statistics.arithmeticMeans: error: trying to calculate"
	  <<" average from less than 1 data point!\n";
      exit(-1);
  }

  double sum=0;
  for (int i=0; i<ndata; i++){sum+=data[i];}
  //cout <<"Statistics.arithmeticMeans:  ndata="<<ndata<<" sum="<<sum<<endl;
  return(((double)(sum))/ndata);
}

//#######################################################################
/// Arithmetic Mean
/// take data[imin] ... data[imax]
//#######################################################################

double Statistics::arithmeticMeans(const double data[], 
				   int imin, int imax) const {
  int ndata=imax-imin+1;
  if(ndata<1)
    {
      cerr<<"Statistics.arithm_avg: error: trying to calculate"
	  <<" average from less than 1 data point!\n";
      exit(-1);
    }
  double sum=0;
  for (int i=imin; i<imax+1; i++){sum+=data[i];}
  //cout <<"Statistics.arithmeticMeans:  ndata="<<ndata<<" sum="<<sum<<endl;
  return(sum/ndata);
}


//#######################################################################
// distribution function of a set of double values in x[]
// output in xSorted[], distr[]
//#######################################################################


void Statistics::distrF(const double x[], int nData,
			  double xSorted[], double distr[]) const{

  int* indexSorted = new int[nData];
  quicksort(x, nData, indexSorted);

  for (int i=0; i<nData; i++){
    xSorted[i]=x[indexSorted[i]];
    distr[i]=(i+1.)/nData;
  }
  delete[] indexSorted;
}


//#######################################################################
// quantiles, incl median MT sep2014
// returns the highest value of the ndata*q lowest values
//#######################################################################



double Statistics::quantile(const double x[], int nData, double q) const{
  if( (q<0) || (q>1)){
    cerr<<"error: Statistics::quantile: quantil argument must be in range [0,1]"
	<<endl;
    exit(-1);
  }

  int* indexSorted = new int[nData];

  quicksort(x, nData, indexSorted);
  int iQuantil=(int)( (nData-1)*q);
  double rQuantil=nData*q-iQuantil;
  double xLower=x[indexSorted[iQuantil]];
  double xHigher=(iQuantil<nData-1) ? x[indexSorted[iQuantil+1]] : xLower;
  delete[] indexSorted;
  return (1-rQuantil)*xLower + rQuantil*xHigher;
}

/** quantile for integer incl median: 
 fits with def median by step function
 */

int Statistics::quantile(const int x[], int nData, double q) const{
  if( (q<0) || (q>1)){
    cerr<<"error: Statistics::quantile: quantil argument must be in range [0,1]"
	<<endl;
    exit(-1);
  }
  int* indexSorted = new int[nData];
  quicksort(x, nData, indexSorted);
  int iQuantile=(int)(nData*q);
  int output=x[indexSorted[iQuantile]];
  delete[] indexSorted;
  return output;
}



////////////////////////////////////////////////////////////////////////
/// take data[0] ... data[ndata-1]
/// use definition of descr. statistics (n, not (n-1) in denominator)
////////////////////////////////////////////////////////////////////////

double Statistics::varianceDescr(const double data[], int ndata) const {
  if(ndata<2){
    cerr<<"Statistics.variance: error: trying to calculate variance"
	<<" from less than 2 data points!\n";
    exit(-1);
  }

  double avg=arithmeticMeans(data, ndata);

  double sumxx=0;
  for (int i=0; i<ndata; i++){sumxx+=(data[i]-avg)*(data[i]-avg);}
  double varianceDescr = sumxx/(ndata);
  return(varianceDescr);
}

// inductive variance

double Statistics::variance(const double data[], int ndata) const {
  return ndata/(ndata-1)*varianceDescr(data, ndata);
}

//#######################################################################
/// variance for indices imin...imax (inclusive) descr. statistics
//#######################################################################

double Statistics::varianceDescr(const double data[],
			    int imin, int imax) const{
  int ndata=imax-imin+1;
  if(ndata<2){
    cerr<<"Statistics.variance: error: trying to calculate variance"
	<<" from less than imax-imin+1=2 data points!\n";
    exit(-1);
  }

  double avg=arithmeticMeans(data, imin, imax);

  double sumxx=0;
  for (int i=imin; i<imax+1; i++){sumxx+=(data[i]-avg)*(data[i]-avg);}
  double varianceDescr = sumxx/(ndata);
  //double varianceInd = sumxx/(ndata-1); //!!

  //cout <<"Statistics::variance: ndata="<<ndata
  //     <<" varianceDescr="<<varianceDescr<<endl;
  return(varianceDescr);
}

// inductive variance

double Statistics::variance(const double data[], 
			    int imin, int imax) const {
  return (imax-imin+1)/(imax-imin)*varianceDescr(data, imin, imax);
}


////////////////////////////////////////////////////////////////////////
/// Descriptive Covariance
/// take data[0] ... data[ndata-1]
////////////////////////////////////////////////////////////////////////

double Statistics::covariance(const double xdata[], const double ydata[],int ndata)
{
  if(ndata<2)
    {
      cerr<<"Statistics.variance: error: trying to calculate variance"
	  <<" from less than 2 data points!\n";
      exit(-1);
    }

  double avgx = arithmeticMeans(xdata, ndata);
  double avgy = arithmeticMeans(ydata, ndata);
  double sumxy=0;

  for (int i=0; i<ndata; i++)
    {
      sumxy += (xdata[i]-avgx)*(ydata[i]-avgy);
    }
  //double covarianceInd = sumxy/(ndata-1);
  double covarianceDescr = sumxy/(ndata);
  //cout <<"Statistics::sigmaxy: ndata="<<ndata 
  //<<" covarianceDescr="<<covarianceDescr<<endl;
  return(covarianceDescr);
}

////////////////////////////////////////////////////////////////////////
///  Descriptive correlation take data[0] ... data[ndata-1]
////////////////////////////////////////////////////////////////////////

double Statistics::correlation(const double xdata[], const double ydata[], int ndata){

  if(ndata<3){
    cerr<<"Statistics.correlation: error: trying to calculate correlation"
	<<" from less than 3 data points!\n";
    exit(-1);
  }

  double sxx=covariance(xdata, xdata, ndata);
  double sxy=covariance(xdata, ydata, ndata);
  double syy=covariance(ydata, ydata, ndata);
  double r=sxy/sqrt(sxx*syy);
  return r;
}

////////////////////////////////////////////////////////////////////////
/// calculates autocorrelation function from xdata 
/// for values 0 ... nshiftmax-1 of the shift index 
/// output: auto-correlation function "corrfun" with 
/// corrfun[0]=1, ..., corrfun[nshiftmax-1] in [-1,1] 
/// NOTICE: not coded efficiently, scales as ndata*nshiftmax (quick hack)
////////////////////////////////////////////////////////////////////////

void Statistics::calcCorrelationFun(const double xdata[], int ndata,
				    int nshiftmax, double corrfun[])
{

  if(ndata-nshiftmax<3)
    {
      cerr<<"Statistics.calcCorrelationFun: error: trying to calculate"
	  <<" autocorrelation function with ndata-nshiftmax<3"
	  <<" data points!\n";
      exit(-1);
    }

  //arne (9-9-05): avoid static memory assignment!)
  //  if(ndata>=NMAX){
  //  cerr<<"Statistics.calcCorrelationFun: error: ndata = "
  //<<ndata<<" exceeds constants NMAX = "<<NMAX<<"! Exit\n";
  // exit(-1);
  //}
  //double  data[NMAX];

  double* ydata = new double[ndata];  //(nmax would be sufficient) 

  for (int nshift=0; nshift<nshiftmax; nshift++)
    {
      int nmax=ndata-nshift;
      for (int i=0; i<nmax; i++)
	{
	  ydata[i]=xdata[i+nshift];
	}
      corrfun[nshift]=correlation(xdata,ydata,nmax);
    }
  delete[] ydata; 
}

////////////////////////////////////////////////////////////////////////
/// Martin 17-2-11 calculates cross correlation function from xdata, ydata 
/// for values -nshiftmax ... nshiftmax of the shift index 
/// output: cross-correlation function "corrfun" with 
/// corrfun[nshiftmax]=CCF with shift=0
/// if i>nshiftmax, then ydata is advanced with respect to xdata
/// NOTICE: not coded efficiently, scales as ndata*nshiftmax (quick hack)
////////////////////////////////////////////////////////////////////////

void Statistics::calcCrossCorrelationFun(const double xdata[], const double ydata[],
					 int ndata,
				    int nshiftmax, double corrfun[])
{

  if(ndata-2*nshiftmax-1<3)
    {
      cerr<<"Statistics.calcCrossCorrelationFun: error: trying to calculate"
	  <<" cross correlation function with ndata-2*nshiftmax-1<3"
	  <<" data points!\n";
      exit(-1);
    }

  //arne (9-9-05): avoid static memory assignment!)
  //  if(ndata>=NMAX){
  //  cerr<<"Statistics.calcCorrelationFun: error: ndata = "
  //<<ndata<<" exceeds constants NMAX = "<<NMAX<<"! Exit\n";
  // exit(-1);
  //}
  //double  data[NMAX];

  double* udata = new double[ndata];  //(nmax would be sufficient) 
  double* vdata = new double[ndata];  //(nmax would be sufficient) 

  for (int di=-nshiftmax; di<nshiftmax+1; di++){
    int ixmin=(di<0) ? 0 : di;
    int iymin=(di<0) ? -di : 0;
    int nmax=ndata-abs(di);
    for (int i=0; i<nmax; i++){
      udata[i]=xdata[i+ixmin];
      vdata[i]=ydata[i+iymin];
    }
    corrfun[nshiftmax+di]=correlation(udata,vdata,nmax);
  }
  delete[] udata; 
  delete[] vdata; 
}


////////////////////////////////////////////////////////////////////////
/// arne 19.1.2005
/// calculates autocorrelation function from xdata 
/// for values 0 ... nshiftmax-1 of the shift index 
/// but detrends the xdata first with using the symmetrical (ordinary)
/// arithmetic moving average of periodMA.
/// output: auto-correlation function "corrfun" with 
/// corrfun[0]=1, ..., corrfun[nshiftmax-1] in [-1,1] 
/// if length of xdata isn't sufficient the success flag will
/// be false and corrfun will be left unchanged!
/// NOTICE: not coded efficiently, scales as ndata*nshiftmax (quick hack)
////////////////////////////////////////////////////////////////////////

void Statistics::calcCorrelationFunDetrend(const double xdata[], int ndata,
					   int nshiftmax, int periodMA,
					   double corrfun[], bool& success)
{
  success = true;

  //(i) calc moving average:

  if(ndata<=periodMA){
    cerr<<"Statistics.calcCorrelationFunDetrend: error: "
	<<" cannot calc moving average with ndata<=periodMA";
    cerr<<"   ndata = "<<ndata<<" and periodMA = "<<periodMA<<endl;
    cerr<<"return with success=false flag and unchanged parameters\n\n";
    success=false;
    return;
  }

  //arne (9-9-05): avoid static memory assignment!) subsitute NMAX --> ndata
  //double movAvg[NMAX]; 
  double* movAvg = new double[ndata];
  for(int i=0;i<ndata;i++){movAvg[i]=0;}
  //calcMA(xdata, ndata, periodMA, movAvg);
  //arne & martin 21-1-2005:
  //gauss is better (no edge)
  //FT of gauss not oscillatory  
  smoothGauss(xdata, ndata, periodMA, movAvg);
  //(ii) detrend data
  int iBoundary = periodMA/2; //integer division!
  int nDetrendedData = ndata-2*iBoundary; //array size decreased by MA!!!

  if(nDetrendedData-nshiftmax<3){
    cerr<<"Statistics.calcCorrelationFunDetrend: error: "
	<<" cannot calc correlation function from  from (nDetrendedData-nshiftmax) = "
	<<nDetrendedData<<"-"<<nshiftmax<<" datapoints!"<<endl;
    cerr<<"return with success=false flag and unchanged parameters";
    success=false;
    delete[] movAvg;
    return;
  }

  for (int i=0;i<nDetrendedData; i++){
    //for efficiency use already allocated array to
    //save detrended xdata!!!
    movAvg[i] = xdata[i+iBoundary] - movAvg[i+iBoundary];
  }

  calcCorrelationFun(movAvg, nDetrendedData, nshiftmax, corrfun); 

  delete[] movAvg;
}

////////////////////////////////////////////////////////////////////////
/// calculates linear regression y=a+bx and "Bestimmtheitsmass" B
/// take data[0] ... data[ndata-1]
////////////////////////////////////////////////////////////////////////

void Statistics::linRegression(const double xdata[], const double ydata[], 
			       int ndata,
			       double& a, double& b, 
			       double& B, double& sumdev)
{
  linRegression(xdata,ydata,ndata,0,ndata-1,a,b,B,sumdev);
}


////////////////////////////////////////////////////////////////////////
/// calculates linear regression y=a+bx and "Bestimmtheitsmass" B
/// from a time series (xdata=0,1,2,3,...)
////////////////////////////////////////////////////////////////////////

void Statistics::linRegression(const double ydata[], 
			       int ndata,
			       double& a, double& b, 
			       double& B, double& sumdev)
{
  double* xdata = new double[ndata];
  for (int i=0; i<ndata; i++){
    xdata[i]=i;
  }
  linRegression(xdata,ydata,ndata,0,ndata-1,a,b,B,sumdev);
  delete[] xdata;
}


////////////////////////////////////////////////////////////////////////
/// calculates linear regression y=a+bx and "Bestimmtheitsmass" B
/// of a section of (xdata, ydata) in the index range [imin, imax]
/// take data[imin] ... data[imax] -> (imax-imin+1) data points
/// sumdev=sum(square deviations)
////////////////////////////////////////////////////////////////////////

void Statistics::linRegression(const double xdata[], const double ydata[], 
			       int ndata, int imin, int imax, 
			       double& a, double& b, 
			       double& B, double& sumdev)
{
  int nSec=imax+1-imin;
  if(nSec<3)
    {
      cerr<<"Statistics.linRegression: error: trying to calculate\n "
	  <<" linear regression"
	  <<" from less than 3 data points!\n";
      exit(-1);
    }

  if(imin<0)
    {
      cerr<<"Statistics.linRegression: Error: imin<0"<<endl; 
      exit(-1);
    }

  //if((imax>NMAX)||(ndata>NMAX)||(imax>ndata-1)){
  //   cerr<<"linRegression: Error: imax or ndata >NMAX="<<NMAX
  // <<" or imax>ndata-1"<<endl; exit(-1);
  //  }

  //arne (9-9-05): avoid static memory assignment!)
  //subsitute NMAX --> nSec
  //double xdataSec[NMAX+1];
  //double ydataSec[NMAX+1];


  double* xdataSec = new double[nSec];
  double* ydataSec = new double[nSec];

  for (int i=imin; i<=imax; i++)
    {
      xdataSec[i-imin] = xdata[i];
      ydataSec[i-imin] = ydata[i];
    }

  double xbar=arithmeticMeans(xdataSec, nSec);
  double ybar=arithmeticMeans(ydataSec, nSec);

  cout <<"Statistics.linRegression: nSec="<<nSec<<" xbar="<<xbar<<" ybar="<<ybar<<endl;
  
  // do not use variance,covariance to be independent of 
  // use of induct/descript. stat in variance
  // use sum(x_i^2) etc, not sum(x_i-xbar)^2 to check hand-written calc.

  double sumxx=0;
  double sumxy=0;
  double sumyy=0;
  for (int i=0; i<nSec; i++){
    sumxx += (xdataSec[i])*(xdataSec[i]);
    sumxy += (xdataSec[i])*(ydataSec[i]);
    sumyy += (ydataSec[i])*(ydataSec[i]);
  }

  cout<<"Statistics.linRegression: sumxx="<<sumxx
      <<" sumxy="<<sumxy<<" sumyy="<<sumyy<<endl;

  double sxx=sumxx/nSec - xbar*xbar;
  double sxy=sumxy/nSec - xbar*ybar;
  double syy=sumyy/nSec - ybar*ybar;

  double r=sxy/sqrt(sxx*syy);

  b = sxy/sxx;
  a = ybar - b*xbar;
  B = r*r;
  sumdev=0;
  //cout <<"Lin Regression: imin="<<imin<<" imax="<<imax<<endl;

  for (int i=0; i<nSec; i++)
    {
      sumdev += SQR(a+b*xdataSec[i]-ydataSec[i]);
      //cout <<"i+imin="<<(i+imin)<<" sumdev="<<sumdev<<endl;
    }

  if(false){
    cout<<" linRegression: "
        <<" xdata["<<imin<<"]="<<xdata[imin]
        <<" xdata["<<imax<<"]="<<xdata[imax]
        <<" ydata["<<imax<<"]="<<ydata[imax]
	<< endl
        <<" xbar="<<xbar
        <<" ybar="<<ybar
        <<" sxx="<<sxx
        <<" syy="<<syy
        <<" sxy="<<sxy
        <<" r="<<r
	<<endl;
  }
  delete[] xdataSec;
  delete[] ydataSec;
}


//###########################################################
// MT 2022-09: Smoothing with triang kernel with vectors
// dnSmooth: half-width
//###########################################################

vector<double> Statistics::smoothTriang(const vector<int> data, 
					int dnSmooth)
{
  vector <double> dataDouble(data.begin(), data.end());
  return smoothTriang(dataDouble,dnSmooth);
}


vector<double> Statistics::smoothTriang(const vector<double> data, 
					int dnSmooth)
{

  int ndata=data.size();
  if(dnSmooth<0){
    cerr<<"Statistics.smoothTriang: error: kernel has negative half-width"
	<<dnSmooth<<"!\n";
    
    cerr<<" doing nothing"<<endl; return data;
    //exit(-1);
  }
    
  if(ndata<2*dnSmooth+1){
    if(false){
      cerr<<"Statistics.smoothTriang: warning: ndata="<<ndata
	  <<" 2*dnSmooth+1="<<2*dnSmooth+1
	  <<"\n trying to calculate"
	  <<" smoothed values from less than 2*dnSmooth+1 data points!";
      cerr<<" doing nothing"<<endl;
    }
    
    return data;
    //exit(-1);
  }


  // create the return vector (by value)
  
  vector<double>smoothedData(ndata,0);  // initialized with zeroes

  
  // create non-normalied kernel (dnSmooth=1: [1,2,1] )
  
  vector<double> kernel(2*dnSmooth+1);
  for(int di=-dnSmooth; di<=dnSmooth; di++){
    kernel[di+dnSmooth]=1+dnSmooth-fabs(di);
  }

  // calculate regular norm for points inside and normalize kernel

  double norm=0;
  for(int j=0; j<int(kernel.size()); j++){
    norm+=kernel[j];
  }
  for(int j=0; j<int(kernel.size()); j++){
    kernel[j]/=norm;
  }

  if(false){
    cout<<"Statistics.smoothTriang: dnSmooth= "<<dnSmooth<<endl;
    for(int j=0; j<=2*dnSmooth; j++){
      cout<<"j="<<j<<" kernel[j]="<<kernel[j]<<endl;
    }
  }

  // calculate left and right boundary values

  smoothedData[0]=data[0];
  smoothedData[ndata-1]=data[ndata-1];

  for(int i=1; i<dnSmooth; i++){
    double boundaryNorm=0;
    for(int j=0; j<=2*i; j++){
      smoothedData[i]+=kernel[dnSmooth-i+j]*data[j];
      smoothedData[ndata-1-i]+=kernel[dnSmooth-i+j]*data[ndata-1-2*i+j];
      boundaryNorm+=kernel[dnSmooth-i+j];
    }
    
    smoothedData[i]/=boundaryNorm;
    smoothedData[ndata-1-i]/=boundaryNorm;

    //cout <<"Statistics.smoothTriang: smooth boundary data"
    //	 <<" i="<<i<<" ndata-1-i="<<ndata-1-i
    //	 <<" boundaryNorm="<<boundaryNorm<<endl;
  }

  
  // calculate regular middle values where full kernel is availaible
  // max... because otherwise double count for dnSmooth=0
  
  for(int i=max(1,dnSmooth); i<min(ndata-1, ndata-dnSmooth); i++){
    for(int j=-dnSmooth; j<=dnSmooth; j++){
      smoothedData[i]+=kernel[dnSmooth+j]*data[i+j];
    }
  }

  return smoothedData;
}
  





////////////////////////////////////////////////////////////////////////
/// static NMAX is not very elegant
////////////////////////////////////////////////////////////////////////

void Statistics::smoothGauss(const double data[], int ndata, 
			     double dnSmooth, double smoothedData[])
{
  
  static const int NMAX=100000;
  int nCenter=NMAX/2;
  int maxWidth=min( static_cast<int>(2.*dnSmooth), static_cast<int>(NMAX/2));
  
  if(ndata<3)
    {
      cerr<<"Statistics.smoothGauss: error: trying to calculate smoothed values"
	  <<" from less than 3 data points!\n";
      exit(-1);
    }
  
  if(ndata>=NMAX)
    {
      cerr<<"Statistics.smoothGauss: error: data field larger than"<<NMAX<<endl;
      exit(-1);
    }
  
  // create table of values of (not normalized) Gaussian values
  double* gaussfactor = new double[NMAX+1];
  
  for (int j=nCenter-maxWidth; j<=nCenter+maxWidth; j++){
    double dn=j-nCenter;
    gaussfactor[j] = exp(-SQR(dn/dnSmooth));
    // cout <<"gaussfactor["<<j<<"]="<<gaussfactor[j]<<endl;
  }
  
  // make normalization

  double* norm = new double[NMAX+1];
  for (int i=0; i<ndata; i++){
    norm[i]=0;
    int actWidth=min(maxWidth, min(i, ndata-1-i));
    for (int j=i-actWidth; j<=i+actWidth; j++){
      norm[i] += gaussfactor[nCenter + i-j];
    }
    if(norm[i]==0){
      cerr<<"Statistics.smoothGauss:error in normalization: norm[i]=0\n!"; exit(-1);
    }
  }

  // make weighted sum

  double* wsum = new double[NMAX+1];
  for (int i=0; i<ndata; i++){
    wsum[i]=0;
    int actWidth=min(maxWidth, min(i, ndata-1-i));
    for (int j=i-actWidth; j<=i+actWidth; j++){
      wsum[i] += gaussfactor[nCenter + i-j]*data[j];
    }
  }

  // make result

  for (int i=0; i<ndata; i++){
    if(i==(ndata-1)){
      cout<<"i = "<<i<<" wsum = "<<wsum[i]<<" norm = "<<norm[i]<<endl;
    }
    smoothedData[i] = wsum[i]/norm[i];
  }

  delete[] gaussfactor;
  delete[] norm;
  delete[] wsum;
} 

////////////////////////////////////////////////////////////////////////
/// arne 28.9.2004:
/// general gauss smoothing: calc gaussian weighted average for each data point from ydata.
/// xdata is needed to define a distance metrics.
/// calc convolution for each datapoint 
////////////////////////////////////////////////////////////////////////

void Statistics::gSmoothGauss(const double xdata[],
			      const double ydata[], 
			      const int ndata, 
			      const double sigma,
			      double smoothedData[])
{
  cout<<"Statistics::gSmoothGauss with ndata = "<<ndata<<endl;

  if(ndata<3)
    {
      cerr<<"Statistics.gSmoothGauss: error: trying to calculate smoothed values"
	  <<" from less than 3 data points! Exit(-1)\n";
      exit(-1);
    } 
  
  for(int i=0; i<ndata; i++)
    {
      double xCenter = xdata[i]; 
      double norm=0;
      double weightedSum=0;
      double gaussfactor=0;
      //do not cut off here! --> not very efficient (see calcSymmetricEMA)
      for(int j=0; j<ndata;j++) 
	{
	  double dx = xdata[j]-xCenter;
	  gaussfactor    = exp(-SQR(dx/(2*sigma))); //arne: faktor 2 erst am 14.8.2006!!!
	  norm          += gaussfactor;
	  weightedSum   += gaussfactor*ydata[j];
	  if(false)
	    {
	      cout<<"i="<<i<<",j="<<j
		  <<" xdata[j]="<<xdata[j]
		  <<" ydata[j]="<<ydata[j]
		  <<" dx = "<<dx
		  <<", gaussfactor="<<gaussfactor<<endl;
	    }
	}
    if(norm==0)
      {
	cerr<<"Statistics.gSmoothGauss:error: norm[i]=0\n!"; exit(-1);
      }
    smoothedData[i] = weightedSum/norm;
    
    if(false)
      {
	cout<<"weightedSum = "<<weightedSum<<", norm = "<<norm
	    <<", gaussfactor(ndata-1)="<<gaussfactor<<endl;
      }
    } //of for(i)
} 

////////////////////////////////////////////////////////////////////////
/// arne 8-9-05
/// <ul>
/// <li>return average with gaussian weight of width sigma around xCenter for ydata
/// xdata is needed to define a distance metrics.
/// <li>Cut-off is 4*sigma to enhance performance.
/// </ul>
////////////////////////////////////////////////////////////////////////

double Statistics::gGaussAverage(const double xdata[],
				 const double ydata[],
				 const int ndata, 
				 const double xCenter,
				 const double sigma)
{
  double cutOff = 4*sigma;
  
  if(false) cout<<"Statistics::gGaussAverage with ndata = "<<ndata<<endl;
  
  if(ndata<3)
    {
      cerr<<"Statistics.gGaussAverage: error: trying to calculate smoothed values"
	  <<" from less than 3 data points! Exit(-1)\n";
      exit(-1);
    } 
  
  double norm=0;
  double weightedSum=0;
  double gaussfactor=0;
  double gaussAverage=0;

  for(int i=0; i<ndata;i++)
    {
      double dx = xdata[i]-xCenter;
      if(fabs(dx)<=cutOff)
	{
	  gaussfactor    = exp(-SQR(dx/(2*sigma)));//arne: faktor 2 erst am 14.8.2006!!!
	
	  norm          += gaussfactor;
	  weightedSum   += gaussfactor*ydata[i];
	}
    }

  if(norm==0){
    cerr<<"Statistics.gSmoothGauss:error: norm[i]=0\n!";
    exit(-1);
  }

  gaussAverage = weightedSum/norm;

  if(false)
    {
      cout<<"gGaussAverage::weightedSum = "<<weightedSum<<", norm = "<<norm
	  <<", gaussfactor(ndata-1)="<<gaussfactor<<endl;
      cout<<"return gaussAverage = "<<gaussAverage<<endl;
    }
  
  return(gaussAverage);
} 

////////////////////////////////////////////////////////////////////////
/// arne (september 05)
/// <ul>
/// <li> calculates linear regression y=a+bx and "Bestimmtheitsmass" B 
/// of a section of (xdata, ydata) in the index range [imin, imax].
/// <li> take data[imin] ... data[imax] -> (imax-imin+1) data points
/// <li> sumdev=sum(square deviations)
/// <li> Used in "bollinger" tool.
/// </ul>
////////////////////////////////////////////////////////////////////////

void Statistics::gaussLinRegression(const double xdata[], const double ydata[], 
				    int ndata, 
				    double xCenter, double sigma,
				    double& a, double& b, 
				    double& B, double& sumdev)
{ 
  double cutOff = 4*sigma;

  if(ndata<3){
    cerr<<"Statistics.gaussLinRegression: error: trying to calculate\n "
	<<" linear regression"
	<<" from less than 3 data points!\n";
    exit(-1);
  }
  
  double xbar=gGaussAverage(xdata, xdata, ndata, xCenter, sigma);
  double ybar=gGaussAverage(xdata, ydata, ndata, xCenter, sigma);

  if(false)cout <<"Statistics.gaussLinRegression: xCenter="<<xCenter<<", xbar="<<xbar<<" ybar="<<ybar<<endl;

  //here gauss-weighted sum 
  double sumxx=0;
  double sumxy=0;
  double sumyy=0;
  double norm=0;
  for (int i=0; i<ndata; i++)
    {
      double dx = xdata[i]-xCenter; //distance relative to xCenter
      if(fabs(dx)<=cutOff)
	{
	  double gaussfactor = exp(-SQR(dx/(2*sigma)));//arne: faktor 2 erst am 14.8.2006!!!
	  norm  += gaussfactor;
	  sumxx += gaussfactor*(xdata[i])*(xdata[i]);
	  sumxy += gaussfactor*(xdata[i])*(ydata[i]);
	  sumyy += gaussfactor*(ydata[i])*(ydata[i]);
	}
    }

  if(norm==0){
    cerr<<"Statistics.gSmoothGauss:error: norm[i]=0\n!"; exit(-1);
  }

  if(false)
    {
      cout <<"Statistics.linRegression: sumxx="<<sumxx
	   <<" sumxy="<<sumxy<<" sumyy="<<sumyy<<endl;
    }

  double sxx=sumxx/norm - xbar*xbar;
  double sxy=sumxy/norm - xbar*ybar;
  double syy=sumyy/norm - ybar*ybar;

  double r=sxy/sqrt(sxx*syy);

  b = sxy/sxx;
  a = ybar - b*xbar;
  B = r*r;
  sumdev=0;
  //cout <<"Lin Regression: imin="<<imin<<" imax="<<imax<<endl;
  for (int i=0; i<ndata; i++){
    sumdev += SQR(a+b*xdata[i]-ydata[i]);
    //cout <<"i+imin="<<(i+imin)<<" sumdev="<<sumdev<<endl;
  }

  if(false){
    cout<<" gaussLinRegression: "
        <<" xdata["<<0<<"]="<<xdata[0]
        <<" xdata["<<(ndata-1)<<"]="<<xdata[ndata-1]
        <<" ydata["<<(ndata-1)<<"]="<<ydata[ndata-1]
	<< endl
        <<" xbar="<<xbar
        <<" ybar="<<ybar
        <<" sxx="<<sxx
        <<" syy="<<syy
        <<" sxy="<<sxy
        <<" r="<<r
	<<endl;
  }
}



////////////////////////////////////////////////////////////////////////
/// asymmetrical exponential moving average
/// Arne, 20.10.2003
/// <ul>
/// <li> Input: time-series data with constant dt and averaging period
/// in units of dt (=tau),
/// <li> output: moving average. Notice: data size decreased by one period!
/// boundary regions filled by zeros instead producing shorter
/// time series to avoid hidden shifts!
/// </ul>
////////////////////////////////////////////////////////////////////////

void Statistics::calcEMA(const double data[], int ndata,
			 double tau, double movEAvg[])
{
  calcEMA(data,ndata,data[0],tau,movEAvg);
}

////////////////////////////////////////////////////////////////////////
/// asymmetrical exponential moving average
/// Arne, 20.10.2003
/// <ul>
/// <li> Input: time-series data with constant dt and averaging period
/// in units of dt (=tau),
/// <li> user defines an init. value:  movEAvg[0]=initVal;
/// <li> output: moving average. Notice: data size decreased by one period!
/// boundary regions filled by zeros instead producing shorter
/// time series to avoid hidden shifts!
/// </ul>
////////////////////////////////////////////////////////////////////////

void Statistics::calcEMA(const double data[], int ndata, double initVal,
			 double tau, double movEAvg[]){
  if(ndata-tau<1){
    cerr<<"Statistics.calcEMA: error: data size is less than"
	<<" period+1; provide more data!"<<endl;
    exit(-1);
  }

  double alpha = 1-exp(-1/tau);

  ///assumption for first entry: \hat{y}_0 = y_0 

  movEAvg[0]=initVal;

  for (int i=1; i<ndata; i++){
    movEAvg[i] = alpha*data[i] + (1-alpha)*movEAvg[i-1];
  }

  // test output

  if(false){
    cout <<"Testing Statistics.calcEMA: ndata="<<ndata<<endl;
    for (int i=0; i<ndata; i++){
      cout<<"i="<<i<<" data[i]="<<data[i]<<" movEAvg[i]="<<movEAvg[i]<<endl;
    }
  }
}



///#####################################################################
/// symmetrical exponential moving average 
/// version for data with constant sampling intervals dt
/// <ul>
/// <li> Input: time-series data with constant dt and averaging period
/// in units of dt (tau=dnSmooth*dt),
/// <li> output: symmetrical moving exponential average.
/// <li> numerical cut off is 7*dnSmooth
/// <li> Boundary is automatically corrected by calculating norm 
///      for each data point separatly
/// </ul>
///#####################################################################


void Statistics::calcSymmetricEMA(const double data[], int ndata,
				  double dnSmooth, double smoothedData[]){

  //cout<<"Statistics::calcSymmetricEMA with ndata = "<<ndata<<endl;

  if(ndata<3){
    cerr<<"Statistics.calcSymmetricEMA: error: trying to calculate smoothed values"
	<<" from less than 3 data points! Exit(-1)\n";
    exit(-1);
  }

  //numerical cut-off: dnSmooth is 1/e and  e^{-7} = 0.000091
  int cutOff = static_cast<int>(7*dnSmooth);
  int maxWidth=min( cutOff, ndata); 

  for(int i=0; i<ndata; i++){
    double norm=0;
    double weightedSum=0;
    double expfactor=0;

    for(int j=max(0,i-maxWidth); j<=min(ndata-1, i+maxWidth); j++){
      double dn = i-j; //equidistant data
      expfactor = exp( -fabs(dn/dnSmooth) );
      norm          += expfactor;
      weightedSum   += expfactor*data[j];
      if(false){
	      cout<<"i="<<i<<",j="<<j
		  <<" data[j]="<<data[j]
		  <<" dn = "<<dn
		  <<", expfactor="<<expfactor<<endl;
      }
    }

    if(norm==0){
 	cerr<<"Statistics.calcSymmetricEMA:error: norm[i]=0\n!"; exit(-1);
    }
    smoothedData[i] = weightedSum/norm;
    
    if(false){
	cout<<"i="<<i
	    <<", weightedSum = "<<weightedSum<<", norm = "<<norm
	    <<", expfactor(ndata-1)="<<expfactor
	    <<", smoothedData="<<smoothedData[i]<<endl;
    }
  } //of for(i)

}


///#####################################################################
/// symmetrical exponential moving average 
/// version for data with nonequal sampling intervals dt or if resampling wanted
/// <ul>
/// <li> Input: (i)  times tdata with increasing values, size ndata
///            (ii)  time-series data for each tdata point
///            (iii) smoothing scale dtSmooth and output sampling interval dtout
/// <li> output: resampled tdataOut, and EMA dataOut
/// <li> if only resampling or filling of missing lines wanted, 
///      set dtOut=normal dtIn and dtSmooth=0
/// <li> It is in the responsibility of the caller to check array bounds of dataOut
/// <li> Boundary is automatically corrected by calculating norm 
///      for each data point separatly
/// </ul>
///#####################################################################


void Statistics::calcSymmetricEMA(const double tdata[],const double data[],int ndata,
				  double dtSmooth, double dtOut,
				  double tdataOut[], double dataOut[]){

  //cout<<"Statistics::calcSymmetricEMA with resampling, ndata = "<<ndata<<endl;

  if(ndata<3){
    cerr<<"Statistics.calcSymmetricEMA: error: trying to calculate smoothed values"
	<<" from less than 3 data points! Exit(-1)\n";
    exit(-1);
  }

  // check input

  for(int j=1; j<ndata; j++){
    if(tdata[j]<tdata[j-1]){
      cerr<<"Statistics.calcSymmetricEMA: Error: tdata are not in increasing order"
	  <<" j="<<j<<" tdata[j]="<<tdata[j]<<" tdata[j-1]="<<tdata[j-1]<<endl;
      exit(-1);
    }
  }

  // prepare loop

  double tmax=tdata[ndata-1];
  double tmin=tdata[0];
  int ndataOut=(int)((tmax-tmin)/dtOut)+1;
  double dtCutOff = 5*dtSmooth;  //numerical cut-off at  e^{-5} = 0.007

  // if only resampling wanted

  double dtAvgIn=(tmax-tmin)/(ndata-1);
  if(dtSmooth<0.01*dtAvgIn){
    dtSmooth=0.01*dtAvgIn;
    dtCutOff=dtAvgIn;
  }

  for(int i=0; i<ndataOut; i++){           // i relates to output
    double t=tmin+i/(ndataOut-1.)*(tmax-tmin);
    double norm=0;
    double weightedSum=0;
    double expfactor=0;
    bool goOn=true;

    for(int j=0; (j<ndata)&&goOn; j++){ // j relates to input

      if(tdata[j]>=t-dtCutOff){
        expfactor = exp( -fabs((tdata[j]-t)/dtSmooth) );
        norm          += expfactor;
        weightedSum   += expfactor*data[j];
	goOn           = (tdata[j]<=t+dtCutOff);
        if(false){
	      cout<<"i="<<i<<",j="<<j
		  <<" data[j]="<<data[j]
		  <<" t="<<t<<"  tdata[j]="<<tdata[j]
		  <<" expfactor="<<expfactor<<endl;
	}
      }
    }

    if(norm==0){
      cerr<<"Statistics.calcSymmetricEMA:error: i="<<i<<" t="<<t
	  <<" norm[i]=0\n!"; exit(-1);
    }

    tdataOut[i]=t;
    dataOut[i] = weightedSum/norm;

    
    if(false){
	cout<<"i="<<i
	    <<", weightedSum = "<<weightedSum<<", norm = "<<norm
	    <<" tdataOut[i]="<<tdataOut[i]
	    <<" dataOut[i]="<<dataOut[i]
	    <<endl;
    }

  } //outer loop output sampling i

}




////////////////////////////////////////////////////////////////////////
/// symmetrical moving average
/// <ul>
/// <li> Input: time-series data with constant dt and averaging period
/// in units of dt,
/// <li> output: moving average. Notice: data size decreased by one period!
/// boundary regions filled by zeros instead producing shorter
/// time series to avoid hidden shifts!
/// </ul>
////////////////////////////////////////////////////////////////////////

void Statistics::calcMA(const double data[], int ndata,
			int period, double movAvg[])
{
  if(ndata-period<1){
    cerr<<"Statistics.calcMA: error: data size is less than"
	<<" period+1; provide more data!"<<endl;
    exit(-1);
  }

  bool even=(period%2==0);
  int m=(even) ? period/2 : (period-1)/2;    // (period/2 tut's auch)

  // set boundaries to zero

  for(int i=0;i<m; i++){ 
    movAvg[i]=0;
    movAvg[ndata-1-i]=0;
  }

  // determine first nontrivial value movAvg[m]

  movAvg[m]=0;
  for (int j=0; j<2*m+1; j++){
    movAvg[m] += data[j];
  }
  if(even){
    movAvg[m] -= 0.5*(data[0]+data[2*m]);
  }

  // determine the rest by adding most recent + subtracting 1 period old data

  for (int i=m+1; i<ndata-m; i++){
    double change=(even)
      ? 0.5*(data[i+m]-data[i-m]+data[i+m-1]-data[i-m-1])
      : data[i+m]-data[i-m-1];
    movAvg[i] =movAvg[i-1]+change;
  }

  // divide everything by period

  for (int i=0; i<ndata; i++){
    movAvg[i] /= period;
  }

  // test output

  if(false){
    cout <<"Testing Statistics.calcMA: ndata="<<ndata<<endl;
    for (int i=0; i<ndata; i++){
      cout<<"i="<<i<<" data[i]="<<data[i]<<" movAvg[i]="<<movAvg[i]<<endl;
    }
  }

}


////////////////////////////////////////////////////////////////////////
/// "Saisonbereinigung", stationary case
/// <ul>
/// <li> Input: 
/// <ul> <li> stationary (trend free) time-series 
///      with constant dt and known integer period (units of dt), 
///      <li> the number of available periods
///      <li> The first index taken from the data
/// </ul>
/// <li> output: "Saisonfigur" of length=period.
/// </ul>
////////////////////////////////////////////////////////////////////////

void Statistics::calculate_detr_saison(const double data[], int ndata,
				       int period, int nPeriods, int istart,
				       double data_saison[])
{
  int n=static_cast<int>((ndata-istart)/period);  // number of usable periods
  if(n<nPeriods){
    cerr<<"Statistics.calculate_saison: error: data size is less than"
	<<nPeriods<<" periods+initial offset; provide more data!"<<endl;
    exit(-1);
  }

  double avg=arithmeticMeans(data,istart,istart+nPeriods*period-1);
  for (int k=0; k<period; k++){
    data_saison[k]=0;
    for (int l=0; l<nPeriods; l++){
      data_saison[k] +=data[istart+l*period+k]-avg;
    }
    data_saison[k] /= nPeriods;
  }
  if(true){
    cout <<"Statistics.calculate_detr_saison: ndata="<<ndata
	 <<" istart="<<istart<<" number of periods="<<n
	 <<" avg="<<avg<<endl;
  }

}


////////////////////////////////////////////////////////////////////////
/// "Saisonbereinigung", general case
/// <ul>
/// <li> Input: time-series data with constant dt with known integer period
/// in units of dt,
/// <li> output: "Saisonfigur".
/// </ul>
////////////////////////////////////////////////////////////////////////

void Statistics::calculate_saison(const double data[], int ndata,
				  int period, double data_saison[])
{
  int n = static_cast<int>(ndata/period)-1;  // number of usable periods (need 1 for avg)
  bool even = (period%2==0);
  int istart = (even) ? period/2 : (period-1)/2;    // (remark period/2 would work as well)

  if(n<1){
    cerr<<"Statistics.calculate_saison: error: data size is less than"
	<<" 2*period; provide more data!"<<endl;
    exit(-1);
  }
 
  //detrend with symmetric ordinary moving average
  //arne (9-9-05: avoid static memory assignment!) substitute NMAX --> ndata
  //  double maData[NMAX+1];    // moving average
  //  double detrData[NMAX+1];  // detrended data
 
  double* maData   = new double[ndata+1];  // moving average
  double* detrData = new double[ndata+1];  // detrended data

  calcMA(data, ndata, period, maData);

  for (int i=0; i<ndata; i++)
    {
      detrData[i]=data[i]-maData[i];
    }

  // apply special case for stationary data

  //double data_prelim_saison[NMAX+1];    // preliminary saisonal component
  double* data_prelim_saison = new double[ndata+1];    // preliminary saisonal component

  calculate_detr_saison(detrData, ndata, period, n, istart, data_prelim_saison);

  // undo the shift neeeded by moving average
  for (int k=0; k<period; k++)
    {
      int ktrue=(k+period-istart)%period;
      data_saison[k]=data_prelim_saison[ktrue];
    }

  // apply correction to make saison average of data_saison=0

  double avg=arithmeticMeans(data_prelim_saison, period);
  for (int k=0; k<period; k++)
    {
      data_saison[k] -= avg;
    }
  delete[] maData;
  delete[] detrData;
  delete[] data_prelim_saison;
}



//###################################################
/*
  quicksort routine
  input: array of doubles and number of data elements (index running
         from 0 to (n-1) to be sorted
  output: integer array giving the indices for a sort to ascending
          values
  Algorithm:
  sort first pairs of neighboring values, then create sorted 4-element chunks
  by sorting neighboring pairs of these pairs using the information that the
  original pairs are already sorted. Then created sorted 8-element chunks by
  sorting pairs of the 4-element chunks and so on 
  details: 
  If n is not a power of 2, some chunks have not the maximum length, but
  this does not disturb essentially the mechanism

  Distinction between complete and uncomplete groups => 25% saving
  NOTICE: stdlib has void qsort (void* base, size_t num, size_t size,
            int (*compar)(const void*,const void*));
*/
//###################################################

void Statistics::quicksort(const int x[], int nData, int indexOrdered[]) const{ 
  double xDouble[nData];
  for(int i=0; i<nData; i++){
    xDouble[i]=(double)x[i];
  }
  quicksort(xDouble, nData, indexOrdered);
}

  
void Statistics::quicksort(const double x[], int nData, int indexOrdered[]) const{ 

  
  int* iOld = new int[nData];
  int* iNew = new int[nData];
  

  // initialize index array with result of first sorting round

  for (int j = 0; j < (int) nData / 2; j++) {
    int i = 2 * j;
    if (x[i] < x[i + 1]){
      iNew[i] = i;
      iNew[i + 1] = i + 1;
    }
    else {
      iNew[i] = i + 1;
      iNew[i + 1] = i;
    }
  }

  if (nData % 2 == 1){
    iNew[nData - 1] = nData - 1;
  }



  // determine number of sorting rounds:
  // n=7 => pmax=3, n=8 => pmax=3; n=9 => pmax=4

  int pmax = (int) (log(nData - 0.5) / log(2)) + 1;


  // main sorting routine (p=size of NEW group)

  int lgOld = 1; // length of old group
  int lg = 2; // length of new group

  for (int p = 2; p <= pmax; p++){
    //cout <<"p="<<p<<" pmax="<<pmax<<endl;
  

    // set up new iteration for  group size 2^p

    for(int i=0; i<nData; i++){
        iOld[i]= iNew[i];  // result iNew of last iteration = new iOld
    }


    lgOld = lg;
    lg *= 2;
    int ng = (nData - 1) / lg + 1; // number of new groups


    // iterate over groups of size 2^p

    for (int ig = 0; ig < ng; ig++){
      int m1 = ig * lg;
      int m2 = ig * lg + lgOld;
      int i1 = 0;
      int i2 = 0;
      int j1 = m1; // in inner loop: m1+i1
      int j2 = m2; // in inner loop: m2+i2
      int i = m1;
      bool bothGroupsComplete = (ig < ng - 1);
      int imax = min(m1 + lg, nData);


      // ################## bothGroupsComplete ##############

      if (bothGroupsComplete){

	for (i = m1; (i < imax) && (i1 < lgOld) && (i2 < lgOld); i++) {
	  if ((x[iOld[j1]] <= x[iOld[j2]])){
	    iNew[i] = iOld[j1];
	    i1++;
	    j1++;
	  } 
	  else {
	    iNew[i] = iOld[j2];
	    i2++;
	    j2++;
	  }
	}
      
	// treat the rest which is already sorted

	if (i1 >= lgOld){
	  int jmax = m1 + lg - i;
	  for (int j = 0; j < jmax; j++){

	    iNew[i + j] = iOld[j2 + j];
                                                }
	}

	if (i2 >= lgOld){
 
	  int jmax = m1 + lg - i; // NOT m2!
	  for (int j = 0; j < jmax; j++){
	    iNew[i + j] = iOld[j1 + j];
	  }
                                        }
      } // if bothGroups complete
    

      // ################## at least one groups not complete
      else{
	for (i = m1; (i < imax) && (i1 < lgOld) && (i2 < lgOld) && (j1 < nData)
					       && (j2 < nData); i++) {

	  if ((x[iOld[j1]] <= x[iOld[j2]])){
	    iNew[i] = iOld[j1];
	    i1++;
	    j1++;
	  } 
	  else {
	    iNew[i] = iOld[j2];
	    i2++;
	    j2++;
	  }
	}

	// treat the rest which is already sorted

	if (i1 >= lgOld){
	  int jmax = min(m1 + lg, nData) - i;
	  for (int j = 0; j < jmax; j++){
	    iNew[i + j] = iOld[j2 + j];

	  }
	}
	if (i2 >= lgOld){
	  int jmax = min(m1 + lg, nData) - i; // NOT m2!
	  for (int j = 0; j < jmax; j++){ 

	    iNew[i + j] = iOld[j1 + j];
	  }
	}
 
	if (j2 >= nData){
	  for (int j = 0; j < nData - i; j++){
	    iNew[i + j] = iOld[j1 + j];
	  }
	}

      } // end treatment of not complete groups

    } // end sorting through elements of groups at level p

    if(false){
      for (int i=0; i<nData; i++){
	cout <<"p="<<p<<" i="<<i<<" x[iNew[i]]="<<x[iNew[i]]<<endl;
      }
    }
  }

 // end iteration through levels p

  // output

  for (int i=0; i<nData; i++){indexOrdered[i]=iNew[i];}
  
  delete[] iOld; 
  delete[] iNew;
  
} // end method quicksort



////////////////////////////////////////////////////////////////////////
/// bubble sort algorithm
/// <ul>
/// <li> Input: array data and length of array,
/// <li> optional bool argument: true=ascending, false=descending order (default is true)
/// <li> output: array in sorted order
/// <li> Remark: QUICK SORT IS FASTER, IMPLEMENT QUICKSORT WHEN NECESSARY!
/// </ul>
////////////////////////////////////////////////////////////////////////

void Statistics::bubble_sort(double array[], int ndata, bool ascending)
{
  int i, j, flag = 1;    // set flag to 1 to begin initial pass
  double temp;           // holding variable
  for(i=1; (i<ndata) && flag; i++)
    {
      flag = 0;
      for(j=0; j < (ndata-1); j++)
	{
	  if( (ascending && array[j+1] < array[j]) 
	      ||  (!ascending && array[j+1] > array[j]))      // ascending order simply changes to <
	    { 
	      temp       = array[j];      // swap elements
	      array[j]   = array[j+1];
	      array[j+1] = temp;
	      flag       = 1;            // indicates that a swap occurred.
	    }
	}
    }
  return;
}


////////////////////////////////////////////////////////////////////////
/// bubble sort algorithm for 2d data (two arrays)
/// which are sorted according first array
/// <ul>
/// <li> Input: array data, secondary data and length of array,
/// <li> optional bool argument: true=ascending, false=descending order (default is true)
/// <li> output: array in sorted order
/// <li> Remark: QUICK SORT IS FASTER, IMPLEMENT QUICKSORT WHEN NECESSARY!
/// </ul>
////////////////////////////////////////////////////////////////////////

void Statistics::bubble_sort_2d(double array[], double array2[], int ndata, bool ascending)
{
  int i, j, flag = 1;    // set flag to 1 to begin initial pass
  double temp, temp2;           // holding variable
  for(i=1; (i<ndata) && flag; i++)
    {
      flag = 0;
      for(j=0; j < (ndata-1); j++)
	{
	  if( (ascending && array[j+1] < array[j]) 
	      ||  (!ascending && array[j+1] > array[j]))      // ascending order simply changes to <
	    { 
	      temp        = array[j];      // swap elements
	      array[j]    = array[j+1];
	      array[j+1]  = temp;
	      temp2       = array2[j];    // swap second array, too
	      array2[j]   = array2[j+1];
	      array2[j+1] = temp2;
	      flag        = 1;            // indicates that a swap occurred.
	    }
	}
    }
  return;
}




////////////////////////////////////////////////////////////////////////
/// histogram
/// arne: have changed datatype of histogram[] to integer! (6-10-05)
////////////////////////////////////////////////////////////////////////

void  Statistics::calculate_histogram(const double data[], int ndata,
				      double classMin, double classWidth, 
				      int nClass, int histogram[]){
  for(int k=0; k<nClass; k++){histogram[k]=0;} //!! need to init to zero!!
  for (int i=0; i<ndata; i++){
      //int k=static_cast<int>((data[i]-classMin)/(classWidth*(1.00001)));
      // changed by Christian: removed the 1.00001 factor (also in the
      // functions below)
      int k=static_cast<int>((data[i]-classMin)/classWidth);
      
      //cout<<"statistics::calculate_histogram: k="<<k<<endl;
      
      if((k>=0) && (k<nClass)){
	  histogram[k]++;
      }
  }
}


// MT 2023-09: vector argument (double for compatibility write_array)
// NEW: With extreme values gathered at the two ends
void  Statistics::calculate_histogram(const vector<double> data,
				      double classMin, double classWidth, 
				      int nClass, double histogram[]){
  for(int k=0; k<nClass; k++){histogram[k]=0;} //!! need to init to zero!!

  for (int i=0; i<int(data.size()); i++){
    int k=static_cast<int>((data[i]-classMin)/classWidth);

    // !! MT 2023-09 values beyond boundaries added to the boundaries

    int keff=(k<0) ? 0 : (k>nClass-1) ? nClass-1 : k;
    histogram[keff]++;
  }
}


void  Statistics::calculate_histogram(const double data[], int ndata,
				      int imin, int nInterval,
				      double classMin, double classWidth, 
				      int nClass, int histogram[]){
  for(int k=0; k<nClass; k++){histogram[k]=0;} //!! need to init to zero!!
  for (int i=imin; i<min(ndata, imin+nInterval); i++){

      int k=static_cast<int>((data[i]-classMin)/classWidth);
      
      //cout<<"statistics::calculate_histogram: k="<<k<<endl;
      
      if((k>=0) && (k<nClass)){
	  histogram[k]++;
      }
  }
}

void  Statistics::calculate_histogram(const int data[], int ndata,
				      double classMin, double classWidth, 
				      int nClass, int histogram[]){
  for(int k=0; k<nClass; k++){histogram[k]=0;} //!! need to init to zero!!
  for (int i=0; i<ndata; i++){
      //int k=static_cast<int>((data[i]-classMin)/(classWidth*(1.00001)));
      // changed by Christian: removed the 1.00001 factor (also in the
      // functions below)
      int k=static_cast<int>((data[i]-classMin)/classWidth);
      
      //cout<<"statistics::calculate_histogram: k="<<k<<endl;
      
      if((k>=0) && (k<nClass)){
	  histogram[k]++;
      }
  }
}

////////////////////////////////////////////////////////////////////////
/// returned array histogram[] is of type double !!!
/// provided to be compatible with other versions
////////////////////////////////////////////////////////////////////////

void  Statistics::calculate_histogram(const double data[], int ndata,
				      double classMin, double classWidth, 
				      int nClass, double histogram[])
{
  for(int k=0; k<nClass; k++){histogram[k]=0;} //!! need to init to zero!!
  for (int i=0; i<ndata; i++)
    {
      int k=static_cast<int>((data[i]-classMin)/classWidth);
      
      //cout<<"statistics::calculate_histogram: k="<<k<<endl;
      
      if((k>=0) && (k<nClass))
	{
	  histogram[k]++;
	}
    }
}


////////////////////////////////////////////////////////////////////////
/// 2d histogram
/// Added by Christian, March 29, 2007
////////////////////////////////////////////////////////////////////////

void  Statistics::calculate_histogram_2d(
    const double data1[], const double data2[], int ndata,
    double classMin1, double classWidth1, int nClass1,
    double classMin2, double classWidth2, int nClass2,
    int histogram[]){

  for(int k=0; k<nClass1*nClass2; k++){histogram[k]=0;} //!! need to init to zero!!


  for (int i=0; i<ndata; i++){
      int k1 = (int)floor((data1[i]-classMin1)/classWidth1);
      int k2 = (int)floor((data2[i]-classMin2)/classWidth2);
      
      if ((k1 >= 0) && (k2 >= 0) && (k1 < nClass1) && (k2 < nClass2))
	{
	  histogram[k1*nClass2+k2]++;
	}
    }
}

////////////////////////////////////////////////////////////////////////
/// 2d histogram with conditional means output
/// Added by Christian, July 21, 2007
////////////////////////////////////////////////////////////////////////

void  Statistics::calculate_histogram_2d(
    const double data1[], const double data2[], int ndata,
    double classMin1, double classWidth1, int nClass1,
    double classMin2, double classWidth2, int nClass2,
    int histogram[],
    double means1[], double means2[],
    int sums1[], int sums2[])
{
  for (int i=0; i<ndata; i++)
    {
      int k1 = (int)floor((data1[i]-classMin1)/classWidth1);
      int k2 = (int)floor((data2[i]-classMin2)/classWidth2);
      
      if ((k1 >= 0) && (k2 >= 0) && (k1 < nClass1) && (k2 < nClass2))
	{
	  histogram[k1*nClass2+k2]++;
          sums1[k1]++;
          means1[k1] += data2[i];
          sums2[k2]++;
          means2[k2] += data1[i];
	}
    }
  for (int k1 = 0; k1 < nClass1; k1++)
    means1[k1] /= sums1[k1];
  for (int k2 = 0; k2 < nClass2; k2++)
    means2[k2] /= sums2[k2];
}



//////////////////////////////////////////////////////////////////////////////

//##################################################
// Cumulative standard normal distribution
// from http://www.sitmo.com/doc/Calculating_the_Cumulative_Normal_Distribution#C.2B.2B_code
// added by Treibi, Jan2010
//##################################################

double Statistics::Phi(const double x){
  const double b1 =  0.319381530;
  const double b2 = -0.356563782;
  const double b3 =  1.781477937;
  const double b4 = -1.821255978;
  const double b5 =  1.330274429;
  const double p  =  0.2316419;
  const double c  =  0.39894228;

  if(x >= 0.0) {
      double t = 1.0 / ( 1.0 + p * x );
      return (1.0 - c * exp( -x * x / 2.0 ) * t *
      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
  }
  else {
      double t = 1.0 / ( 1.0 - p * x );
      return ( c * exp( -x * x / 2.0 ) * t *
      ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }
}

