#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// c++ 
#include <iostream>
#include <sstream>
//#include <fstream>
//#include <istream>

using namespace std;

 
// own

#define SQR(x) ((x)*(x))
#include "Math.h"

//##################################################
// Densities, distribution functions and inverse (qualtile) functions
// of standard distributions =>
// StatisticFunctions.cpp, .h
//##################################################

//#######################################
/// Cumulative standard normal distribution
//#######################################

// http://www.richelbilderbeek.nl/CppGetCumulativeDensityNormal.htm

double Math::norm(const double x)
{
  const double c0 = 0.2316419;
  const double c1 = 1.330274429;
  const double c2 = 1.821255978;
  const double c3 = 1.781477937;
  const double c4 = 0.356563782;
  const double c5 = 0.319381530;
  const double c6 = 0.398942280401;
  const double negative = (x < 0 ? 1.0 : 0.0);
  const double xPos = (x < 0.0 ? -x : x);
  const double k = 1.0 / ( 1.0 + (c0 * xPos));
  const double y1 = (((((((c1*k-c2)*k)+c3)*k)-c4)*k)+c5)*k;
  const double y2 = 1.0 - (c6*exp(-0.5*xPos*xPos)*y1);
  return ((1.0-negative)*y2) + (negative*(1.0-y2));
}




//#######################################
/// Digital Fourier transformation (computing time not relevant, therefore not FFT)
//#######################################

void Math::calcDFT(const double input_re[NMAX], const double input_im[NMAX], const int n,
		   double output_re[NMAX], double output_im[NMAX], bool inverse){
  if(n>NMAX){
    cerr<<"Error: more datapoints than Math.NMAX"<<endl; exit(-1);
  }
  
  double prefactor=(inverse) ? 1 : 1./n;   // M_PI = pi from c's math.h
  for (int j=0; j<n; j++){
    output_re[j]=0;
    for (int k=0; k<n; k++){
       double arg=2*M_PI*j*k/n;
       double exp_re=cos(arg);
       double exp_im=(inverse) ? sin(arg) : -sin(arg);
       output_re[j] += exp_re*input_re[k]  - exp_im*input_im[k];
       output_im[j] += exp_re*input_im[k] + exp_im*input_re[k];
    }
    output_re[j] *=prefactor;
    output_im[j] *=prefactor;
    // if(j>n-3)cout <<"in Math.calcDFT: j="<<j<<" input_im[0]= "<<input_im[0]<<endl;
  }
}

  
//#######################################
/// Inverse of regular quadratic matrix (elements a[0][] ... a[n-1][n-1])
//#######################################

void Math::calcMatrixInverse(const double matrix[NMATRIX][NMATRIX], const int n,
			     double inv[NMATRIX][NMATRIX], bool test){

  if(n>NMAX){
    cerr<<"Math.calcMatrixInverse: error:\nToo big matrices to invert: n="<<n
	<<" larger than NMAX="<<NMAX<<" increase NMAX\n";
    exit(-1);
  }
  
  // Initialisieren der Inversen als Einheitsmatrix

  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      inv[i][j]=0;
    }
    inv[i][i]=1;
  }


  // Kopieren von matrix in temporaere Matrix a (diese wird veraendert!)

  double a[NMATRIX][NMATRIX];
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      //cout <<"i="<<i<<" j="<<j<<" matrix[i][j]="<<matrix[i][j]<<endl;
      a[i][j]=matrix[i][j];
      //cout <<"i="<<i<<" j="<<j<<" a[i][j]="<<a[i][j]<<endl;
      
    }
  }


  //cout <<"hier4"<<endl;
  //cout <<"a[4][4]="<<a[4][4]<<endl;
 
  // Invertieren durch Zeilenmanipulationen an der Gleichung
  // A.X=I   => I.X=X (lhs) = manip(I)=inv (rhs)
  // Das "X" brauch explizit nicht eingefuehrt werden
  // (kein Zeilentausch, um kleine Pivotelemente zu verhindern!!)

  
  for (int k=0; k<n; k++){
    //cout <<"a[4][4]="<<a[4][4]<<endl;

     double piv=a[k][k];
     if(fabs(piv)<1.e-10){
       cerr <<"Error: diagonal element a["<<k<<"]["<<k<<"]="<<a[k][k]<<
  	 " of matrix to be inverted too small!\n";
        exit(-1);
     }
     //cout <<"k="<<k<<" pivot element piv="<<piv<<endl;


     for (int i=0; i<n; i++){
       

       // Manipulieren der linken Seite a_{ij}

       if(i != k) for (int j=k+1; j<n; j++){
	a[i][j]-=a[i][k]/piv*a[k][j];
       }

      // Manipulieren der rechten Seite inv_{ij}
       
       if(i != k) for (int j=0; j<=k; j++){
	inv[i][j]-=a[i][k]/piv*inv[k][j];
       }
     }
  }

  // Dividieren der rechten Zeilen  durch die verbleibenden
  // Diagonalelemente links:

  for (int k=0; k<n; k++){
     double piv=a[k][k];
     if(fabs(piv)<1.e-10){
        cerr <<"Error: diagonal element "<<(k+1)<<
  	 " of matrix to be inverted too small!\n";
        exit(-1);
     }
     //cout <<"final diagonal division: k="<<k<<" pivot element piv="<<piv<<endl;
     for (int j=0; j<n; j++){
       inv[k][j]/=piv;
     }
  }

 if(test){
   test_inversion(matrix,n,inv);
 }
  
}

void Math::test_inversion(const double matrix[NMATRIX][NMATRIX], const int n,
			  const double inv[NMATRIX][NMATRIX]) const {
  
   cout <<"in Math.test_inversion: Testing inverse of the provided matrices ..."
	<<endl;
   for (int l=0; l<n; l++){
      for (int m=0; m<n; m++){
	cout << "inv["<<l<<"]["<<m<<"]="<<inv[l][m]<<endl;
      }
   }
   
  // test by making product (Inverted matrix).matrix
   
   double prod1[NMATRIX][NMATRIX];

   for (int i=0; i<n; i++){
     for (int j=0; j<n; j++){
      prod1[i][j]=0;
      for (int k=0; k<n; k++){
	prod1[i][j] +=inv[i][k]*matrix[k][j];
      }
     }
   }

   cout <<"Test product (Inverted matrix).matrix :\n";
   for (int i=0; i<n; i++){
     for (int j=0; j<n; j++){
      cout <<"prod1["<<i<<","<<j<<"] = "<<prod1[i][j]<<endl;
     }
   }
} // end test_inversion





//#######################################
/// general-purpose operations
//#######################################


double Math::getmax(const double data[], int ndata)const
{
  double x=data[0];
  for(int i=1; i<ndata; i++){
    if(data[i]>x){x=data[i];}
  }
  return x;
}

double Math::getmin(const double data[], int ndata)const
{
  double x=data[0];
  for(int i=1; i<ndata; i++){
    if(data[i]<x){x=data[i];}
  }
  return x;
}

int Math::getmax(const int data[], int ndata)const
{
  int x=data[0];
  for(int i=1; i<ndata; i++){
    if(data[i]>x){x=data[i];}
  }
  return x;
}

int Math::getmin(const int data[], int ndata)const
{
  int x=data[0];
  for(int i=1; i<ndata; i++){
    if(data[i]<x){x=data[i];}
  }
  return x;
}

//###################################################


