////////////////////////////////////////////////////////////////////////////////////
//  Example program that shows how to use levmar in order to fit the three-
//  parameter exponential model x_i = p[0]*exp(-p[1]*i) + p[2] to a set of
//  data measurements; example is based on a similar one from GSL.
//
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

#include "levmar.h"

#ifndef LM_DBL_PREC
#error Example program assumes that levmar has been compiled with double precision, see LM_DBL_PREC!
#endif


/* the following macros concern the initialization of a random number generator for adding noise */
#undef REPEATABLE_RANDOM
#define DBL_RAND_MAX (double)(RAND_MAX)

#ifdef _MSC_VER // MSVC
#include <process.h>
#define GETPID  _getpid
#elif defined(__GNUC__) // GCC
#include <sys/types.h>
#include <unistd.h>
#define GETPID  getpid
#else
#warning Do not know the name of the function returning the process id for your OS/compiler combination
#define GETPID  0
#endif /* _MSC_VER */

#ifdef REPEATABLE_RANDOM
#define INIT_RANDOM(seed) srandom(seed)
#else
#define INIT_RANDOM(seed) srandom((int)GETPID()) // seed unused
#endif

/* Gaussian noise with mean m and variance s, uses the Box-Muller transformation */
double gNoise(double m, double s)
{
double r1, r2, val;

  r1=((double)random())/DBL_RAND_MAX;
  r2=((double)random())/DBL_RAND_MAX;

  val=sqrt(-2.0*log(r1))*cos(2.0*M_PI*r2);

  val=s*val+m;

  return val;
}

/* model to be fitted to measurements: haty_i = p[0]*exp(-p[1]*i) + p[2]+...,  i=0...n-1 */
void ymodelfunc(double *p, double *haty, int m, int n, double *data)
{
register int i;

  for(i=0; i<n; ++i){
    //haty[i]=p[0]*exp(-p[1]*i) + p[2];  //m=3
    haty[i]=p[0]*exp(-p[1]*i) + p[2] + p[3]/(i*i+1);  //m=4
  }
}

/* Jacobian of ymodelfunc() */
void jacymodelfunc(double *p, double *jac, int m, int n, double *data)
{   
register int i, j;
  
  /* fill Jacobian row by row */
  for(i=j=0; i<n; ++i){
    jac[j++]=exp(-p[1]*i);
    jac[j++]=-p[0]*i*exp(-p[1]*i);
    jac[j++]=1.0;
  }
}

//#############################################
int main()
{
  const int n=1000; //n measurements 
  const int m=4; // m parameters
  double p[m]; // parameter vector
  double ydata[n]; // data vector
  double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
  const int n_iter_max=1000;

  register int i;
  int ret; // return signal of the dlevmar functions

  /* generate some measurement using the exponential model with
   * parameters (5.0, 0.1, 1.0), corrupted with zero-mean
   * Gaussian noise of s=0.1
   */
  INIT_RANDOM(0);
  for(i=0; i<n; ++i)
    ydata[i]=(5.0*exp(-0.01*i) + 1.0) + exp(-0.01*i)*gNoise(0, 1);

  /* initial parameters estimate: (1.0, 0.0, 0.0) */
  //p[0]=1.0; p[1]=0.0; p[2]=0.0;  //m=3
  p[0]=0.0; p[1]=0.0; p[2]=0.0; p[3]=0;  //m=4

  /* optimization control parameters; passing to levmar NULL instead of opts reverts to defaults */
  opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
  opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used 

  double covar[m*m];


  /* invoke the optimization function */
  //ret=dlevmar_der(ymodelfunc, jacymodelfunc, p, ydata, m, n, n_iter_max, opts, info, NULL, NULL, NULL); // with analytic Jacobian
  ret=dlevmar_dif(ymodelfunc, p, ydata, m, n, n_iter_max, opts, info, NULL, covar, NULL); // without Jacobian



  printf("Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [initially %g]\n", info[5], info[6], info[1], info[0]);
  printf("  reason 1: small gradient;  2: small param change; 3: itmax reached;\n");
  printf("         4: singular matrix; 7: NaN or Inf\n");
  printf(" number of function evaluations: %g\n", info[7]);
  printf(" number of Jakobian evaluations: %g\n", info[8]);
  printf("Best fit parameters: ");
  int im=0; for (im=0; im<m; im++){printf(" %.7g", p[im]);} printf("\n");
  printf("Parameter covariance Matrix:\n");
  //int im2=0;
  for (im=0; im<m; im++){
    for (int im2=0; im2<m; im2++){
      printf(" %.7g\t", covar[m*im+im2]);
    }
    printf("\n");
  }
  exit(0);
}
