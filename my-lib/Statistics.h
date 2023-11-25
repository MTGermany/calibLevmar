#ifndef STATISTICS_H
#define STATISTICS_H

#include <vector>

/*! \class Statistics
 *  \brief Small collection of some numerical helpers
 *  \author Martin Treiber and Arne Kesting
 *  \version 1.0
 *  \date    2000-2005
 *  \bug (unknown)
 *  \warning (???)
 */

//##################################################
// (feb17) Densities, distribution functions and inverse (quantile) functions
// of standard distributions =>
// StatisticalFunctions.cpp, .h
//##################################################

class Statistics{ 

 public:

  Statistics(){;} 

  double sum(const double data[], int ndata) const;
  int sum(const int data[], int ndata) const;
  double arithmeticMeans(const double data[], int ndata) const;
  double arithmeticMeans(const int data[], int ndata) const;
  double arithmeticMeans(const double data[], int imin, int imax) const;
  double getmax(const double data[], int ndata) const;
  double getmin(const double data[], int ndata) const;
  int getmax(const int data[], int ndata) const;
  int getmin(const int data[], int ndata) const;
  double variance(const double data[], int ndata) const;
  double variance(const double data[], int imin, int imax) const;
  double varianceDescr(const double data[], int ndata) const;
  double varianceDescr(const double data[], int imin, int imax) const;
  void distrF(const double x[], int nData,
	      double xSorted[], double distr[]) const;
  double quantile(const double data[], int ndata, double q) const;
  int quantile(const int data[], int ndata, double q) const;
  int median(const int data[], int ndata){return quantile(data,ndata,0.5);}
  double covariance(const double xdata[], const double ydata[], int ndata);
  double correlation(const double xdata[], const double ydata[], int ndata);

  void calcCorrelationFun(const double xdata[], int ndata, 
			  int nshiftmax, double corrfun[]);

  void calcCrossCorrelationFun(const double xdata[], const double ydata[], int ndata,
			  int nshiftmax, double corrfun[]);

  void calcCorrelationFunDetrend(const double xdata[], int ndata,
				 int nshiftmax, int periodMA,
				 double corrfun[], bool& success);

  void linRegression(const double xdata[], const double ydata[], int ndata,
		     double& a, double& b, 
		     double& B, double& sumdev);

  void linRegression(const double ydata[], int ndata,
		     double& a, double& b, 
		     double& B, double& sumdev);

  void linRegression(const double xdata[], const double ydata[], 
		     int ndata, int imin, int imax, 
		     double& a, double& b, 
		     double& B, double& sumdev);


   // cumulative standard normal distribution
  double Phi(const double x);

  void smoothGauss(const double data[], int ndata, double dnSmooth,
		   double smoothedData[]);

  void gSmoothGauss(const double xdata[],
		    const double ydata[], 
		    const int ndata, 
		    const double sigma,
		    double smoothedData[]);

  
  // smooth with triang kernel (other interface!) [MT 2022-09]
  vector<double> smoothTriang(const vector<double> data, int dnSmooth);
  vector<double> smoothTriang(const vector<int> data, int dnSmooth);
  

  // time series
  void calcEMA(const double data[], int ndata, double initVal,
	       double tau, double movEAvg[]);

  void calcEMA(const double data[], int ndata,
	       double tau, double movEAvg[]);

  void calcSymmetricEMA(const double data[], int ndata,
			double dnSmooth, double smoothedData[]);

  void calcSymmetricEMA(const double tdata[], const double data[], int ndata,
			double dtSmooth, double dtOut,
			double tdataOut[], double dataOut[]);

  void calcMA(const double data[], int ndata,
	      int period, double movAvg[]);

  void calculate_detr_saison(const double data[], int ndata,
			     int period, int nPeriods, int istart,
			     double data_saison[]);

  void calculate_saison(const double data[], int ndata,
			int period, double data_saison[]);

   double gGaussAverage(const double xdata[],
		       const double ydata[],
		       const int ndata, 
		       const double xCenter,
		       const double sigma);

  void gaussLinRegression(const double xdata[], const double ydata[], 
			  int ndata, 
			  double xCenter, double sigma,
			  double& a, double& b, 
			  double& B, double& sumdev);


  void bubble_sort(double array[], int ndata, bool ascending=true);

  void bubble_sort_2d(double array[], double array2[], int ndata, bool ascending=true);

  void quicksort(const double x[], int nData, int indexOrdered[]) const;
  void quicksort(const int x[], int nData, int indexOrdered[]) const;


  void calculate_histogram(const double data[], int ndata,
                           double classMin, double classWidth, int nClass, 
                           int histogram[]);

  // MT 2023-09: vector argument (double for compatibility write_array)
  void calculate_histogram(const vector<double> data, 
                           double classMin, double classWidth, int nClass, 
                           double histogram[]); 

  void calculate_histogram(const double data[], int ndata,
			   int imin, int nInterval,
                           double classMin, double classWidth, int nClass, 
                           int histogram[]); 

  void calculate_histogram(const int data[], int ndata,
                           double classMin, double classWidth, int nClass, 
                           int histogram[]); 

  void calculate_histogram(const double data[], int ndata,
                           double classMin, double classWidth, int nClass, 
                           double histogram[]);

  void calculate_histogram_2d(const double data1[], const double data2[], int ndata,
                              double classMin1, double classWidth1, int nClass1, 
                              double classMin2, double classWidth2, int nClass2, 
                              int histogram[]);

  void calculate_histogram_2d(const double data1[], const double data2[], int ndata,
                              double classMin1, double classWidth1, int nClass1, 
                              double classMin2, double classWidth2, int nClass2, 
                              int histogram[],
                              double means1[], double means2[],
                              int sums1[], int sums2[]);


  ////////////////////////////////////////////////////////////////////////////

 private:

  
};

#endif // STATISTICS_H
