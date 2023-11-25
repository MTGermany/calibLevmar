#ifndef MATH_H
#define MATH_H

/*! \class Math
 *  \brief Small collection of some numerical helpers
 *  \author Martin Treiber and Arne Kesting
 */

//##################################################
// (feb17) Densities, distribution functions and inverse (quantile) functions
// of standard distributions =>
// StatisticalFunctions.cpp, .h
//##################################################

class Math{ 

 public:

  Math(){;}
  static const int NMAX=2000; // greater than 1440! bei 5000 spinnt Notebook!
  static const int NMATRIX=100; // Notebook spinnt bei mehr als etwa 1 MB Daten

  double norm(const double x); // cumulative standard normal distribution


  void calcMatrixInverse(const double matrix[NMATRIX][NMATRIX], const int n,
			 double inv[NMATRIX][NMATRIX], bool test);
  void calcDFT(const double input_re[NMAX],
	       const double input_im[NMAX], const int n, 
	       double output_re[NMAX], double output_im[NMAX], bool inverse);
 
  void test_inversion(const double matrix[NMATRIX][NMATRIX], const int n,
		      const double inv[NMATRIX][NMATRIX]) const;
  double getmax(const double data[], int ndata) const;
  double getmin(const double data[], int ndata) const;
  int getmax(const int data[], int ndata) const;
  int getmin(const int data[], int ndata) const;

 
};

#endif // MATH_H
