/*
 * -------------------------------------------------------------------
 *
 * Copyright 2004,2005,2006,2007 Anthony Brockwell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the Free
 * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * -------------------------------------------------------------------
 *
 *  Here are declarations of matrix and vector
 *  classes.  Started Sep. 1992, by Anthony Brockwell.
 *
 *  Use of LAPACK routines incorporated around late 2002.
 */

#ifndef alreadymatrix
#define alreadymatrix

#include <complex.h>
#include <iostream.h>

#define complex std::complex<double>

namespace mslib {

const int MATLAB_MODE=0, BINARY_MODE=1, DISP_MODE=2, ROW_MODE=3; // used for I/O

class Double_Store {
 protected:
   double *stuff;
   int num_refs;
 public:
   Double_Store(int len);
   ~Double_Store();
   double& operator[](int i)
     {return stuff[i];}
};

class Vector;

class Matrix {
protected:
  int xsize, ysize, maxels;
  double *stuff;
  int iomode;

public:
  Matrix();         // no param. constructor
  Matrix(const Matrix&); // copy constructor
  Matrix(int,int,double **dp=NULL);  // constructors
  ~Matrix();        // destructor

  Matrix& operator=(const Matrix&);
  Matrix& operator=(double **dp);
  Matrix& operator=(const double);
  Matrix& operator+=(const Matrix&);
  Matrix& operator-=(const Matrix&);
  Matrix& operator+=(double other);
  Matrix& operator-=(double other);
  Matrix operator+(Matrix);
  //friend const Matrix operator*(Matrix& left, Matrix& right);
  Matrix operator-(Matrix);
  Matrix operator*(Matrix);
  Vector operator*(Vector);
  Matrix operator*(double);
  Matrix operator/(double);
  Matrix& operator*=(double);
  Matrix& operator/=(double);
  bool operator<=(const Matrix&);	// elementwise comparison
  inline double *operator[](int i) const
    { return (&stuff[i*xsize]); }
  Vector extract_row(int rownum);
  Vector extract_column(int colnum);

  Matrix LU_decomp(int *retsign=NULL, Matrix *permutations=NULL,
		  int *indices=NULL);
  Matrix inverse();
  Matrix g4_inverse(double epsilon);

  double symm_cond_num();   // for symmetric A, returns condition number
  void solve_system(Vector *b);     // solve Ax=b, (ret. soln in b)
  void symm_solve_system(Vector *b); // same but for symmetric A
  void robust_symm_solve_system(Vector *b);  // handles near-singular matrices

  Matrix transpose();
  Matrix submatrix(int r1, int c1, int r2p1, int c2p1);
  bool singular(double *);	// returns determinant and flag
  int nrows() { return ysize; };
  int ncols() { return xsize; };
  double *get_stuff() { return stuff; };
  complex *eigenvalues();
  Matrix symm_eigenvectors(double *evals);  // a crude routine: LAPACK dsyev should be used when possible
  double norm();
  double max();
  double min();
  Vector mean();        // returns mean down the columns
  Matrix covariance();  // returns sample covariance (each row is a datapoint)
  double kernel_function(double arg, int n, double sd);  // used by quantiles function below
  double kernel_cdf(double arg, int n, double sd);       // also used by quantiles f-n
  Vector quantiles(double alpha, bool kernel_smooth=false); 
                        // returns vector of quantiles down columns,
                        // alpha is between 0 and 1

  void identity();
  void zeroes();
  void resize(int,int);
  void set_data(double *, int, int, int);
  double symm_quad_form(Vector& v);

  int read_file(const char *);      // read data set as col. vec. and return # points
  void write_file(const char *);    // write file
  void set_iomode(int im) {iomode=im;}
  void append_bottom(Matrix& newstuff, int trans=0);

  friend ostream& operator<<(ostream&, Matrix&);
  friend istream& operator>>(istream&, Matrix&);
};

// Vectors are implemented as a special case of matrices.

class Vector : public Matrix {
public:
  Vector();
  Vector(int);

  double& operator[](int i) {return stuff[i];}
  Vector& operator=(const Matrix&);
  Vector operator-(Vector);
  Vector operator*(const double&);
  double dot(const Vector&);
  double norm();
  double sum();
  double mean();
  double var();
  double angle();  // returns angle for 2-d vectors
  Matrix outer(const Vector&);
  int size() {return ysize;} 
  void resize(int nrows) {Matrix::resize(nrows, 1);}
  void append(double);
  Vector subvector(int r1, int r2);
  void sort();     // heapsort elements of vector

  friend ostream& operator<<(ostream&, Vector&);
  friend istream& operator>>(istream&, Vector&);
};

// Here are some miscellaneous functions

Matrix blockmatrix(Matrix& A, Matrix& B, Matrix& C, Matrix& D);
Matrix vblockmatrix(Matrix& A, Matrix& B);
Matrix hblockmatrix(Matrix& A, Matrix& B);

}

#endif
