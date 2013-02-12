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
 *  This module defines functions for the "Matrix" class,
 *  which include the basic operations of matrix addition,
 *  subtraction, multiplication, as well as system solving,
 *  finding eigenvalues/vectors, etc.
 *
 *  (Originally this code was based largely on algorithms in the book
 *   "Numerical Recipes for ...", but now the algorithms have been
 *   largely replaced with more efficient algorithms, including
 *   some from the BLAS/LAPACK libraries.)
 *
 *  Author: Anthony Brockwell.
 *
 *  Initially created 1992.
 *
 *  Update Mar. 97 to fix oversight where LU decomposition
 *  did not use pivoting and to allow calculation of eigenvalues.
 *
 *  Updated May 2002 to use LAPACK, which is much
 *  faster than numerical recipes routines.  This change 
 *  probably also eliminated many bugs.
 *
 *  Use the -DUSE_LAPACK compile command-line option
 *  to use LAPACK (much faster).  Without this option, slow
 *  relatively inefficient routines will be used instead.
*/

#include <iostream>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <string>
#include "matrix.h"
#include "poly.h"		// used for root-finding
#include "heapsort.h"

using namespace mslib;

#define TINY 1e-20

#ifdef USE_LAPACK
extern "C" {
  void dgesv_(int *n, int *ncolsb, double *a, int *lda, int *ipiv,
	      double *b, int *ldb, int *info);
  void dsysv_(char *uplo, int *n, int *ncolsb, double *a, int *lda, int *ipiv,
	      double *b, int *ldb, double *work, int *lword, int *info);
  void dsymv_(char *uplo, int *n, double *alpha, double *a, int *lda, double *x,
	      int *incx, double *beta, double *y, int *incy);
  void dsyev_(char *job, char *uplo, int *n, double *a,
	      int *lda, double *evals, double *work, int *lwork,
	      int *output);                         // finds eigenvalues/eigenvectors of a symmetric matrix
  void dgemm_(char *, char *, int *m, int *n,
	      int *k, double *alpha, double *a,
	      int *lda, double *b, int *ldb, double *beta,
	      double *c, int *ldc);
  void dgemv_(char *, int *m, int *n, double *alpha, double *a, int *lda,
	      double *x, int *incx, double *beta, double *y, int *incy);
}
#endif


namespace mslib {


double const_one=1.0,const_zero=0.0;
int const_int_one=1;


Double_Store::Double_Store(int len)
{
  stuff = new double[len];
  num_refs = 1;
}

Double_Store::~Double_Store()
{
  delete[] stuff;
}


double complex_abs(complex& c)
{
  double tx = c.real()*c.real() + c.imag()*c.imag();
  return sqrt(tx);
}

// Constructors First
Matrix::Matrix()    // constructor if no parameters (init. to 1 by 1)
{ 
  xsize = 1; ysize = 1;
  iomode = MATLAB_MODE;
  maxels = xsize*ysize;
  stuff = (double *) new double[maxels];
}

Matrix::Matrix(int y,int x,double **dp)
{
  iomode = MATLAB_MODE;
  xsize = x; ysize = y;
  maxels = x*y;
  stuff = new double[maxels];
  if (dp!=NULL)
    (*this)=dp;
}

Matrix::~Matrix()
{ delete[] stuff; }

// Then Member Functions

bool Matrix::singular(double *dp)
  // see Num. Methods p. 43
{
  int i,sign;
  double dx=1.0;
  Matrix temp=LU_decomp(&sign);
  for (i=0 ; i<xsize ; ++i)
    dx *= temp[i][i];
  dx *= sign;
  if (dp!=NULL)
    *dp=dx;
  return(fabs(dx)<TINY);
}

double Matrix::symm_quad_form(Vector& x)
{
  Vector temp(x.nrows());
  cout << "temp = " << temp << endl;
  dsymv_("U",&xsize,&const_one,stuff,&xsize,x.stuff,&const_int_one,&const_zero,temp.stuff,&const_int_one);
  cout << "then temp = " << temp << endl;
  double rval = x.dot(temp);
  return rval;
}

Matrix Matrix::LU_decomp(int *rsign, Matrix *perms, int *indices)
{
  int i,j,k,m,*index,sign=1;
  double biggest,dx;
  // this procedure returns the LU decomposition (in a single
  // Matrix - a(ii) = 1 as per Num. Recipes pg. 41
  Matrix temp(ysize, xsize);
  if (xsize!=ysize) printf ("LU decomposing non-square Matrix!\n");
  temp=*this;
  index=new int[xsize];
  for (i=0 ; i<xsize ; ++i)
    index[i] = i;
  
  for (j=0 ; j<xsize ; ++j)
    {
      // calculate the betas in column j
      for (i=0 ; i<j ; ++i) // solve for b(ij) first
	for (k=0 ; k<=i-1 ; ++k)
	  temp[i][j] -= temp[i][k]*temp[k][j];
      // calculate the alphas in column j
      biggest=0.0;   m=j;
      for (i=j ; i<ysize ; ++i) // solve for a(ij) next
	{
	  for (k=0 ; k<=j-1 ; ++k)
	    temp[i][j] -= temp[i][k]*temp[k][j];
	  if ((dx=fabs(temp[i][j]))>biggest) {
	    biggest=dx;
	    m=i;				// m is row with largest abs. value of alpha
	  }
	}
      // pivot: swap row j for row m, with m>j
      i=index[j];
      index[j]=index[m];
      index[m]=i;
      if (j!=m) {
	for (i=0 ; i<xsize ; ++i) {
	  dx=temp[j][i];
	  temp[j][i]=temp[m][i];
	  temp[m][i]=dx;
	}
	sign=-sign;
      }
      // do the division for the alphas
      for (i=j+1 ; i<ysize ; ++i)
	temp[i][j] /= temp[j][j];
    }
  if (rsign!=NULL)
    *rsign = sign;
  if (perms!=NULL)
    {
      Matrix tm(ysize, xsize);
      tm.zeroes();
      for (i=0 ; i<xsize ; ++i)
	tm[i][index[i]]=1.0;
      *perms = tm;
    }
  if (indices!=NULL)
    for (i=0 ; i<xsize ; ++i)
      indices[i] = index[i];
  delete[] index;
  return temp;
}

#ifdef USE_LAPACK

// The following routines use BLAS and LAPACK and are much more efficient
// than the crappy ones lower down which use numerical recipes algorithms.

void Matrix::solve_system(Vector *vb)
  // solves Ax=b, (*this == A)
  // result is stored in b again
{
  // Now it uses LAPACK!
  int info,n=xsize,one=1;
  int *ipiv=new int[n];

  Matrix acopy = this->transpose();  // acopy is same as a, stored in
    // column-major order
  dgesv_(&n, &one, acopy.get_stuff(), &n, ipiv, vb->get_stuff(), &n, &info);

  delete[] ipiv;
}

void Matrix::symm_solve_system(Vector *vb)
  // exactly like solvesystem but makes use of the fact that
  // the Matrix is symmetric
{
  int info,n=xsize,one=1,lwork=2048;
  int *ipiv=new int[n];
  char uplo='U';
  static double work[2048];

  Matrix acopy = *this;  // acopy is same as a
  dsysv_(&uplo, &n, &one, acopy.get_stuff(), &n, ipiv, vb->get_stuff(), &n, 
	 work, &lwork, &info);

  delete[] ipiv;
}

Matrix Matrix::g4_inverse(double epsilon)
{
  int i,j,d = nrows(),output;
  Matrix temp(d,d);
  int scratchsize = 8*ncols();
  double *scratch = new double[scratchsize];

  Vector evals(d);
  Matrix evecst,evecs;

  evecst = *this;
  
  dsyev_("V","U",&d,&evecst[0][0],&d,&evals[0],scratch,&scratchsize,
	 &output);
  evecs=evecst.transpose();

  for (j=0 ; j<d ; ++j)
    if (fabs(evals[j])>epsilon)
      evals[j]=1.0/evals[j];
    else
      evals[j]=0.0;

  // now construct inverse
  for (i=0 ; i<d ; ++i)
    for (j=0 ; j<d ; ++j)
      evecst[i][j] *= evals[i];
  temp = evecs*evecst;

  delete[] scratch;
  return temp;  
}

void Matrix::robust_symm_solve_system(Vector *vb)
{
  double tx;
  tx = symm_cond_num();
  if (tx<1e12)
    solve_system(vb);
  else
    {
      Matrix tm = g4_inverse(0.00001);
      *vb = tm*(*vb);
    }
}

double Matrix::symm_cond_num()
{
  // returns condition number of symmetric Matrix
  static double scratch[2048];
  int scratchsize=2048;
  int output, d=nrows();
  Vector evals(d);
  Matrix evecs(d,d);
  evecs = *this;
  dsyev_("N","U",&d,&evecs[0][0],&d,&evals[0],scratch,&scratchsize,
	 &output);
  return (fabs(evals[d-1]/evals[0]));
}

Matrix Matrix::inverse()
  {
    // Now it uses LAPACK!
    int info,n=xsize;
    int *ipiv=new int[n];
    Matrix b(n,n);
    b.identity();

    Matrix acopy = this->transpose();  // acopy is same as a, stored in
    // column-major order
    dgesv_(&n, &n, acopy.get_stuff(), &n, ipiv, b.get_stuff(), &n, &info);
    
    delete[] ipiv;

    return (b.transpose());
  }

#else

// The following routines are only compiled if USE_LAPACK is not defined.
// They are highly inefficient!!! Only use them if absolutely necessary!

void Matrix::solve_system(Vector *vb)
  // solves Ax=b, (*this == A)
  // result is stored in b again
{
  Matrix b = this->inverse();
  (*vb) = b*(*vb);
}

void Matrix::symm_solve_system(Vector *vb)
  // exactly like solvesystem but makes use of the fact that
  // the Matrix is symmetric
{
  solve_system(vb);
}

Matrix Matrix::g4_inverse(double epsilon)
{
  cout << "Warning: inverse() used instead of g4_inverse()." << endl;
  return inverse();
}

void Matrix::robust_symm_solve_system(Vector *vb)
{
  double tx;
  tx = symm_cond_num();
  if (tx<1e12)
    solve_system(vb);
  else
    {
      Matrix tm = g4_inverse(0.00001);
      *vb = tm*(*vb);
    }
}

double Matrix::symm_cond_num()
{
  // returns condition number of symmetric Matrix
  static double scratch[2048];
  int scratchsize=2048;
  int output, d=nrows();
  Matrix evecs(d,d);
  evecs = *this;
  // get eigenvalues (which will be real since it's symmetric)
  complex *evals = eigenvalues();
  double scn = (complex_abs(evals[0]))/(complex_abs(evals[d-1]));
  delete[] evals;
  return scn;
}

Matrix Matrix::inverse()
  {
    // returns inverse of *this
    int sign;
    Matrix perms,temp((*this).LU_decomp(&sign, &perms)),
           temp2(ysize,xsize), temp3(ysize,xsize);
    int i,j,k;

    // first find solns to Lx=b, b is [1 0 0 0]', [0 1 0 0]', etc.
    for (k=0 ; k<xsize ; ++k)    // do this many vector solns
      for (i=0 ; i<xsize ; ++i)  // find this many elements in each vector
        {
          temp2[i][k] = (i==k);
          for (j=0 ; j<i ; ++j)
            temp2[i][k] -= temp[i][j]*temp2[j][k];
        }

    // now find solns to Uy=c, c is result of previous stage
    for (k=0 ; k<xsize ; ++k)   // this many vector solns again
      for (i=xsize-1 ; i>=0 ; --i)
        {
          temp3[i][k] = temp2[i][k];
          for (j=i+1 ; j<xsize ; ++j)
            temp3[i][k] -= temp[i][j]*temp3[j][k];
          temp3[i][k] /= temp[i][i];
        }
    temp3 = temp3 * perms;
    return temp3;
  }

#endif

Matrix blockmatrix(Matrix& A, Matrix& B, Matrix& C, Matrix& D)
  // returns [[A,B];[C,D]]
{
  int i,j,ca = A.ncols(), ra = A.nrows(),
    cb = B.ncols(), rb = B.nrows(),
    cc = C.ncols(), rc = C.nrows(),
    cd = D.ncols(), rd = D.nrows();
  
  if ((ra!=rb) || (ca!=cc) || (cb!=cd) || (rc!=rd))
    cout << "Error: block Matrix is not right!" << endl;

  Matrix temp(ra+rc,ca+cb);
  for (i=0 ; i<ra ; ++i)
    {
      for (j=0 ; j<ca ; ++j)
	temp[i][j] = A[i][j];
      for (j=0 ; j<cb ; ++j)
	temp[i][j+ca] = B[i][j];
    }
  for (i=0 ; i<rc ; ++i)
    {
      for (j=0 ; j<cc ; ++j)
	temp[i+ra][j] = C[i][j];
      for (j=0 ; j<cd ; ++j)
	temp[i+ra][j+cc] = D[i][j];
    }

  return temp;
}

Matrix vblockmatrix(Matrix& A, Matrix& B)
  // returns [A;B]
{
  int i,j,ca = A.ncols(), ra = A.nrows(),
    cb = B.ncols(), rb = B.nrows();
  
  if (ca!=cb)
    cout << "Error: block Matrix is not right!" << endl;

  Matrix temp(ra+rb,ca);
  for (i=0 ; i<ra ; ++i)
    for (j=0 ; j<ca ; ++j)
      temp[i][j] = A[i][j];
  for (i=0 ; i<rb ; ++i)
    for (j=0 ; j<ca ; ++j)
      temp[i+ra][j] = B[i][j];

  return temp;
}

Matrix hblockmatrix(Matrix& A, Matrix& B)
  // returns [A,B]
{
  int i,j,ca = A.ncols(), ra = A.nrows(),
    cb = B.ncols(), rb = B.nrows();
  
  if (ra!=rb)
    cout << "Error: block Matrix is not right!" << endl;

  Matrix temp(ra,ca+cb);
  for (i=0 ; i<ra ; ++i)
    {
      for (j=0 ; j<ca ; ++j)
	temp[i][j] = A[i][j];
      for (j=0 ; j<cb ; ++j)
	temp[i][j+ca] = B[i][j];
    }

  return temp;
}

Matrix Matrix::transpose()
  {
    int i,j;
    Matrix temp(xsize,ysize);
    for (i=0 ; i<xsize ; ++i)
      for (j=0 ; j<ysize ; ++j)
	temp[i][j] = (*this)[j][i];
    return(temp);
  }

double Matrix::min()
{
  double minval=stuff[0];
  int i;
  for (i=1 ; i<xsize*ysize ; ++i)
    if (stuff[i]<minval)
      minval=stuff[i];
  return minval;
}

double Matrix::max()
{
  double maxval=stuff[0];
  int i;
  for (i=1 ; i<xsize*ysize ; ++i)
    if (stuff[i]>maxval)
      maxval=stuff[i];
  return maxval;
}

ostream& operator<<(ostream& os, Matrix& m)
  // stream it in Matlab format
{
  int i,j,vecmode=0;
  switch (m.iomode) {
  case MATLAB_MODE:
    os << "[";
    for (i=0 ; i<m.nrows() ; ++i)
      {
	os << "[";
	for (j=0 ; j<m.ncols() ; ++j)
	  {
	    os << m[i][j];
	    if (j<m.ncols()-1)
	      os << ",";
	  }
	os << "]";
	if (i<m.nrows()-1)
	  os << ";";
      }
    os << "]";
    break;
  case DISP_MODE:
    os << m.nrows() << "*" << m.ncols() << " Matrix = " << endl;
  case ROW_MODE:  // run on from DISP_MODE
    for (i=0 ; i<m.nrows() ; ++i)
      for (j=0 ; j<m.ncols() ; ++j)
        if (j<m.ncols()-1)
	  os << m[i][j] << " ";
        else
	  os << m[i][j] << endl;

  default:
    break;
  }
  return(os);
}

istream& operator>>(istream& is, Matrix& m)
  // read it in (Matlab format only)
{
  int i,j,k,done,maxels,rows,cols;
  double *localstuff,*newstuff;
  char c;
  
  localstuff = new double[256];
  maxels = 256;
  is >> ws;   // eat white space
  is.get(c);  // remove "["
  done=0;

  rows = 0;
  i = 0;
  while (!done)
    {
      ++rows;

      // read a row
      c = 0;
      while ((c!='[') && is)
	is.get(c);
      j=i;
      while ((c!=']') && is)
	{
	  is >> localstuff[i++];             // read a double
	  if (i==maxels)
	    {
	      // we need to expand memory
	      maxels *= 2;
	      newstuff = new double[maxels];
	      for (k=0 ; k<maxels/2 ; ++k)
		newstuff[k] = localstuff[k];
	      delete[] localstuff;
	      localstuff = newstuff;
	    }
	  is >> ws;                         // eat white space
	  is.get(c);   
	  if ((c!=',') && (c!=']'))
	    is.putback(c);
	}
      cols=i-j;
      is >> ws;
      is.get(c);
      done=((c==']') || (!is));
    }

  m.set_data(localstuff,rows,cols,maxels);

  return is;
}

void Matrix::set_data(double *newdata, int r, int c, int me)
{
  delete[] stuff;
  stuff = newdata;
  xsize = c;   ysize = r;   maxels = me;
}

void Matrix::identity()
  {
    int i,j;
    for (i=0 ; i<ysize ; ++i)
      for (j=0 ; j<xsize ; ++j)
         if (i==j) stuff[i*xsize+j]=1; else stuff[i*xsize+j]=0;
  }

void Matrix::zeroes()
{
  int i;
  for (i=0 ; i<xsize*ysize ; ++i)
    stuff[i]=0;
}

Matrix::Matrix(const Matrix& m1) :
  maxels(m1.maxels), xsize(m1.xsize), ysize(m1.ysize), iomode(m1.iomode)
  // copy constructor
  {
    int i;
    // Copy constructor
    stuff = new double[maxels];
    for (i=0 ; i<m1.xsize*m1.ysize ; ++i)
      stuff[i] = m1.stuff[i];
  }

bool Matrix::operator<=(const Matrix& other)
  {
    int i;
    bool islessorequal=true;
    for (i=0 ; i<xsize*ysize ; ++i)
      islessorequal &= (stuff[i] <= other.stuff[i]);
    return islessorequal;
  }

Matrix& Matrix::operator=(const Matrix& m1)
  {
    int i;
    maxels = m1.maxels;    xsize=m1.xsize;   ysize=m1.ysize;
    delete[] stuff;  // get rid of anything that was there
    stuff = new double[maxels];
    for (i=0 ; i<m1.xsize*m1.ysize ; ++i)
      stuff[i] = m1.stuff[i];
    iomode = m1.iomode;
    return *this;
  }

Matrix& Matrix::operator+=(const Matrix& other)
{
  int i, len=xsize*ysize;
  if ((xsize!=other.xsize) || (ysize!=other.ysize))
    { cout << "Cannot add matrices of different sizes!" << endl; }
  else
    for (i=0 ; i<len ; ++i)
      stuff[i] += other.stuff[i];
  return *this;
}

Matrix& Matrix::operator+=(double other)
{
  int i, len=xsize*ysize;
  for (i=0 ; i<len ; ++i)
    stuff[i] += other;
  return *this;
}

Matrix& Matrix::operator-=(double other)
{
  int i, len=xsize*ysize;
  for (i=0 ; i<len ; ++i)
    stuff[i] -= other;
  return *this;
}

Matrix& Matrix::operator-=(const Matrix& other)
{
  int i, len=xsize*ysize;
  if ((xsize!=other.xsize) || (ysize!=other.ysize))
    { cout << "Cannot subtract matrices of different sizes!" << endl; }
  else
    for (i=0 ; i<len ; ++i)
      stuff[i] -= other.stuff[i];
  return *this;
}

Matrix& Matrix::operator=(const double dx)
  {
    int i;
    for (i=0 ; i<xsize*ysize ; ++i)
      stuff[i] = dx;
    return *this;
  }

Matrix& Matrix::operator=(double **dp)
  {
    int i,j;
    for (i=0 ; i<xsize ; ++i)
      for (j=0 ; j<ysize ; ++j)
        stuff[j*xsize+i] = dp[j][i];
    return *this;
  }

void Matrix::resize(int b, int a)
{
  double *xx;
  int i,j,len=a*b;

  xx = new double[len];

  // always clips or pads
  for (i=0 ; i<(ysize > b ? ysize : b) ; ++i)
    for (j=0 ; j<(xsize > a ? xsize : a) ; ++j)
      {
	if ((i<b) && (j<a))
	  {
	    if ((i<ysize) && (j<xsize))
	      xx[i*a+j] = stuff[i*xsize+j];
	    else
	      xx[i*a+j] = 0;
	  }
      }

  delete[] stuff;
  stuff = xx;
  xsize=a;   ysize=b;
  maxels = len;
}

Matrix Matrix::operator+(Matrix m)
  {
    Matrix temp(ysize, xsize);
    // Add this to Matrix m
    int i,index,len=xsize*ysize;
    for (i=0 ; i<len ; ++i)
      temp.stuff[i] = stuff[i] + m.stuff[i];
    return temp;
  }

Matrix Matrix::operator-(Matrix m)
  {
    Matrix temp(ysize, xsize);
    // Subtract m from this
    int i,index,len=xsize*ysize;
    for (i=0 ; i<len ; ++i)
      temp.stuff[i] = stuff[i] - m.stuff[i];
    return temp;
  }

Matrix Matrix::operator*(Matrix other)
  {
    Matrix temp(ysize, other.xsize);
    int r1=nrows(),c1=ncols(),r2=other.nrows(),c2=other.ncols();

    // use BLAS level 3 dgemm routine
    // want m' * this'
    // note: for dgemm, we already store transpose

#ifdef USE_LAPACK
    dgemm_("N","N",&c2,&r1,&r2,&const_one,&other[0][0],&c2,
	   stuff,&c1,
	   &const_zero,
	   &temp[0][0],&c2);
#else
    int i,j,k;
    double tx;
    for (i=0 ; i<ysize ; ++i)
      for (j=0 ; j<other.xsize ; ++j)
	{
	  for (k=0,tx=0 ; k<xsize ; ++k)
	    tx += stuff[i*xsize+k]*other.stuff[k*other.xsize+j];
	  temp[i][j] = tx;
	}
#endif
    return temp;
  }

Vector Matrix::operator*(Vector v)
{
  // now use BLAS level 2 routine
  Vector temp(ysize);

#ifdef USE_LAPACK
  int m=ncols(),n=nrows();
  dgemv_("T",&m,&n,&const_one,stuff,&m,v.stuff,&const_int_one,&const_zero,temp.stuff,
	 &const_int_one);
#else
  int i,j,k,index;
  for (j=0 ; j<temp.ysize ; ++j)
    {
      for (temp[j]=0.0, k=0 ; k<xsize ; ++k)
	temp[j] += stuff[j*xsize + k] * v[k];
    }
#endif
  return temp;
}

Matrix Matrix::operator*(double k)
{
  Matrix temp(ysize, xsize);
  int i,j;
  for (i=0 ; i<temp.xsize ; ++i)
    for (j=0 ; j<temp.ysize ; ++j)
      temp[j][i] = k*(*this)[j][i];
  return(temp);
}

Matrix& Matrix::operator*=(double k)
{
  int i;
  for (i=0 ; i<xsize*ysize; ++i)
    stuff[i] *= k;
  return *this;
}

Matrix Matrix::operator/(double k)
{
  Matrix temp(ysize, xsize);
  int i,j;
  for (i=0 ; i<temp.xsize ; ++i)
    for (j=0 ; j<temp.ysize ; ++j)
      temp[j][i] = (*this)[j][i]/k;
  return(temp);
}

Matrix& Matrix::operator/=(double k)
{
  int i;
  for (i=0 ; i<xsize*ysize; ++i)
    stuff[i] /= k;
  return *this;
}

void Matrix::write_file(const char *c)
{
  std::ofstream outfile(c);
  int i,j;

  for (i=0 ; i<nrows() ; ++i)
    {
      for (j=0 ; j<ncols() ; ++j)
	outfile << (*this)[i][j] << " ";
      outfile << endl;
    }
  outfile.close();
}

int Matrix::read_file(const char *c)  // c is filename (or NULL)
  {
    char fnm[40], *ln;
    int i,j,numcols,dummy_lines,numrows,temprowcount;
    double tx;
    Matrix latestrow;
    int rowunit = 512;

    ln = new char[32768];      // buffer for individual line

    if (c==NULL)
      { cout << "Enter file name: ";   cin >> fnm;  }
    else strcpy(fnm, c);

    // open, etc.
    std::ifstream infile(fnm);

    // read automatically and adjust size
    dummy_lines = -1;
    do {
      ++dummy_lines;
      infile.getline(ln, 32768);
    } while ((ln[0]=='%') || (strlen(ln)==0));

    // trim off white space at end
    i = strlen(ln);
    for (j=i-1 ; j>=0 ; --j)
      if ((ln[j]==' ') || (ln[j]=='\t'))
	ln[j] = 0;
      else
	j=-1;
    std::istrstream ss(ln);
    streambuf *sb = ss.rdbuf();
    numcols = 0;
    while (sb->in_avail()>0)
      {
	ss >> tx;
	++numcols;
      }
    resize(1, numcols);
    latestrow.resize(rowunit, numcols);
    temprowcount = 0;

    // now go back to beginning and start again
    infile.seekg(0);
    for (i=0 ; i<dummy_lines ; ++i) // clear out comments and blank lines at beginning
      infile.getline(ln, 32768);

    for (i=0 ; i<numcols ; ++i)
      infile >> (*this)[0][i];
    while (!infile.eof())
      {
	for (i=0 ; i<numcols ; ++i)
	  infile >> latestrow[temprowcount][i];
	++temprowcount;
	if (temprowcount==rowunit)
	  {
	    temprowcount = 0;
	    // append it
	    append_bottom(latestrow);
	  }
      }
    // finish off latestrow buffer
    if (temprowcount>0)
      {
	latestrow = latestrow.submatrix(0,0,temprowcount,numcols);
	append_bottom(latestrow);
      }
    delete[] ln;

    // now chop off bottom row
    --ysize;

    infile.close();

    return(nrows()*ncols());
  }

complex *Matrix::eigenvalues()
// we do this in two parts
// (1) get characteristic equation
// (2) find roots
  {
    int i,j;
    if (xsize!=ysize)
      { cout << "Eigenvalue dimension mismatch."; return NULL; }
    Matrix temp(ysize, xsize), mp(ysize+1, 1), mcoeff(ysize+1, 1);
    // first evaluate char. poly. at 0,1,2,...,ysize
    temp = *this;
    for (i=0 ; i<=ysize ; ++i)
      {
        temp.singular(&(mp[i][0]));
        for (j=0 ; j<ysize ; ++j)
          --temp[j][j];
      }
    // now create Matrix
    Matrix temp2(ysize+1, xsize+1);
    for (i=0 ; i<=ysize ; ++i)
      {
        temp2[i][0]=1;
        for (j=1 ; j<=xsize ; ++j)
           temp2[i][j] = temp2[i][j-1]*i;
      }
    // and solve for coefficients of ch. poly.
    mcoeff = temp2.inverse()*mp;
    // now create array of complex and store roots of poly
    double *coeffs = new double[xsize+1];
    for (i=0 ; i<=xsize ; ++i)
      coeffs[i] = mcoeff[i][0];
    complex *roots = rootsof(xsize, coeffs);
    delete[] coeffs;
    // simple bubble sort of eigenvalues by magnitude
    // so they end up going from largest to smallest
    complex tc;
    for (i=0 ; i<xsize-1 ; ++i)
      for (j=0 ; j<xsize-1 ; ++j)
        {
          if (abs(roots[j])<abs(roots[j+1]))
            {
              tc=roots[j];
              roots[j] = roots[j+1];
              roots[j+1]=tc;
            }
        }
    return(roots);
  }

Matrix Matrix::symm_eigenvectors(double *evals)
{
  // this crude routine finds eigenvectors for corresponding eigenvalues by
  // computing (A - \lambda I) and solving (A - \lambda I)v = 0,
  // after replacing the last row in (A-\lambda I) by a row of ones,
  // and the last entry in the zero vector by a 1.0 (for scaling)

  // eigenvectors are returned in columns

  int i,j,dim=ncols();
  Matrix tm, retval(dim,dim);
  Vector tv(dim);
  for (i=0 ; i<dim ; ++i)
    {
      tm = *this;
      for (j=0 ; j<dim ; ++j)
	tm[j][j] -= evals[i];
      for (j=0 ; j<dim ; ++j)
	tm[dim-1][j] = 1.0;
      
      tv.zeroes();
      tv[dim-1] = 1.0;

      tm.solve_system(&tv);
      tv = tv * (1.0/tv.norm());

      for (j=0 ; j<dim ; ++j)
	retval[j][i] = tv[j];
    }

  return retval;
}

double Matrix::norm()
  {
    Matrix t=transpose();
    Matrix tx=t*(*this);
    complex *evals = tx.eigenvalues();
    double retval=sqrt(abs(evals[0]));
    delete[] evals;
    return retval;
  }

Vector Matrix::extract_row(int rownum)
{
  Vector tv(ncols());
  double *dp = (*this)[rownum];
  int i;
  for (i=0 ; i<ncols() ; ++i)
    tv[i] = dp[i];
  return tv;
}

Vector Matrix::extract_column(int colnum)
{
  Vector tv(nrows());
  double *dp = &stuff[colnum];
  int i,gap=ncols();
  for (i=0 ; i<nrows() ; ++i)
    tv[i] = dp[i*gap];
  return tv;
}

void Matrix::append_bottom(Matrix& appending, int trans)
{
  int i,j;
  Matrix bottom(appending);
  if (trans)
    bottom = bottom.transpose();

  int nr=nrows(),nc=ncols(),
    nrb=bottom.nrows(), ncb=bottom.ncols(), newlen;
  const int extrapad=32;
  double *newstuff;

  if (nc!=ncb)
    {
      cout << "Matrix::append_bottom() error, mismatching # of columns" << endl;
      return;
    }

  newlen = nc*(nr+nrb);
  if (newlen>maxels)
    {
      newstuff = new double[newlen+nc*extrapad];
      
      for (i=0 ; i<nc*nr ; ++i)
	newstuff[i] = stuff[i];
      for (i=nc*nr ; i<newlen ; ++i)
	newstuff[i] = bottom.stuff[i-nc*nr];

      delete[] stuff;
      stuff = newstuff;
      maxels = newlen + nc*extrapad;
    }
  else
    {
      for (i=nc*nr ; i<newlen ; ++i)
	stuff[i] = bottom.stuff[i-nc*nr];
    }
  ysize += nrb;
}

Matrix Matrix::submatrix(int r1, int c1, int r2p1, int c2p1)
{
  Matrix tm(r2p1-r1, c2p1-c1);
  int i,j;
  double *tx;
  for (i=r1 ; i<r2p1 ; ++i)
    {
      // copy a row
      tx = &(stuff[i*ncols()+c1]);
      for (j=0 ; j<c2p1-c1 ; ++j)
	tm[i-r1][j] = tx[j];
    }
  return tm;
}

Vector Matrix::mean()
{
  Vector tv(ncols());
  int i,j;

  tv.zeroes();
  for (i=0 ; i<nrows() ; ++i)
    for (j=0 ; j<ncols() ; ++j)
      tv[j] += (*this)[i][j];
  tv /= nrows();
  return tv;
}

double Matrix::kernel_function(double x, int n, double sd)
{
  double bw = sd*pow(n,-0.4);
  double retval = (1/sqrt(2*M_PI))*exp(-x*x/(2*bw*bw))/bw;
  return retval;
}

double Matrix::kernel_cdf(double x, int n, double sd)
{
  double bw = sd*pow(n,-0.4);
  double retval = norm_cdf(x, 0, bw*bw);
  return retval;
}

Vector Matrix::quantiles(double alpha, bool kernel_smooth)
{
  int i,j,k,n=nrows();
  int i0 = (int)floor(n*alpha - 0.5), i1 = i0+1;
  double frac = n*alpha - 0.5 - i0, x0, x1, tx, tF, tf, tsd, deltax;
  Vector tv(n), retval(ncols());

  if (i1>=n)    i1=n-1;
  if (i0>=n)    i0=n-1;
  if (i1<0)     i1=0;
  if (i0<0)     i0=0;

  for (j=0 ; j<ncols() ; ++j)
    {
      // extract column
      for (i=0 ; i<n ; ++i)
	tv[i] = (*this)[i][j];
      // and sort it
      sort_objects(n, &tv[0], NULL);

      // then compute quantile
      x0 = tv[i0];   x1 = tv[i1];
      tx = x0*(1-frac) + x1*frac;
      retval[j] = tx;

      if (kernel_smooth)
	{
	  tsd = sqrt(tv.var());

	  // do some Newton-Raphson iteration (6 steps)
	  for (k=0 ; k<6 ; ++k)
	    {
	      // evaluate F and f at current point
	      for (tf=0,tF=0,i=0 ; i<n ; ++i)
		{
		  tf += kernel_function(tx-tv[i], n, tsd);
		  tF += kernel_cdf(tx-tv[i], n, tsd);
		}
	      tf /= n;   tF /= n;
	      
	      // now take a Newton-Raphson step to get to the point where F(tx)=alpha
	      deltax = (alpha-tF)/tf;
	      tx += deltax;
	    }

	  retval[j] = tx;
	}
    }

  return retval;
}

Matrix Matrix::covariance()
{
  Matrix tm(ncols(),ncols());
  Vector tv(ncols());
  int i,j;

  tm.zeroes();
  for (i=0 ; i<nrows() ; ++i)
    {
      tv = extract_row(i);
      tm += tv.outer(tv);
    } 
  tm = tm * 1.0/nrows();  // now we have E(X X^T)
  tv = mean();
  tm -= tv.outer(tv);
  return tm;
}


//-----------------------------------------------
//    Code for class Vector next
//-----------------------------------------------

Vector::Vector()
  : Matrix(1, 1)
{
}

Vector::Vector(int length)
  : Matrix(length, 1)
{
}


Vector Vector::operator-(Vector othervec)
{
  Vector temp(nrows());
  for (int i=0 ; i<nrows() ; ++i)
    temp[i] = (*this)[i]-othervec[i];
  return temp;
}

Vector Vector::operator*(const double& other)
{
  Vector temp(nrows());
  for (int i=0; i<nrows() ; ++i)
    temp[i] = (*this)[i]*other;
  return temp;
}

Vector& Vector::operator=(const Matrix& otherMatrix)
  {
    Matrix *mp = (Matrix *) &otherMatrix;
    if (mp->ncols()!=1)
      cout << "Bad Matrix to vector conversion error!" << endl;

    Matrix::operator=(otherMatrix); 
    return *this;
  }

double Vector::dot(const Vector& othervector)
{
  Vector *mp = (Vector *) &othervector;
  double sum=0.0;
  for (int i=0 ; i<mp->nrows() ; ++i)
    sum += stuff[i]*(*mp)[i];
  return sum;
}

double Vector::norm()
{
  double sum=0.0;
  for (int i=0 ; i<nrows() ; ++i)
    sum += stuff[i]*stuff[i];
  return sqrt(sum);
}

double Vector::sum()
{
  double sum=0.0;
  for (int i=0 ; i<nrows() ; ++i)
    sum += stuff[i];
  return sum;
}

double Vector::mean()
{
  return sum()/nrows();
}

double Vector::var()
{
  int i;
  double tx,ty;
  Vector squares(*this);
  for (i=0 ; i<nrows() ; ++i)
    squares[i] = squares[i]*squares[i];
  tx = squares.mean();
  ty = mean();
  return (tx-ty*ty);
}

Matrix Vector::outer(const Vector& other)
{
  Vector *mp = (Vector *) &other;
  Matrix temp(nrows(), (*mp).nrows());
  int i,j;
  for (i=0 ; i<nrows() ; ++i)
    for (j=0 ; j<mp->nrows() ; ++j)
      temp[i][j] = stuff[i]*(*mp)[j];
  return temp;
}

ostream& operator<<(ostream& os, Vector& m)
  // stream it in Matlab format
{
  int i,j,vecmode=0;
  if (m.iomode==MATLAB_MODE)
    {
      os << "[";
      for (i=0 ; i<m.nrows() ; ++i)
	{
	  os << m[i];
	  if (i<m.nrows()-1)
	    os << ",";
	}
      os << "]'";
    }
  if (m.iomode==DISP_MODE || m.iomode==ROW_MODE)
    {
      for (i=0 ; i<m.nrows() ; ++i)
	os << m[i] << " ";
    }
  return(os);
}

istream& operator>>(istream& is, Vector& m)
{
  int i,j,k,done,maxels,rows,cols;
  double *localstuff,*newstuff;
  char c;
  
  localstuff = new double[100];
  maxels = 100;
  is >> ws;
  is.get(c);  // remove "["
  done=0;

  i = 0;
  c = is.peek();
  if (c==']')
    {
      is.get(c);
      done = true;
    }
  while (!done)
    {
      is >> localstuff[i++];
      if (i==maxels)
	{
	  // we need to expand memory
	  maxels *= 2;
	  newstuff = new double[maxels];
	  for (k=0 ; k<maxels/2 ; ++k)
	    newstuff[k] = localstuff[k];
	  delete[] localstuff;
	  localstuff = newstuff;
	}
      is >> ws;
      is.get(c);   
      if ((c!=',') && (c!=']'))
	is.unget();
      done=(c==']');
    }
  c = is.peek();
  if (c=='\'')
    is.get(c);  // get rid of transpose char. ' at end
  m.set_data(localstuff,i,1,maxels);
  return is;
}  

double Vector::angle()
{
  double ang,n;
  n = norm();
  if (n<1e-8)
    n=1e-8;  // avoid 0 problems
  ang = asin(stuff[1]/n);
  if (stuff[0]<0)
    ang = M_PI - ang;
  while (ang<0)
    ang += 2*M_PI;
  while (ang>2*M_PI)
    ang -= 2*M_PI;
  return ang;
}

void Vector::append(double newone)
{
  int nr = ysize;
  resize(nr+1);
  stuff[nr] = newone;
}

Vector Vector::subvector(int r1, int r2)
{
  int i;
  Vector tv(r2-r1);
  for (i=r1 ; i<r2 ; ++i)
    tv[i-r1] = stuff[i];
  return tv;
}

}

void Vector::sort()
{
  int l,j,ir,i,n=ysize;
  double rra,*dp = stuff;

  if (n<2) return;	       // we don't need to sort if n=1

  l = (n >> 1)+1;
  ir = n;

  for (;;) {
    if (l>1)		// i.e. we are still in hiring phase
      rra = dp[--l-1];
    else {
      rra=dp[ir-1];
      dp[ir-1]=dp[0];
      if (--ir == 1) {	// we are finished
	dp[0]=rra;
	return;
      }
    }
    i=l;		j=l<<1;
    while (j <= ir)
      {
	if ((j < ir) && (dp[j-1]<dp[j])) ++j;	// j is better underling
	if (rra<dp[j-1]) {		// demote rra
	  dp[i-1]=dp[j-1];
	  j += (i=j);
	}
	else j=ir+1;			// finished sifting
      }
    dp[i-1]=rra;
  }
}
