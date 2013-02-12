/*
 * -------------------------------------------------------------------
 *
 * Copyright 2004 Anthony Brockwell
 *
 * This program is free software; you can redistribute it and/or
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
 * This class handles truncated power series splines.
 * Anthony Brockwell, June 11th, 2002.
 */


#include <iostream>
#include <fstream>
#include <iomanip>
#include <strstream>
#include "matrix.h"
#include "cubespline.h"
#include "glm.h"

using namespace std;

cubespline::cubespline(int nknots)
  : numknots(nknots)
{
  if (numknots>0) {
    knotpos.resize(numknots);
  }
}

cubespline::~cubespline()
{
}

double cubespline::clip(double x)
{
  return x > 0 ? x : 0;
}

double cubespline::cube(double x)
{
  return x*x*x;
}

void cubespline::setknots(Vector& knots)
{
  int i, n = knots.nrows();
  if (n!=numknots)
    {
      cout << "Wrong number of knots in cubespline::setknots!" << endl;
      return;
    }
  knotpos = knots;
}

Vector cubespline::knotsforgaps(double min, double max, Vector dirch)
{
  int i,n=dirch.nrows()-1;
  double tx;
  Vector retval(n);

  // 1. make sure that dirch sums to 1.0
  for (i=0 ; i<dirch.nrows() ; ++i)
    dirch[i] = fabs(dirch[i]); // only non-neg. entries are allowed
  tx = dirch.sum();
  dirch = dirch * (1.0/tx);  // normalize it

  // 2. translate into knot positions
  retval[0] = dirch[0]*(max-min)+min;
  for (i=1 ; i<n ; ++i)
    retval[i] = retval[i-1] + dirch[i]*(max-min);

  return retval;
}

Vector cubespline::gapsforknots(double min, double max, Vector knots)
{
  int i,n=knots.nrows()+1;
  double tx=max-min,upper,lower;
  Vector retval(n);

  for (i=0 ; i<n ; ++i)
    {
      upper = i<n-1 ? knots[i] : max;
      lower = i>0 ? knots[i-1] : min;
      retval[i] = (upper-lower)/tx;
    }
  return retval;
}

double cubespline::cost(Vector& knotgaps, void *vp)
  // This function returns sum squared residuals when
  // regression fit for given knot positions is carried out.
  // Hence it can be used to get MLEs for Gaussian free-knot splines.
{
  cubespline *cp = (cubespline *) vp;
  Vector knots;
  knots = cp->knotsforgaps(cp->minx,cp->maxx,knotgaps);
  cp->setknots(knots);
  cp->lsfit(*(cp->tempdata)); // fit it
  return cp->sig2;
}

double cubespline::cost2(Vector& gapsandparms, void *vp)
  // This function returns the negative Poisson log-likelihood.
  // It requires both gaps (numknots+1) and params (3+numknots),
  // and can be used to find MLEs for Poisson free-knot splines.
{
  cubespline *cp = (cubespline *) vp;
  Vector gaps(cp->numknots+1), knots;
  gaps = gapsandparms.submatrix(0,0,cp->numknots+1,1);
  knots = cp->knotsforgaps(cp->minx,cp->maxx,gaps);
  cp->setknots(knots);

  // extract coefficients as well
  cp->beta = gapsandparms.submatrix(cp->numknots+1,0,2*cp->numknots+4,1);

  // now we can compute and return negative Poisson log-likelihood
  return -cp->loglikelihood(*(cp->tempdata),1);
}

double cubespline::cost3(double *xpos, void *vp)
  // This function returns -fittedvalue(x) and is used
  // to find argmax.
{
  cubespline *cp = (cubespline *) vp;
  double tx = 1.0/(1.0+exp(-xpos[0]));  // restrict range to [0,1]
  return -cp->getfitfor(tx);
}

double cubespline::cost4(Vector& betas, void *vp)
  // This function returns -fittedvalue(x) and is used
  // to find argmax.
{
  cubespline *cp = (cubespline *) vp;
  cp->beta = betas;
  return -cp->loglikelihood(*(cp->tempdata),1);
}

void cubespline::fitPn(Matrix& data, int fileit)
{
  int i,j,n = data.nrows();
  Matrix XtX;
  Vector Y(n), betahat(3+numknots);
  double tx;
  PoissonGLM pglm;

  X = getdesignmatrix(data);
  Y = data.submatrix(0,1,n,2);

  /*
  if (fileit)
    {
      X.writefile("xout");
      Y.writefile("yout");
    }
  */

  betahat = pglm.Estimate(Y,X);

  // copy from vector into class members
  beta = betahat;
  betacov = *(pglm.GetBetahatCovInv());

  // finally, sort out variance
  fitted = X*betahat;
  for (i=0 ; i<n ; ++i)
    fitted[i] = exp(fitted[i]);

  sig2 = 0;  // irrelevant for Poisson spline
}

Matrix cubespline::getdesignmatrix(Matrix& data)
{
  int i,j,n=data.nrows();
  double tx;
  Matrix X(n, 3+numknots);
  // fill in design matrix
  for (i=0 ; i<n ; ++i)
    {
      tx = data[i][0];
      X[i][0] = 1.0;
      X[i][1] = tx;
      X[i][2] = tx*tx;
      for (j=0 ; j<numknots ; ++j)
	X[i][3+j] = clip(cube(tx-knotpos[j]));
    }
  return X;
}

double cubespline::getfitfor(double x)
{
  int i;
  double tx;
  tx = beta[0] + beta[1]*x
    + beta[2]*x*x;
  for (i=0 ; i<numknots ; ++i)
    tx += beta[3+i] * clip(cube(x-knotpos[i]));
  return tx;
}


void cubespline::lsfit(Matrix& data)
  // This function fits the spline by least squares,
  // and assumes that knotpos[0..numknots-1] have been
  // filled in.
{
  int i,j,n = data.nrows();
  Matrix XtX;
  Vector Y(n), betahat(3+numknots);
  double tx,ss;

  X = getdesignmatrix(data);
  X.write_file("x.out");
  Y = data.submatrix(0,1,n,2);
  Y.write_file("y.out");
  XtX = X.transpose()*X;
  betahat = X.transpose()*Y;
  XtX.symm_solve_system(&betahat);

  // copy from vector into class members
  beta = betahat;

  // finally, sort out variance
  fitted = X*betahat;
  Vector resids;
  resids = fitted - Y;  // get residuals
  ss = resids.dot(resids);
  sig2 = ss/n;

  // fill in betacov!!!
  // empty for now!!!???
}

double logfact(double x)
{
  int i,ix = floor(x+0.5);
  double tx=0;
  for (i=2 ; i<=ix ; ++i)
    tx += log(((double)i));
  return tx;
}

double cubespline::loglikelihood(Matrix& data, int mode, double disp)
{
  int i,n=data.nrows();
  Matrix X;
  double ll,lambda;

  X = getdesignmatrix(data);
  fitted = X*beta;  // fitted value

  if (mode==0)
    {
      // if mode is 0 we return log Gaussian likelihood
      Vector resids, Y;
      Y = data.submatrix(0,1,n,2);
      resids = fitted - Y;  // get residuals
      sig2 = resids.dot(resids)*disp*disp/n;
      ll = 0;
      for (i=0 ; i<n ; ++i)
	ll += -0.5*log(2*M_PI*sig2)-0.5*
	  (resids[i])*(resids[i])/sig2;
    }
  else if (mode==1)
    {
      // if mode is 1 we return Poisson likelihood
      ll = 0;
      for (i=0 ; i<n ; ++i)
	{
	  fitted[i] = exp(fitted[i]);
	  lambda = fitted[i];
	  ll += -lambda + data[i][1]*log(lambda) - logfact(data[i][1]);
	}
    }
  else return 0;
  return ll;
}

void cubespline::display(FILE *gnupipe, double minx, double maxx,
			 double miny, double maxy, Matrix *data)
{
  int i;
  if (gnupipe==NULL)
    {
      cout << setprecision(10) << "Betas: " << beta << endl;
      cout << "Knotpos: " << knotpos << endl;
      cout << "Sig^2 = " << sig2 << endl;
    }
  else
    {
      // we should plot it in gnuplot
      char cmd[2048];
      ostrstream os(cmd, 2048);

      ofstream knotsout("knots.out");
      for (i=0 ; i<numknots ; ++i)
	{
	  knotsout << knotpos[i] << " " << miny << endl;
	  knotsout << knotpos[i] << " " << maxy << endl << endl;
	}
      knotsout.close();

      // basic definitions
      os << "g(x)=x>0?x:0" << endl;
      os << "c3(x,k) = g((x-k)**3)" << endl;
      os << ends;
      fputs(cmd, gnupipe);
      fflush(gnupipe);

      // now the actual spline
      os.seekp(0);
      os << setprecision(14);
      os << "f(x)=" << beta[0] << "+" << beta[1] << "*x+"
	 << beta[2] << "*x*x + "
	 << beta[3] << "*c3(x," << knotpos[0] << ")";
      for (i=1 ; i<numknots ; ++i)
	os << "+" << beta[3+i] << "*c3(x," << knotpos[i] << ")";
  
      os << endl;
      os << "plot [" << minx << ":" << maxx << "] [";
      os << miny << ":" << maxy << "] ";
      if (data!=NULL)
	{
	  ofstream outfile("temp.out");
	  int ndata = data->nrows();
	  for (i=0 ; i<ndata ; ++i)
	    outfile << (*data)[i][0] << " " << (*data)[i][1] << endl;
	  outfile.close();
	  os << "\"temp.out\",";
	}
      os << "exp(f(x)), \"knots.out\" with lines" << endl << ends;
      fputs(cmd, gnupipe);
      fflush(gnupipe);
    }
}

Vector cubespline::getbeta(Matrix& cov)
{
  cov = betacov;
  return beta;
}
