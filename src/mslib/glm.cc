/*
 * -------------------------------------------------------------------
 *
 * Copyright 2004 Anthony Brockwell
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
 *
 * This module contains code for GLM fitting.
 * It uses iteratively reweighted least squares (see McCullagh & Nelder).
 * Nov. 2002, Anthony Brockwell.
 *
 * Base class need functions linkfn, linkderiv,
 *   variancefn and dispersion filled in. 
 *
 */

#include <math.h>
#include "matrix.h"
#include "abrand.h"
#include "glm.h"

namespace mslib {

GenLinearModel::GenLinearModel()
{
}

Vector GenLinearModel::Estimate(Vector& y, Matrix& X, Vector *initial,
				 int maxits)
{
  int i, j, iter, done=0;
  double epsilon[2]={-1,-1},newepsilon,det;
  Vector z,tv;
  Matrix WX,tm;
  int p = X.ncols(), n = X.nrows();
  Vector eta0(n), eta1(n);

  wts.resize(n);

  if (y.nrows()!=n)
    {
      cout << "GLM error: y does not have the right number of rows." << endl;
      return tv;
    }
  WX.resize(n,p);
  betahat.resize(p);
  tv.resize(p);

  // find some kind of initializer = LS fit for inv. transformed y[i]
  if (initial==NULL)
    {
      for (i=0 ; i<n ; ++i)
	eta0[i] = invlinkfn(y[i]);
      tm = X.transpose()*X;
      betahat = X.transpose()*eta0;
      tm.robust_symm_solve_system(&betahat);
    }
  else
    betahat = *initial;

  for (iter=0 ; (iter<maxits) && (!done) ; ++iter)
    {
      // compute z,W and update betahat
      //      cout << "[" << iter << "]" <<  betahat << 
      //	" : " << epsilon[1] << endl;

      eta0 = X*betahat;
      tv = y;
      for (i=0 ; i<n ; ++i)
	{
	  tv[i] -= linkfn(eta0[i]);
	  tv[i] /= linkderiv(eta0[i]);
	}
      z = eta0 + tv;
      for (i=0 ; i<n ; ++i)
	{
	  wts[i] = linkderiv(eta0[i]);
	  wts[i] = wts[i]*wts[i]/(dispersion(i)*variancefn(eta0[i]));
	}

      for (i=0 ; i<n ; ++i)
	for (j=0 ; j<p ; ++j)
	  WX[i][j] = X[i][j]*wts[i];

      for (i=0 ; i<n ; ++i)
	z[i] *= wts[i];
      z = X.transpose()*z;
      tm = X.transpose()*WX;
      tm.robust_symm_solve_system(&z);

      eta1 = X*z;
      newepsilon = (eta1-eta0).norm();
      epsilon[0]=epsilon[1];
      epsilon[1]=newepsilon;
      if (epsilon[0]!=-1)
	if (epsilon[1]>epsilon[0])
	  done=1;

      betahat = z;
    }

  // compute approximate cov. matrix for betahat
  eta0 = X*betahat;
  for (i=0 ; i<n ; ++i)
    {
      wts[i] = linkderiv(eta0[i]);
      wts[i] = wts[i]*wts[i]/(dispersion(i)*variancefn(eta0[i]));
    }
  for (i=0 ; i<n ; ++i)
    for (j=0 ; j<p ; ++j)
      WX[i][j] = X[i][j]*wts[i];
  betacov = (X.transpose()*WX);

  return betahat;
}

Vector GenLinearModel::Fitted(Matrix& X)
  // returns E[Y|X] based on last fit
{
  int i;
  Vector y = X*betahat;
  for (i=0 ; i<y.nrows() ; ++i)
    y[i] = linkfn(y[i]);
  return y;
}

/* 
** Next we have class for logistic regression
** or more generally, for binomial response
** with the logistic link function, scaled to
** be between 0 and 1.
*/

BinomialGLM::BinomialGLM(int m)
  : binm(m), GenLinearModel()
{
}

double BinomialGLM::dispersion(int i)
{
  return 1.0/binm;
}

double BinomialGLM::linkfn(double x)
  // $\mu = linkfn(X \beta)$
  // inverse of link fn as defined by Firth
{
  return exp(x)/(1.0+exp(x));  // default is identity
}

double BinomialGLM::invlinkfn(double x)
{
  double eps=0.001;
  if (x<eps)
    x=eps;
  if (x>1-eps)
    x=1-eps;
  return log(x/(1.0-x));
}

double BinomialGLM::linkderiv(double x)
{
  return linkfn(x)/(1.0+exp(x));
}

double BinomialGLM::variancefn(double etai)
  // returns variance of y_i for eta_i, divided by dispersion
{
  double tx = linkfn(etai);
  return tx*(1.0-tx);
}

PoissonGLM::PoissonGLM()
  : GenLinearModel()
{
}

double PoissonGLM::linkfn(double x)
{
  return exp(x);
}

double PoissonGLM::invlinkfn(double x)
{
  double eps=0.1;
  if (x<eps)
    x=eps;
  return log(x);
}

double PoissonGLM::linkderiv(double x)
{
  return exp(x);
}

double PoissonGLM::variancefn(double etai)
{
  double tx = linkfn(etai);
  return tx;
}

double PoissonGLM::dispersion(int i)
{
  return 1.0;
}

double PoissonGLM::scaleddeviance(Vector& y, Matrix& X)
{
  // assumes betahat has already been filled in
  int i;
  double part1,part2;
  Vector fitted = X*betahat;
  for (i=0 ; i<fitted.nrows() ; ++i)
    fitted[i] = exp(fitted[i]);

  // 1. work out log-likelihood of exactly fitted values
  for (i=0,part1=0 ; i<y.nrows() ; ++i)
    part1 += logpois_density(y[i], y[i]);
  // 2. work out other bit
  for (i=0,part2=0 ; i<y.nrows() ; ++i)
    part2 += logpois_density(y[i], fitted[i]);
  return (2*(part1-part2));
}

}
