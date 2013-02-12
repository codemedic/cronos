/* -------------------------------------------------------------------
 *
 * Copyright 2004,2005,2006 Anthony Brockwell
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
 * This module contains the source code defining
 * the GARCHModel class, derived from TimeSeriesModel.
 */

#include <iomanip.h>
#include <strstream>
#include <complex>
#include <string>
#include <functional>
#include <gsl/gsl_specfunc.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include "matrix.h"
#include "abrand.h"
#include "poly.h"
#include "garch.h"

using namespace std;

const int alloc_unit = 256, KFburnin=12;
extern double accuracy;  // used in FindMinimum
const double zero=1e-10;

GARCHModel::GARCHModel()
  : TimeSeriesModel()
{
  p=q=0;
  as.resize(1);
  as[0] = 1.0;
  bs.resize(1);
  mean = 0.0;
  loglikelihood_is_valid = false;
}

GARCHModel::~GARCHModel()
{
}

void GARCHModel::SetParameters(int p0, int q0, Vector& as0, Vector& bs0, double m)
{
  int i;
  p=p0;   q=q0;
  as.resize(p+1);
  bs.resize(q);
  for (i=0 ; i<=p ; ++i)
    as[i] = as0[i];
  for (i=0 ; i<q ; ++i)
    bs[i] = bs0[i];
  mean=m;
}

void GARCHModel::GetParameters(int& p1, int& q1, Vector& as1, Vector& bs1, double& mu)
{
  p1=p;   q1=q;
  as1=as;   bs1=bs;
  mu=mean;
}

double GARCHModel::LogPrior(bool& isok)
{
  int i;
  double tx = 0.0;
  isok = as[0]>0.0;
  for (i=1 ; i<=p ; ++i)
    {
      isok &= as[i] >= 0.0;
      tx += as[i];
    }
  for (i=0 ; i<q ; ++i)
    {
      isok &= bs[i] >= 0.0;
      tx += bs[i];
    }
  isok &= (tx<1.0);
  return 0.0;
}

Vector GARCHModel::ParameterBundle()
{
  int i;
  Vector pvec(p+q+2);
  pvec[0] = mean;
  for (i=0 ; i<=p ; ++i)
    pvec[i+1] = as[i];
  for (i=0 ; i<q ; ++i)
    pvec[i+p+2] = bs[i];
  return pvec;
}

void GARCHModel::UnbundleParameters(Vector &pvec)
{
  int i;
  mean = pvec[0];
  for (i=0 ; i<=p ; ++i)
    as[i] = pvec[i+1];
  for (i=0 ; i<q ; ++i)
    bs[i] = pvec[i+p+2];
}

void GARCHModel::StreamParameterName(ostringstream& os, int pnum)
{
  if (pnum==0)
    os << "\u03bc = mean"; 
  else if (pnum<=(p+1))
    os << "\u03b1(" << (pnum-1) << ")";
  else
    os << "\u03b2(" << (pnum-p-1) << ")";
}

double GARCHModel::GetVar()
{
  int i;
  double tx = as[0], beta=0.0;
  for (i=1 ; i<=p ; ++i)
    beta += as[i];
  for (i=0 ; i<q ; ++i)
    beta += bs[i];
  tx /= (1.0-beta);
  return tx;
}

void GARCHModel::BayesianFit(int p, int q, TimeSeries *ts,
			     void (*callback)(int,void *), void *cb_parms)
{
  int i,j,iteration;
  double tx,ty,delta,scale;
  bool accept,isok;
  Matrix allparms;

  // start with IIDN model
  double mu = ts->SampleMean();
  Vector as(p+1),bs(q+1),acf(1),pvec;

  ts->ComputeSampleACF(&acf,NULL,0);
  as.zeroes();
  bs.zeroes();
  as[0] = acf[0];
  SetParameters(p,q,as,bs,mu);
  pvec = ParameterBundle();

  // create copy of timeseries with missing values filled in
  TimeSeries filled;
  for (i=0 ; i<ts->GetN() ; ++i)
    {
      if (!ts->IsMissing(i))
	filled.Append((*ts)[i]);
      else
	filled.Append(norm_random(mu,as[0]));
    }

  // now iterate
  allparms = pvec.transpose();
  tx = ComputeLogLikelihood(&filled) + LogPrior(isok);
  for (iteration=0 ; iteration < iterationlimit ; ++iteration)
    {
      // register the iteration
      if (callback!=NULL)
	(*callback)(2,cb_parms);

      // do an update for each of the parameters
      // make sure they preserve variance!
      for (i=1 ; i<=p+q ; ++i)
	{
	  scale = GetVar();
	  delta = norm_random(0,0.0001);
	  pvec[i+1] += delta;
	  UnbundleParameters(pvec);
	  scale /= GetVar();
	  pvec[1] *= scale;        // scale by changing a[0]
	  UnbundleParameters(pvec);
	  LogPrior(accept);  // if it's invalid, don't accept
	  if (accept)
	    {
	      ty = ComputeLogLikelihood(&filled) + LogPrior(isok);
	      accept = (log(unif_random())<(ty-tx));
	    }
	  if (accept)
	    tx = ty;
	  else
	    {
	      pvec[i+1] -= delta;
	      pvec[1] /= scale;
	      UnbundleParameters(pvec);
	    }
	}

      // do a variance update (pvec[0] = sigma^2)
      delta = norm_random(0,0.0001);
      delta = exp(delta);
      pvec[1] *= delta;
      UnbundleParameters(pvec);
      ty = ComputeLogLikelihood(&filled) + LogPrior(isok);
      if (log(unif_random())<(ty-tx))
	// accept
	tx = ty;
      else
	{
	  pvec[1] /= delta;
	  UnbundleParameters(pvec);
	}

      // update the missing values in filled
      // ???

      // record results
      allparms.append_bottom(pvec, 1);
    }

  // now fitted model is posterior mean
  pvec.zeroes();
  for (i=100 ; i<iterationlimit ; ++i)
    for (j=0 ; j<pvec.size() ; ++j)
      pvec[j] += allparms[i][j];
  pvec = pvec*(1.0/(iterationlimit-100));


  UnbundleParameters(pvec);
}

int GARCHModel::FitModel(TimeSeries *ts, const int method, const int ilimit,
			 void (*itercallback)(int,void *), void *cb_parms,
			 ostringstream& msg, ostringstream& supplemental,
			 bool get_parameter_cov)
{
  int i,j,nobs = ts->GetN(),pstar=max(p,q);
  Vector localas(p+1),localbs(q);
  Vector acf(1),mask;
  double mu = ts->SampleMean();
  Matrix parmcovs;

  iterationlimit = ilimit;
  mask = GetMask();  // get current mask (pad if necessary with default values)

  if (method==METHOD_MLE)
    {
      if (ts->HasAnyMissing())
	{
	  msg << "Unable to fit GARCH model when time series has missing observations."
	      << ends;
	  return UNABLE;
	}
      
      parmcovs = MLEFit(ts, mask, itercallback, cb_parms, get_parameter_cov);
      StreamMLESummary(supplemental, parmcovs, get_parameter_cov);

      msg << "GARCH model fit in " << iterationcounter << " iterations." << ends;

      CheckModel(false);
      return SUCCESS;
    }
  if (method==METHOD_BAYESIAN)
    {
      BayesianFit(p, q, ts, itercallback, cb_parms);
      msg << "MCMC estimation carried out with " << ilimit
	  << " iterations." << ends;
      supplemental << "Analysis of MCMC output not yet implemented." << ends;
      return SUCCESS;
    }
}

void GARCHModel::SimulateInto(TimeSeries *ts, int n, bool ovw)
{
  int i,j,nn;
  double sum;

  // clean up time series
  if (!ovw)
    {
      nn = n;
      ts->ClearOut();
      for (i=0 ; i<nn ; ++i)
        ts->Append(0);
    }
  else
    nn = ts->GetN();

  // simulate
  Vector acf(1);
  ComputeACFPACF(acf,NULL,0);     // get variance
  Vector nu(nn);
  nu[0] = acf[0];
  (*ts)[0] = norm_random(mean, acf[0]);
  for (i=1; i<nn ; ++i)
    {
      nu[i] = OneStepVariance(&(*ts)[i-1], &nu[i-1], i-1, acf[0]);
      (*ts)[i] = norm_random(mean, nu[i]);
    }
}

Matrix GARCHModel::SimulatePredictive(TimeSeries &ts, int nsteps, int nsims)
{
  // note: handling of missing values is not optimal.
  //       they are just filled in with draws from one-step predictive d-ns
  //       (instead of the full-conditional distributions)

  int i,j,n0 = ts.GetN(),nn = n0+nsteps, sim_num;
  double sum;
  Matrix retval(nsims, nsteps);

  // simulate
  Vector acf(1);
  ComputeACFPACF(acf,NULL,0);     // get variance
  Vector nu(nn), new_ts(nn);

  for (sim_num=0 ; sim_num<nsims ; ++sim_num)
    {
      nu[0] = acf[0];
      if (ts.IsMissing(0))
	new_ts[0] = norm_random(mean, acf[0]);
      else
	new_ts[0] = ts[0];
      for (i=1; i<nn ; ++i)
	{
	  nu[i] = OneStepVariance(&new_ts[i-1], &nu[i-1], i-1, acf[0]);
	  if (ts.IsMissing(i) || i>=n0)
	    new_ts[i] = norm_random(mean, nu[i]);
	  else
	    new_ts[i] = ts[i];
	}
      for (i=0 ; i<nsteps ; ++i)
	retval[sim_num][i] = new_ts[i+n0];
    }

  // return tail subvector
  return retval;
}

void GARCHModel::RenderInto(ostringstream& os)
{
  int i;
  os << "\\b";  // make it bold
  os << "GARCH(" << p << "," << q << ") Model" << endl;
  os << "\\n";  // make it normal
  os << "X(t) = \\gs\\n(t) Z(t)" << endl;
  os << "\\gs\\n(t)\\s2\\n = ";
  os << setiosflags(ios::fixed) << setprecision(3) << as[0];
  for (i=1 ; i<=p ; ++i)
    os << " + " << setiosflags(ios::fixed) << setprecision(3) << as[i] << "X(t-" << i << ")\\s2\\n";
  for (i=1 ; i<=q ; ++i)
    os << " + " << setiosflags(ios::fixed) << setprecision(3) << bs[i-1] << "\\gs(t-" << i << ")\\s2\\n";
  os << endl << "Mean = " << setprecision(5) << mean << ",  " << endl;
  if (loglikelihood_is_valid)
    {
      os << endl << endl << "Log Likelihood = " << GetLogLikelihood() << endl;
      os << "AICC = " << GetAICC() << endl;
    }
  else
    {
      os << endl << endl << "Log Likelihood = <Not Computed>" << endl;
      os << "AICC = <Not Computed>" << endl;
    }
}

bool GARCHModel::RenderCtsVersionInto(ostringstream &output, ostringstream &error, double delta)
{
  error << "This feature is not yet supported for GARCH models." << ends;
  return false;
}

void GARCHModel::ComputeACFPACF(Vector& acf, Vector* pacf, bool normalize)
{
  int i;
  double tx;
  acf.zeroes();
  if (pacf!=NULL)
    {
      (*pacf)[0] = 1.0;
      pacf->zeroes();
    }
  if (normalize)
    acf[0] = 1.0;
  else
    {
      tx = 1.0 - as.sum() + as[0];
      if (bs.nrows()>0)
	tx -= bs.sum();
      tx = as[0] / tx;
      acf[0] = tx;
    }
}

void GARCHModel::Forecast(TimeSeries *tp, int nsteps, 
			  double *fret, double *fmse)
{
  int i, j, n = tp->GetN();
  Vector acf(1);
  double st;

  ComputeACFPACF(acf,NULL,0);  // get variance only

  // compute predictive variances
  Vector nu(n+nsteps);
  nu[0] = acf[0];
  for (i=1 ; i<n+nsteps ; ++i)
    nu[i] = OneStepVariance(&(*tp)[i-1], &nu[i-1], i-1, acf[0]);


  // fill in forecasts
  for (i=0 ; i<nsteps ; ++i)
    {
      if (fret!=NULL)
	fret[i] = mean;
      if (fmse!=NULL)
	fmse[i] = nu[i+n];
    }
}

void GARCHModel::ComputeStdResiduals(TimeSeries *tp, TimeSeries *resids)
{
  int i, j, n = tp->GetN();
  Vector acf(1);
  double st;

  ComputeACFPACF(acf,NULL,0);  // get variance only

  // compute predictive variances
  Vector nu(n);
  nu[0] = acf[0];
  for (i=1 ; i<n ; ++i)
    nu[i] = OneStepVariance(&(*tp)[i-1], &nu[i-1], i-1, acf[0]);

  // resids are one-step predictors divided by their predictive variances
  resids->ClearOut();
  for (i=0 ; i<n ; ++i)
    resids->Append(((*tp)[i]-mean)/sqrt(nu[i]));
}

double GARCHModel::OneStepVariance(double *xt, double *st, int t, double gamma0)
{
  // This function returns conditional variance of x[t+1], given x[t]=*xt and sigma^2[t]=*st,
  // ..., x[0]=*(xt-t) and sigma^2[0]=*(st-t).  We assume that there is no valid
  // data stored at *(xt-t-k) or *(st-t-k) for k>=1.  Note how this function
  // is overridden in the derived class EGARCHModel.

  double sum=0;
  int i,j,m = max(p,q);
  if (t+1>=m)
    {      
      sum=as[0];
      for (j=1 ; j<=p ; ++j)
	sum += as[j]*(xt[1-j]-mean)*(xt[1-j]-mean);
      for (j=0 ; j<q ; ++j)
	sum += bs[j]*st[-j];
      return sum;
    }
  else
    {
      sum=as[0];
      for (j=1 ; j<=p ; ++j)
	if (1-j >= -t)
	  sum += as[j]*(xt[1-j]-mean)*(xt[1-j]-mean);
	else
	  sum += as[j]*gamma0;
      for (j=0 ; j<q ; ++j)
	if (j >= -t)
	  sum += bs[j]*st[-j];
	else
	  sum += bs[j]*gamma0;
      return sum;
    }
}

double GARCHModel::ComputeLogLikelihood(TimeSeries *tp)
{
  // computes conditional log-likelihood
  int i,j,m = max(p,q), nn=tp->GetN();
  double tx = 0, sum;
  Vector st(nn),acf(1);

  if (nn==0)
    return 0.0;

  ComputeACFPACF(acf,NULL,0);
  for (i=0 ; i<m ; ++i)
    st[i] = acf[0];
  for (i=m ; i<nn ; ++i)
    {
      st[i] = OneStepVariance(&(*tp)[i-1], &st[i-1], i-1, acf[0]);
      tx += lognorm_density((*tp)[i], mean, st[i]);
    }

  tx *= ((double)nn)/(nn-m);   // correct for length of time series (so that AICC makes sense)

  loglikelihood = tx;
  AICC = -2*loglikelihood + 2*(p+1+q);
  loglikelihood_is_valid = true;

  return loglikelihood;
}

bool GARCHModel::CheckModel(bool flipit)
{
  bool retval = true;
  if (!flipit)
    {
      Vector tv = ValidityConstraint();
      retval = (tv[0] >= 0) && (tv[1] >= 0); 
    }
  else // actually fix it
    {
      int i;
      double sum=0.0;
      for (i=0 ; i<=p ; ++i)
	{
	  as[i] = fabs(as[i]);
	  if (i>0)
	    sum += as[i];
	}
      for (i=0 ; i<q ; ++i)
	{
	  bs[i] = fabs(bs[i]);
	  sum += bs[i];
	}
      if (sum>(1-1e-10))
	{
	  double scale = (1-1e-10) / sum;
	  for (i=1 ; i<=p ; ++i)
	    as[i] *= scale;
	  for (i=1 ; i<=q ; ++i)
	    bs[i-1] *= scale;
	}
      retval = true;
    }
  return retval;
}

Vector GARCHModel::ValidityConstraint()
{
  int i;
  Vector retval(2);

  // first constraint: all parms non-neg.
  // put in minimum
  retval[0] = as[0]-zero;
  for (i=1 ; i<=p ; ++i)
    if (as[i]<retval[0])
      retval[0] = as[i];
  for (i=0 ; i<q ; ++i)
    if (bs[i]<retval[0])
      retval[0] = bs[i];

  // 2nd: 1-\sum
  retval[1] = 1.0;
  for (i=1 ; i<=p ; ++i)
    retval[1] -= as[i];
  for (i=0 ; i<q ; ++i)
    retval[1] -= bs[i];

  return retval;
}

Vector GARCHModel::GetMask()
{
  int i;
  Vector retval(p+q+2);
  for (i=0 ; i<min(estimation_mask.size(),p+q+2) ; ++i)
    retval[i] = estimation_mask[i];
  for ( ; i<p+q+2 ; ++i)
    retval[i] = 1;
  return retval;
}

Vector GARCHModel::GetDefaultParameters(Vector& mask, TimeSeries *tsp)
// can modify mask if necessary, for instance so as not to use search for best mean
{
  int i,j,pstar=max(p,q),nobs=tsp->GetN();
  Vector parms = ParameterBundle(), acf(1);

  // mean is handled differently
  if (mask[0]) // mean is to be estimated
    {
      mask[0] = 0;
      parms[0] = tsp->SampleMean();
    }

  // initialize all non-masked parameters now
  tsp->ComputeSampleACF(&acf,NULL,0);

  // first find a reasonable starting point
  Vector localas(p+1),localbs(q+1);
  if (p+q>0)
    {
      // set up a simple regression
      Matrix X(nobs-pstar,pstar+1),XtX;
      Vector y(nobs-pstar);
      for (i=pstar ; i<nobs ; ++i)
	{
	  y[i-pstar] = (*tsp)[i]*(*tsp)[i];
	  for (j=0 ; j<pstar ; ++j)
	    X[i-pstar][j+1] = (*tsp)[i-j-1]*(*tsp)[i-j-1];
	  X[i-pstar][0] = 1.0;
	}
      y = X.transpose()*y;
      XtX = X.transpose()*X;
      XtX.symm_solve_system(&y);
      
      if (y[0]<zero)
	y[0] = zero;
      for (i=1 ; i<y.nrows() ; ++i)
	if (y[i]<0)
	  y[i] = 0.0;
      localas[0] = y[0];
      for (i=0 ; i<p ; ++i)
	if (i<q)
	  localas[i+1] = y[i+1]/2.0;
	else
	  localas[i+1] = y[i+1];
      for (i=0 ; i<q ; ++i)
	if (i<p)
	  localbs[i] = y[i+1]/2.0;
	else
	  localbs[i] = y[i+1];
    }
  else
    localas[0] = acf[0];

  // copy all free parameters into return vector
  for (i=0 ; i<=p ; ++i)
    if (mask[i+1])
      parms[i+1] = localas[i];
  for (i=0 ; i<q ; ++i)
    if (mask[i+p+2])
      parms[i+p+2] = localbs[i];

  return parms;
}

Vector GARCHModel::ComputeSpectralDensity(Vector& omegas)
// GARCH is special case of WN->hence has flat spectrum
{
  int i,j;
  Vector retval(omegas.size()),acf(1);
  complex tc1,tc2,eil,tsum,thetapart,phipart;

  ComputeACFPACF(acf,NULL,false);  // get the variance
  for (i=0 ; i<omegas.size() ; ++i)
    retval[i] = acf[0]/(2*M_PI);

  return retval;
}

//---------------------------------------------------+
//       EGARCH stuff follows                        |
//---------------------------------------------------+
EGARCHModel::EGARCHModel()
  : GARCHModel()
{
  gs.resize(100);
  gs.zeroes();
}

EGARCHModel::~EGARCHModel()
{
}

int EGARCHModel::FitModel(TimeSeries *ts, const int method, const int ilimit,
			 void (*itercallback)(int,void *), void *cb_parms,
			 ostringstream& msg, ostringstream& supplemental,
			 bool get_parameter_cov)
{
  int i,j,nobs = ts->GetN(),pstar=max(p,q);
  Vector localas(p+1),localbs(q);
  Vector acf(1),mask;
  double mu = ts->SampleMean();
  Matrix parmcovs;

  iterationlimit = ilimit;
  mask = GetMask();  // get current mask (pad if necessary with default values)

  if (method==METHOD_MLE)
    {
      if (ts->HasAnyMissing())
	{
	  msg << "Unable to fit EGARCH model when time series has missing observations."
	      << ends;
	  return UNABLE;
	}
      
      parmcovs = MLEFit(ts, mask, itercallback, cb_parms, get_parameter_cov);
      StreamMLESummary(supplemental, parmcovs, get_parameter_cov);

      msg << "EGARCH model fit in " << iterationcounter << " iterations." << ends;

      CheckModel(false);
      return SUCCESS;
    }
  if (method==METHOD_BAYESIAN)
    {
      msg << "MCMC estimation of EGARCH models is not yet supported." << ends;
      return UNABLE;
    }
}

void EGARCHModel::RenderInto(ostringstream& os)
{
  int i;
  os << "\\b";  // make it bold
  os << "EGARCH(" << p << "," << q << ") Model" << endl;
  os << "\\n";  // make it normal
  os << "X(t) = \\gs\\n(t) Z(t)" << endl;
  os << "log(\\gs\\n(t)\\s2\\n) = ";
  os << setiosflags(ios::fixed) << setprecision(3) << as[0];
  for (i=1 ; i<=p ; ++i)
    os << " + " << setiosflags(ios::fixed) << setprecision(3) << as[i] << "[|Z(t-" << i << ")| + "
       << gs[i] << "Z(t-" << i << ")]";
  for (i=1 ; i<=q ; ++i)
    os << " + " << setiosflags(ios::fixed) << setprecision(3) << bs[i-1] << "log(\\gs(t-" << i << ")\\s2\\n)";
  os << endl << "Mean = " << setprecision(5) << mean << ",  " << endl;
  if (loglikelihood_is_valid)
    {
      os << endl << endl << "Log Likelihood = " << GetLogLikelihood() << endl;
      os << "AICC = " << GetAICC() << endl;
    }
  else
    {
      os << endl << endl << "Log Likelihood = <Not Computed>" << endl;
      os << "AICC = <Not Computed>" << endl;
    }
}

double EGARCHModel::OneStepVariance(double *xt, double *st, int t, double gamma0)
{
  double sum=0, zt, tx;
  int j,m = max(p,q);
  if (t+1>=m)
    {      
      sum=as[0];
      for (j=1 ; j<=p ; ++j)
	{
	  zt = xt[1-j]/sqrt(st[1-j]);
	  sum += as[j]*(fabs(zt) + gs[j]*zt);
	}
      for (j=0 ; j<q ; ++j)
	sum += bs[j]*log(st[-j]);
    }
  else
    {
      // find E(log(sigma_t^2)) first
      sum = as[0];
      for (j=1 ; j<=p ; ++j)
	sum += 0.797*as[j];
      for (j=1,tx=1.0 ; j<=q ; ++j)
	tx -= bs[j-1];
      tx = sum/tx;  // tx is expected value

      sum=as[0];
      for (j=1 ; j<=p ; ++j)
	{
	  if (1-j>=-t)
	    {
	      zt = xt[1-j]/sqrt(st[1-j]);
	      sum += as[j]*(fabs(zt) + gs[j]*zt);
	    }
	  else
	    sum += as[j]*0.797885;
	}
      for (j=0 ; j<q ; ++j)
	if (-j>=-t)
	  sum += bs[j]*log(st[-j]);
	else
	  sum += bs[j]*tx;
    }
  // simple thresholding to prevent numerical over/under-flows
  if (sum>20)
    return exp(20);
  else
    if (sum<-200)
      return exp(-200);
    else
      return exp(sum);
}

void EGARCHModel::ComputeACFPACF(Vector& acf, Vector* pacf, bool normalize)
{
  int i, j;
  double tx, sum;

  acf.zeroes();
  if (pacf!=NULL)
    {
      (*pacf)[0] = 1.0;
      pacf->zeroes();
    }
  if (normalize)
    acf[0] = 1.0;
  else
    {
      // find E(log(sigma_t^2)) first
      sum = as[0];
      for (j=1 ; j<=p ; ++j)
	sum += 0.797*as[j];
      for (j=1,tx=1.0 ; j<=q ; ++j)
	tx -= bs[j-1];
      tx = sum/tx;  // tx is expected value

      acf[0] = exp(tx); // incorrect!
    }
}

void EGARCHModel::StabilizeBs()
{
  int i,q2=0;
  double min, tx, ty;
  complex *roots;
  double *cp;

  // first construct moving average polynomial
  for (i=0 ; i<q ; ++i)
    if (fabs(bs[i])>zero)
      q2 = i+1;					// p2 is the effective MA order
  if (q2==0)
    return;

  cp = new double[q+1];
  cp[0] = 1.0;
  for (i=1 ; i<=q ; ++i)
    cp[i] = -bs[i-1];

  roots = rootsof(q2, cp);

  // now find minimum magnitude
  min = abs(roots[0]);
  for (i=1 ; i<=q2 ; ++i)
    if (abs(roots[i-1])<min)
      min = abs(roots[i-1]);
  delete[] cp;
  if (min>1+root_margin)
    {
      delete[] roots;
      return;
    }

  // otherwise we force it to have roots outside the unit circle
  complex *recon;
  for (i=1 ; i<=q2 ; ++i)
    {
      tx = abs(roots[i-1]);
      ty = arg(roots[i-1]);
      if (tx<1+root_margin)
	tx = 1+root_margin;
      roots[i-1] = polar(tx,ty);
    }
  recon = coeffsof(q2, roots);
  for (i=0 ; i<q2 ; ++i)
    bs[i] = -real(recon[i+1]/recon[0]);

  delete[] recon;
  delete[] roots;
}

bool EGARCHModel::CheckModel(bool flipit)
{
  bool retval = true;
  Vector tv = ValidityConstraint();
  if (!flipit)
    {
      retval = (tv[0] >= 0); 
    }
  else // actually fix it
    {
      if (tv[0] < 0) // then it's not OK
	StabilizeBs();
      retval = true;
    }
  return retval;
}

Vector EGARCHModel::ValidityConstraint()
{
  Vector retval(1);

  // first constraint: beta_i roots outside unit circle
  int i, q2=0;
  double min, tx;
  double *cp;
  complex *roots;

  // first construct polynomial from b coefficients
  for (i=0 ; i<q ; ++i)
    if (fabs(bs[i])>zero)
      q2 = i+1;					// q2 is the effective order, is < q if bs_q=0
  if (q2==0)
    min = 1e6;                                  // a bad approximation to \infty
  else
    {
      cp = new double[q2+1];
      cp[0] = 1.0;
      for (i=1 ; i<=q2 ; ++i)
	cp[i] = -bs[i-1];
      roots = rootsof(q2, cp);
      
      // now find minimum magnitude
      min = abs(roots[0]);
      for (i=1 ; i<=q2 ; ++i)
	{
	  tx = abs(roots[i-1]);
	  if (tx<min)
	    min = tx;
	}
      delete[] cp;
      delete[] roots;
    }

  retval[0]=min-1-root_margin;
  return retval;
}

Vector EGARCHModel::ParameterBundle()
{
  int i;
  Vector pvec(2*p+q+2);
  pvec[0] = mean;
  pvec[1] = as[0];
  for (i=1 ; i<=p ; ++i)
    {
      pvec[i+1] = as[i];
      pvec[i+p+1] = gs[i];
    }
  for (i=0 ; i<q ; ++i)
    pvec[i+2*p+2] = bs[i];
  return pvec;
}

void EGARCHModel::UnbundleParameters(Vector &pvec)
{
  int i;
  mean = pvec[0];
  as[0] = pvec[1];
  for (i=1 ; i<=p ; ++i)
    {
      as[i] = pvec[i+1];
      gs[i] = pvec[i+p+1];
    }
  for (i=0 ; i<q ; ++i)
    bs[i] = pvec[i+2*p+2];
}

void EGARCHModel::StreamParameterName(ostringstream& os, int pnum)
{
  if (pnum==0)
    os << "\u03bc = mean"; 
  else if (pnum<=(p+1))
    os << "\u03b1(" << (pnum-1) << ")";
  else if (pnum<=(2*p+1))
    os << "\u03b6(" << (pnum-1-p) << ")";
  else
    os << "\u03b2(" << (pnum-2*p-1) << ")";
}

  
void EGARCHModel::SetParameters(int p0, int q0, Vector& as0, Vector& gs0, Vector& bs0, double m)
{
  int i;
  p=p0;   q=q0;
  as.resize(p+1);
  gs.resize(p+1);
  bs.resize(q);
  as[0] = as0[0];
  for (i=1 ; i<=p ; ++i)
    {
      as[i] = as0[i];
      gs[i] = gs0[i];
    }
  for (i=0 ; i<q ; ++i)
    bs[i] = bs0[i];
  mean=m;
}

void EGARCHModel::GetParameters(int& p1, int& q1, Vector& as1, Vector& gs1, Vector& bs1, double& mu)
{
  p1=p;   q1=q;
  as1=as;   bs1=bs;   gs1=gs;
  mu=mean;
}


Vector EGARCHModel::GetMask()
// returns estimation_mask, padded with defaults
{
  int i;
  Vector retval(2*p+q+2);
  for (i=0 ; i<min(estimation_mask.size(),2*p+q+2) ; ++i)
    retval[i] = estimation_mask[i];
  for ( ; i<2*p+q+2 ; ++i)
    retval[i] = 1;
  return retval;
}

Vector EGARCHModel::GetDefaultParameters(Vector& mask, TimeSeries *tsp)
// returns initial guesses, in form of bundled parameter vector
{
  int i,j,pstar=max(p,q),nobs=tsp->GetN();
  Vector parms = ParameterBundle(), acf(1);
  double tx;

  // mean is handled differently
  if (mask[0]) // mean is to be estimated
    {
      mask[0] = 0;
      parms[0] = tsp->SampleMean();
    }

  // initialize all non-masked parameters now
  tsp->ComputeSampleACF(&acf,NULL,0);

  // first find a reasonable starting point
  Vector localas(p+1),localbs(q+1),localgs(p+1);
  Vector zts(nobs);
  for (i=0 ; i<nobs ; ++i)
    zts[i] = (*tsp)[i];
  tx = zts.var();
  zts = zts*(1.0/sqrt(tx));
  if (p+q>0)
    {
      // set up a simple regression
      Matrix X(nobs-pstar,2*p+q+1),XtX;
      Vector y(nobs-pstar);
      for (i=pstar ; i<nobs ; ++i)
	{
	  y[i-pstar] = log((*tsp)[i]*(*tsp)[i]+1e-6);    // log(X_t^2+1e-6) is a proxy for log(sigma^2)
	  X[i-pstar][0] = 1.0;
	  for (j=0 ; j<p ; ++j)
	    {
	      X[i-pstar][j+1] = fabs(zts[i-j-1]);
	      X[i-pstar][j+p+1] = zts[i-j-1];
	    }
	  for (j=0 ; j<q ; ++j)
	    X[i-pstar][j+2*p+1] = log((*tsp)[i-j-1]*(*tsp)[i-j-1]+1e-6);
	}
      y = X.transpose()*y;
      XtX = X.transpose()*X;
      XtX.symm_solve_system(&y);
      
      localas[0] = y[0];
      for (i=0 ; i<p ; ++i)
	{
	  localas[i+1] = y[i+1];
	  localgs[i+1] = y[i+1+p];
	}
      for (i=0 ; i<q ; ++i)
	localbs[i] = y[i+2*p+1];
    }
  else
    localas[0] = log(acf[0]);

  // copy all free parameters into return vector
  // (note: we already did the mean above)
  for (i=0 ; i<=p ; ++i)
    {
      if (mask[i+1])
	parms[i+1] = localas[i];
      if (mask[i+p+1])
	parms[i+p+1] = localgs[i+1];
    }
  for (i=0 ; i<q ; ++i)
    if (mask[i+2*p+2])
      parms[i+2*p+2] = localbs[i];

  return parms;
}
